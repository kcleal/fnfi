"""
A basic assembler. Takes an overlap graph and merges reads in-place in a pileup style. Different soft-clipped regions
are then overlapped of 'linked'.
"""

import networkx as nx
import difflib
import numpy as np
from collections import defaultdict
import click


def echo(*args):
    click.echo(args, err=True)


def to_prob(q):
    return pow(10, (-1*q/10))


def to_phred(p):
    return int(-10*np.log10(p + 1e-9))  # Ovoid overflow error


def update_edge(u, v, qual, G, kind, strand):
    if G.has_node(v):
        G.node[v]["w"] = to_phred(to_prob(G.node[v]["w"]) * to_prob(qual))  # A and B
        G.node[v]["n"] += 1
        G.node[v]["strand"] += strand
        if not G.has_edge(u, v):
            G.add_edge(u, v)
    else:
        G.add_node(v, w=qual, kind=kind, n=1, strand=strand, rid=[])
        G.add_edge(u, v)


def get_ref_pos(cigar, pos, seq):
    # pysam get_reference_positions(full_length=True) does'nt work for last aligner-styled cigars, hence rolled my own
    p = []  # A list of positions
    o = []  # A list of offsets from last seen position (only insertions or soft-clips have offsets)
    t = []  # The block type, insertion = i, left clip = l, right clip = r, else None
    current_pos = pos + 1
    start = True
    for opp, length in cigar:
        if opp == 4:
            p += [current_pos] * length

            if start:
                t += ["l"] * length
                o += range(length, 0, -1)
            else:
                t += ["r"] * length
                o += range(length)

        elif opp == 0 or opp == 7 or opp == 8 or opp == 3:  # All match, match (=), mis-match (X), N's
            # N's are still in reference
            p += [j for j in range(current_pos, current_pos + length)]
            o += [0] * length
            t += [None] * length
            current_pos += length

        elif opp == 1:  # Insertion
            p += [current_pos] * length
            o += range(1, length + 1)
            t += ["i"] * length
            current_pos += 1  # Advance to the next alignment block

        elif opp == 5 or opp == 2:  # Hard clip or deletion
            current_pos += length

        start = False

    return zip(seq, p, o, t)


def find_cycle(G, source=None, orientation=None):
    """Stolen from networkX 2
    Returns a cycle found via depth-first traversal.

    The cycle is a list of edges indicating the cyclic path.
    Orientation of directed edges is controlled by `orientation`.

    Parameters
    ----------
    G : graph
        A directed/undirected graph/multigraph.

    source : node, list of nodes
        The node from which the traversal begins. If None, then a source
        is chosen arbitrarily and repeatedly until all edges from each node in
        the graph are searched.

    orientation : None | 'original' | 'reverse' | 'ignore' (default: None)
        For directed graphs and directed multigraphs, edge traversals need not
        respect the original orientation of the edges.
        When set to 'reverse' every edge is traversed in the reverse direction.
        When set to 'ignore', every edge is treated as undirected.
        When set to 'original', every edge is treated as directed.
        In all three cases, the yielded edge tuples add a last entry to
        indicate the direction in which that edge was traversed.
        If orientation is None, the yielded edge has no direction indicated.
        The direction is respected, but not reported.

    Returns
    -------
    edges : directed edges
        A list of directed edges indicating the path taken for the loop.
        If no cycle is found, then an exception is raised.
        For graphs, an edge is of the form `(u, v)` where `u` and `v`
        are the tail and head of the edge as determined by the traversal.
        For multigraphs, an edge is of the form `(u, v, key)`, where `key` is
        the key of the edge. When the graph is directed, then `u` and `v`
        are always in the order of the actual directed edge.
        If orientation is not None then the edge tuple is extended to include
        the direction of traversal ('forward' or 'reverse') on that edge.

    Raises
    ------
    NetworkXNoCycle
        If no cycle was found.

    Examples
    --------
    In this example, we construct a DAG and find, in the first call, that there
    are no directed cycles, and so an exception is raised. In the second call,
    we ignore edge orientations and find that there is an undirected cycle.
    Note that the second call finds a directed cycle while effectively
    traversing an undirected graph, and so, we found an "undirected cycle".
    This means that this DAG structure does not form a directed tree (which
    is also known as a polytree).

    >>> import networkx as nx
    >>> G = nx.DiGraph([(0, 1), (0, 2), (1, 2)])
    >>> try:
    ...    nx.find_cycle(G, orientation='original')
    ... except:
    ...    pass
    ...
    >>> list(nx.find_cycle(G, orientation='ignore'))
    [(0, 1, 'forward'), (1, 2, 'forward'), (0, 2, 'reverse')]

    """
    if not G.is_directed() or orientation in (None, 'original'):
        def tailhead(edge):
            return edge[:2]
    elif orientation == 'reverse':
        def tailhead(edge):
            return edge[1], edge[0]
    elif orientation == 'ignore':
        def tailhead(edge):
            if edge[-1] == 'reverse':
                return edge[1], edge[0]
            return edge[:2]

    explored = set()
    cycle = []
    final_node = None
    for start_node in G.nbunch_iter(source):
        if start_node in explored:
            # No loop is possible.
            continue

        edges = []
        # All nodes seen in this iteration of edge_dfs
        seen = {start_node}
        # Nodes in active path.
        active_nodes = {start_node}
        previous_head = None

        for edge in nx.edge_dfs(G, start_node, orientation):
            # Determine if this edge is a continuation of the active path.
            tail, head = tailhead(edge)
            if head in explored:
                # Then we've already explored it. No loop is possible.
                continue
            if previous_head is not None and tail != previous_head:
                # This edge results from backtracking.
                # Pop until we get a node whose head equals the current tail.
                # So for example, we might have:
                #  (0, 1), (1, 2), (2, 3), (1, 4)
                # which must become:
                #  (0, 1), (1, 4)
                while True:
                    try:
                        popped_edge = edges.pop()
                    except IndexError:
                        edges = []
                        active_nodes = {tail}
                        break
                    else:
                        popped_head = tailhead(popped_edge)[1]
                        active_nodes.remove(popped_head)

                    if edges:
                        last_head = tailhead(edges[-1])[1]
                        if tail == last_head:
                            break
            edges.append(edge)

            if head in active_nodes:
                # We have a loop!
                cycle.extend(edges)
                final_node = head
                break
            else:
                seen.add(head)
                active_nodes.add(head)
                previous_head = head

        if cycle:
            break
        else:
            explored.update(seen)

    else:
        assert(len(cycle) == 0)
        raise nx.exception.NetworkXNoCycle('No cycle found.')

    # We now have a list of edges which ends on a cycle.
    # So we need to remove from the beginning edges that are not relevant.

    for i, edge in enumerate(cycle):
        tail, head = tailhead(edge)
        if tail == final_node:
            break

    return cycle[i:]


def dag_longest_path(G, weight='weight', default_weight=1):
    """Returns the longest path in a directed acyclic graph (DAG).

    If `G` has edges with `weight` attribute the edge data are used as
    weight values.

    Parameters
    ----------
    G : NetworkX DiGraph
        A directed acyclic graph (DAG)

    weight : str, optional
        Edge data key to use for weight

    default_weight : int, optional
        The weight of edges that do not have a weight attribute

    Returns
    -------
    list
        Longest path

    Raises
    ------
    NetworkXNotImplemented
        If `G` is not directed

    See also
    --------
    dag_longest_path_length

    """
    if not G:
        return []
    dist = {}  # stores {v : (length, u)}
    for v in nx.topological_sort(G):
        us = [(dist[u][0] + data.get(weight, default_weight), u)
              for u, data in G.pred[v].items()]
        # Use the best predecessor if there is one and its distance is
        # non-negative, otherwise terminate.
        maxu = max(us, key=lambda x: x[0]) if us else (0, v)
        dist[v] = maxu if maxu[0] >= 0 else (0, v)
    u = None
    v = max(dist, key=lambda x: dist[x][0])
    path = []
    while u != v:
        path.append(v)
        u = v
        v = dist[v][1]
    path.reverse()
    return path


def base_assemble(g, reads, bam, id=0):
    """
    Assembles reads that have overlaps. Uses alignment positions to determine contig construction
    :param g: The overlap graph
    :param reads: Dict of read_name: flag: alignment
    :param bam: Original bam for header access
    :param id: Unique ID for the event
    :return: Returns None if no soft-clipped portion of the cluster was assembled, otherwise a result dict is returned
    """
    # Note supplementary are included in assembly; helps link regions
    # Get reads of interest

    rd = [reads[n[0]][(n[1], n[2])] for n in g.nodes()]

    # rnames = set([r.qname for r in rd])
    # roi = "simulated_reads.3.10-id293_A_chr21:46699688_B_chr1:38378863-38649"
    G = nx.DiGraph()

    strand_d = {}
    start_end_rids = defaultdict(list)
    for r in rd:
        if r.seq is None:
            continue

        # (base, position, offset)
        seq_pos = iter(zip(get_ref_pos(r.cigartuples, r.pos, r.seq), r.query_qualities))

        rid = (r.qname, r.flag, r.pos)  # Read id
        u, qual_u = next(seq_pos)
        first_node = u
        for v, qual_v in seq_pos:

            if G.has_edge(u, v):
                G[u][v]["weight"] += int(qual_u + qual_v)
            else:
                G.add_edge(u, v, weight=int(qual_u + qual_v))
            u = v
            qual_u = qual_v

        # Add read id to first and last in path
        start_end_rids[first_node].append(rid)
        start_end_rids[v].append(rid)
        strand_d[rid] = -1 if r.flag & 16 else 1

    try:
        path = dag_longest_path(G, weight="weight")
    except:
        cy = find_cycle(G, orientation="original")
        echo("Cycle in graph", cy)
        for r in rd:
            echo(bam.get_reference_name(r.rname))
            echo(str(r).split("\t"))
        quit()

    longest_left_sc = path[0][2]
    longest_right_sc = path[-1][2]

    if longest_left_sc == 0 and longest_right_sc == 0:
        return None  # No soft-clips, so not overlapping a break

    bases = "".join(i.upper() if k == 0 else i.lower() for i, j, k, l in path)
    read_names = []
    strand_counts = []
    for item in path:
        read_names += start_end_rids[item]
        if item in strand_d:
            strand_counts.append(strand_d[item])
            del strand_d[item]  # only count once

    matches = [i[1] for i in path if i[2] == 0]

    res = {"bamrname": bam.get_reference_name(rd[0].rname),
           "left_clips": longest_left_sc,
           "right_clips": longest_right_sc,
           "strand_l": strand_counts.count(-1),
           "strand_r": strand_counts.count(1),
           "ref_start": matches[0],
           "ref_end": matches[-1] + 1,
           "read_names": set(read_names),
           "contig": bases,
           "id": id}

    return res


def rev_comp(s):
    d = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N", "|": "|", "a": "t", "t": "a", "c": "g", "g": "c",
         "n": "n"}
    return "".join(d[j] for j in s if j != "|")[::-1]


def explore_local(starting_nodes, large_component, color, upper_bound):
    seen = set(starting_nodes)
    found = set([])
    if len(starting_nodes) == 0:
        return set([])
    while True:
        nd = starting_nodes.pop()
        seen.add(nd)
        for edge in large_component.edges(nd, data=True):
            if edge[2]['c'] == color:
                if edge[0] not in seen:
                    starting_nodes.add(edge[0])
                    found.add(edge[0])
                elif edge[1] not in seen:
                    starting_nodes.add(edge[1])
                    found.add(edge[1])
            if len(found) > upper_bound:
                return set([])
        if len(starting_nodes) == 0:
            break
    return found


def link_pair_of_assemblies(a, b, clip_length):
    seqs = []
    for i in (a, b):
        left_clipped = False
        negative_strand = True if i["strand_r"] < 0 else False
        sc_support = i["right_clips"]

        if i["left_clips"] > i["right_clips"]:
            left_clipped = True
            negative_strand = True if i["strand_l"] < 0 else False
            sc_support = i["left_clips"]

        seqs.append((i["contig"], sc_support, i["bamrname"], i["ref_start"], i["ref_end"] + 1, left_clipped,
                     negative_strand))

    ainfo, binfo = seqs
    aseq, bseq = ainfo[0], binfo[0]
    if ainfo[2] == binfo[2]:  # If seqs are on the same chrom
        if ainfo[5] == binfo[5]:  # If clips are on same side, rev comp one of them
            bseq = rev_comp(bseq)
    else:
        if ainfo[6]:  # negative strand
            aseq = rev_comp(aseq)
        if binfo[6]:
            bseq = rev_comp(bseq)

    # See https://docs.python.org/2/library/difflib.html
    m = difflib.SequenceMatcher(a=aseq.upper(), b=bseq.upper(), autojunk=None)
    longest = m.find_longest_match(0, len(aseq), 0, len(bseq))

    a_align = [i.islower() for i in aseq[longest[0]:longest[0] + longest[2]]]
    b_align = [i.islower() for i in bseq[longest[1]:longest[1] + longest[2]]]

    sc_a = sum(a_align)
    sc_b = sum(b_align)
    # non_sc_a = len(a_align) - sc_a
    # non_sc_b = len(b_align) - sc_b

    best_sc = max([sc_a, sc_b])
    # best_non_sc = max([non_sc_a, non_sc_b])
    a["linked"] = "weak link"
    a["best_sc"] = best_sc
    if best_sc > clip_length:  # and best_non_sc >= 5:
        a["linked"] = "strong link"

    return a, b
