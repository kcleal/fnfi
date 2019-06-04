# Networkx functions, from v1/2. Saved here for compatibility
from __future__ import absolute_import
import networkx as nx
import click
from . import c_cluster_funcs
import itertools
from . import data_io
from collections import defaultdict, deque
import time


class Alignment(object):
    """Picklable struct to hold the contents of pysam alignment"""
    __slots__ = ["reference_end", "cigar", "pos", "flag", "rname", "qname", "rnext", "pnext", "seq", "cigartuples",
                 "query_qualities", "has_SA"]

    def __init__(self, a):
        self.reference_end = a.reference_end
        self.cigar = a.cigar
        self.pos = a.pos
        self.flag = a.flag
        self.rname = a.rname
        self.qname = a.qname
        self.rnext = a.rnext
        self.pnext = a.pnext
        self.seq = a.seq
        self.cigartuples = self.cigar
        self.query_qualities = a.query_qualities
        self.has_SA = a.has_tag("SA")


class Scoper(object):
    """Keeps track of which reads are in scope. Maximum distance depends on the template insert_median"""
    def __init__(self, max_dist, clip_length=20, max_depth=60):  # Todo add parameter for this
        self.max_dist = max_dist
        self.max_depth = max_depth
        self.scope = deque([])
        self.clip_length = clip_length

    def update(self, input_read):

        current_read = Alignment(input_read)  # Convert to friendly format
        if len(self.scope) == 0:
            self.scope.append(current_read)
            return current_read

        elif len(self.scope) > 0 and current_read.rname != self.scope[-1].rname:
            self.scope = deque([current_read])  # Empty scope on new chromosome
            return current_read

        current_pos = current_read.pos

        while True:
            if len(self.scope) == 0:
                break
            if len(self.scope) > self.max_depth:
                self.scope.popleft()
            else:
                break

        while True:
            if len(self.scope) == 0:
                break
            if abs(self.scope[0].pos - current_pos) > self.max_dist:
                self.scope.popleft()
            else:
                break

        self.scope.append(current_read)
        return current_read

    @staticmethod
    def overlap(start1, end1, start2, end2):
        return max(0, min(end1, end2) - max(start1, start2))

    def iterate(self):
        if len(self.scope) > 1:

            for i in range(len(self.scope)):
                if i == len(self.scope) - 1:
                    break  # Don't pair with self
                t = self.scope[i]
                if t.qname == self.scope[-1].qname:
                    continue  # Don't pair with same template
                current = self.scope[-1]

                # If current is a supplementary, it MAY be hard-clipped so only look for alignment-overlaps
                if current.flag & 2048 and (current.cigar[0][0] == 5 or current.cigar[-1][0] == 5):

                    # Make sure hard-clip can be overlapped with soft-clip on t
                    if current.cigar[0][0] == 5 and current.cigar[0][1] >= self.clip_length:  # Left HARD-clip current
                        if t.cigar[0][0] == 4 and t.cigar[0][1] >= self.clip_length:  # Left soft-clip on t
                            # Extent of hard and soft-clip overlap
                            c_sc_start = current.pos - current.cigar[0][1]
                            t_sc_start = t.pos - t.cigar[0][1]
                            if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                                # Read, soft-clip overlap primary edge, hard-clip overlap supplementary edge
                                yield t, False, True

                    elif current.cigar[-1][0] == 5 and current.cigar[-1][1] > self.clip_length:  # Right HARD-clip
                        if t.cigar[-1][0] == 4 and t.cigar[-1][1] > self.clip_length:  # Right soft-clip
                            c_sc_end = current.reference_end + current.cigar[-1][1]
                            t_sc_end = t.reference_end + t.cigar[-1][1]
                            if self.overlap(t.reference_end, t_sc_end, current.reference_end,
                                            c_sc_end) >= self.clip_length:
                                yield t, False, True

                #
                if current.cigar[0][0] == 4 and current.cigar[0][1] >= self.clip_length:  # Left soft-clip on current
                    if t.cigar[0][0] == 4:  # and t.cigar[0][1] >= self.clip_length:  # Left soft-clip on t
                        # Extent of sof-clip, make sure-soft-clip portion overlaps not the alignment part
                        c_sc_start = current.pos - current.cigar[0][1]
                        t_sc_start = t.pos - t.cigar[0][1]
                        if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                            yield t, True, False

                    elif t.cigar[0][0] == 5 and t.flag & 2048 and t.cigar[0][1] >= self.clip_length:  # Hard clip on t
                        c_sc_start = current.pos - current.cigar[0][1]
                        t_sc_start = t.pos - t.cigar[0][1]
                        if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                            yield t, False, True

                elif current.cigar[0][0] == 4:
                    if t.cigar[0][0] == 4 and t.cigar[0][1] >= self.clip_length:  # Left soft-clip on t
                        # Extent of sof-clip, make sure-soft-clip portion overlaps not the alignment part
                        c_sc_start = current.pos - current.cigar[0][1]
                        t_sc_start = t.pos - t.cigar[0][1]
                        if self.overlap(t_sc_start, t.pos, c_sc_start, current.pos) >= self.clip_length:
                            yield t, True, False

                elif current.cigar[-1][0] == 4 and current.cigar[-1][1] > self.clip_length:  # Right soft-clip
                    if t.cigar[-1][0] == 4:  # and t.cigar[-1][1] > self.clip_length:
                        c_sc_end = current.reference_end + current.cigar[-1][1]
                        t_sc_end = t.reference_end + t.cigar[-1][1]
                        if self.overlap(t.reference_end, t_sc_end, current.reference_end, c_sc_end) >= self.clip_length:
                            yield t, True, False

                    elif t.cigar[-1][0] == 5 and t.flag & 2048 and t.cigar[-1][1] >= self.clip_length:
                        c_sc_end = current.reference_end + current.cigar[-1][1]
                        t_sc_end = t.reference_end + t.cigar[-1][1]
                        if self.overlap(t.reference_end, t_sc_end, current.reference_end, c_sc_end) >= self.clip_length:
                            yield t, False, True

                elif current.cigar[-1][0] == 4:
                    if t.cigar[-1][0] == 4 and t.cigar[-1][1] > self.clip_length:
                        c_sc_end = current.reference_end + current.cigar[-1][1]
                        t_sc_end = t.reference_end + t.cigar[-1][1]
                        if self.overlap(t.reference_end, t_sc_end, current.reference_end, c_sc_end) >= self.clip_length:
                            yield t, True, False

                # else:
                yield t, False, False

        yield None, None, None  # When not enough reads are in scope


def echo(*args):
    click.echo(args, err=True)


def blockmodel(G, partitions):
    """Stolen from NetworkX 1.9. Returns a reduced graph constructed using the generalized block modeling
    technique. Faster than the quotient graph method, but not by much.

    The blockmodel technique collapses nodes into blocks based on a
    given partitioning of the node set.  Each partition of nodes
    (block) is represented as a single node in the reduced graph.

    Edges between nodes in the block graph are added according to the
    edges in the original graph.  If the parameter multigraph is False
    (the default) a single edge is added with a weight equal to the
    sum of the edge weights between nodes in the original graph
    The default is a weight of 1 if weights are not specified.  If the
    parameter multigraph is True then multiple edges are added each
    with the edge data from the original graph.

    Parameters
    ----------
    G : graph
        A networkx Graph or DiGraph
    partitions : list of lists, or list of sets
        The partition of the nodes.  Must be non-overlapping.
    multigraph : bool, optional
        If True return a MultiGraph with the edge data of the original
        graph applied to each corresponding edge in the new graph.
        If False return a Graph with the sum of the edge weights, or a
        count of the edges if the original graph is unweighted.

    Returns
    -------
    blockmodel : a Networkx graph object

    Examples
    --------
    >>> G=nx.path_graph(6)
    >>> partition=[[0,1],[2,3],[4,5]]
    >>> M=nx.blockmodel(G,partition)

    References
    ----------
    .. [1] Patrick Doreian, Vladimir Batagelj, and Anuska Ferligoj
    	"Generalized Blockmodeling",Cambridge University Press, 2004.
    """
    # Create sets of node partitions, frozenset makes them hashable
    part = [frozenset(i) for i in partitions]

    # Check for overlapping node partitions
    u = set()
    for p1, p2 in zip(part[:-1], part[1:]):
        u.update(p1)
        if len(u.intersection(p2)) > 0:
            raise nx.NetworkXException("Overlapping node partitions.")

    # Initialize blockmodel graph
    M = nx.Graph()

    # Add nodes and properties to blockmodel
    # The blockmodel nodes are node-induced subgraphs of G
    # Label them with integers starting at 0
    for i, p in zip(range(len(part)), part):
        M.add_node(p)
        # The node-induced subgraph is stored as the node 'graph' attribute
        SG = G.subgraph(p)
        M.node[p]['graph'] = SG

    # Create mapping between original node labels and new blockmodel node labels
    block_mapping = {}
    for n in M:
        nodes_in_block = M.node[n]['graph'].nodes()
        block_mapping.update(dict.fromkeys(nodes_in_block, n))

    # Add edges to block graph
    for u, v, d in G.edges(data=True):

        bmu = block_mapping[u]
        bmv = block_mapping[v]
        if bmu == bmv:  # no self loops
            continue

        # For graphs and digraphs add single weighted edge
        weight = d.get('weight', 1.0)  # default to 1 if no weight specified
        if M.has_edge(bmu, bmv):
            M[bmu][bmv]['weight'] += weight
        else:
            M.add_edge(bmu, bmv, weight=weight)
    return M


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


def skip_traingles(G, u, v):
    # Dont induce triangles. Prevents most redundant edges being formed
    if not G.has_node(u) or not G.has_node(v):
        return False

    for out_u in G[u]:
        for neighbor_u in G[out_u]:
            if neighbor_u == v:
                return G[out_u][neighbor_u]["c"]


def construct_graph(infile, max_dist, tree, buf_size=100000, input_windows=()):
    t0 = time.time()
    click.echo("Building graph from {} regions".format(len(input_windows)), err=True)
    all_flags = defaultdict(set)  # Linking clusters together rname: set([(flag, pos)..])

    # Make a buffer of reads to help prevent file lookups later
    read_buffer = dict()  # keys are (rname, flag, pos)
    read_index_buffer = dict()  # keys are int, values are (rname, flag, pos)

    scope = Scoper(max_dist=max_dist)

    count = 0
    buf_del_index = 0

    grey_added = 0
    black_edges = 0

    ii = 0
    G = nx.Graph()
    for widx, window in enumerate(input_windows):

        for r in infile.fetch(*window):

            if len(r.cigar) == 0 or r.flag & 1028:  # Unmapped or duplicate

                continue

            n1 = (r.qname, r.flag, r.pos)

            # if n1 not in name_to_node:
            #     name_to_node[n1] = ii
            ii += 1

            # Add read to buffer
            read_buffer[n1] = r
            read_index_buffer[count] = n1
            count += 1

            # Reduce reads in buffer if too many
            if len(read_buffer) > buf_size:
                if buf_del_index in read_index_buffer:
                    if read_index_buffer[buf_del_index] in read_buffer:
                        del read_buffer[read_index_buffer[buf_del_index]]
                    del read_index_buffer[buf_del_index]
                buf_del_index += 1

            all_flags[r.qname].add((r.flag, r.pos))

            r = scope.update(r)  # Add alignment to scope, convert to a pickle-able object

            ol_include = False
            if tree:
                if data_io.intersecter(tree, infile.get_reference_name(r.rname), r.pos - 150, r.pos + 150) and \
                   data_io.intersecter(tree, infile.get_reference_name(r.rnext), r.pnext - 150, r.pnext + 150):
                    ol_include = True

                    # Only add node if it is a split read
                    if r.has_SA:
                        G.add_node(n1, p=(r.rname, r.pos))

            # Keeps all singletons that dont overlap --include, even if no 'black' or 'grey' edges.
            if not ol_include:
                G.add_node(n1, p=(r.rname, r.pos))

            # targets = set([])
            gg = 0
            bb = 0
            yy = 0
            max_gg = 1  # 6
            max_bb = 1  # 10
            max_yy = 1  # 6
            for t, ol, sup_edge in scope.iterate():  # Other read, overlap current read, supplementary edge
                if not t:
                    break

                n2 = (t.qname, t.flag, t.pos)
                all_flags[t.qname].add((t.flag, t.pos))
                # assert n2 in name_to_node
                # Four types of edges; black means definite match with overlapping soft-clips. Grey means similar
                # rearrangement with start and end coords on reference genome.
                # Yellow means supplementary (with no good soft-clip) matches a primary
                # Finally white edge means a read1 to read2 mate-pair edge

                # inferred_edge = skip_traingles(G, n1, n2)
                # if inferred_edge:
                #     if inferred_edge == "b":
                #         bb += 1
                #     elif inferred_edge == "g":
                #         gg += 1
                #     elif inferred_edge == "y":
                #         yy += 1
                #     continue

                if bb > max_bb and gg > max_gg:
                    break

                if ol:  # and not sup_edge:
                    if bb <= max_bb:  # Stop the graph getting unnecessarily dense
                        prob_same = c_cluster_funcs.alignments_match(r.seq, t.seq,
                                    0 if r.cigartuples[0][0] != 4 else r.cigartuples[0][1],
                                    0 if t.cigartuples[0][0] != 4 else t.cigartuples[0][1],
                                    r.pos, t.pos, r.query_qualities, t.query_qualities)

                        if prob_same > 0.01:  # Add a black edge

                            G.add_node(n1, p=(r.rname, r.pos))
                            G.add_node(n2, p=(t.rname, t.pos))
                            G.add_edge(n1, n2, c="b")
                            black_edges += 1
                            bb += 1
                            continue  # Skip looking for other edges

                # Make sure both are discordant, mapped to same chromosome
                if (r.flag & 2) or (t.flag & 2) or r.rnext == -1 or r.rnext != t.rnext:
                    continue

                if abs(r.pnext - t.pnext) < max_dist:  # Other ends are within clustering distance

                    add_grey = False
                    if not ol_include:
                        # Add a grey edge if they are both discordant reads
                        if not n1[1] & 2 and not n2[1] & 2:
                            add_grey = True

                    if add_grey:
                        if gg <= max_gg:
                            G.add_node(n1, p=(r.rname, r.pos))
                            G.add_node(n2, p=(t.rname, t.pos))
                            G.add_edge(n1, n2, c="g")
                            gg += 1
                            grey_added += 1
                        continue

                elif sup_edge:
                    # Add a yellow edge
                    if yy <= max_yy:
                        G.add_node(n1, p=(r.rname, r.pos))
                        G.add_node(n2, p=(t.rname, t.pos))
                        G.add_edge(n1, n2, c="y")
                        yy += 1

    # Add read-pair information to the graph, link regions together
    new_edges = []
    for g in nx.connected_component_subgraphs(G, copy=False):
        # Add white edges between read-pairs which are NOT joined by edges in the subgraph
        # all_flags: rname: (flag, pos)
        nodes_to_check = [(n, all_flags[n]) for n, f, p in g.nodes()]

        for n, flags in nodes_to_check:  # 2 reads, or 3 if supplementary read

            for f1, f2 in itertools.combinations_with_replacement(flags, 2):  # Check all edges within template
                u, v = (n, f1[0], f1[1]), (n, f2[0], f2[1])
                has_v, has_u = g.has_node(v), g.has_node(u)

                # Only add an edge if either u or v is missing, not if both (or neither) missing
                if has_u != has_v:  # XOR

                    new_edges.append((u, v, {"c": "w"}))

    G.add_edges_from(new_edges)

    t1 = time.time() - t0
    click.echo("Built cluster graph {}s, edges={}, nodes={}".format(round(t1, 1), len(G.edges()), len(G.nodes())),
               err=True)

    return G, read_buffer


def block_model_evidence(bm, parent_graph):
    """
    Quantify the level of support over each edge in the block model
    :param bm: input block model
    :param parent_graph: the big component
    :return: annotated block model, edge attributes give evidence support
    """
    # Go through block model nodes and annotate edges
    seen = set([])

    for node_set in bm.nodes():
        if node_set in seen:

            continue
        seen.add(node_set)
        read_names_a = set([i[0] for i in node_set])

        for neighbor_set in bm[node_set]:
            read_names_b = set([i[0] for i in neighbor_set])

            pe_support = len(read_names_a.intersection(read_names_b))

            # Reads connected with black edges give soft-clip support at each side
            sub = parent_graph.subgraph(node_set.union(neighbor_set))

            black_connected = [(j[0], j[1]) for j in [i for i in sub.edges(data=True) if i[2]['c'] == 'b']]
            black_nodes = len(set(item for sublist in black_connected for item in sublist))

            supplementary = len(set([i[0] for i in sub.nodes() if i[1] & 2048]))

            res = {"pe": pe_support,
                   "sc": black_nodes,  # Todo this is overlapping soft-clips not number of soft-clips
                   "supp": supplementary}

            bm[node_set][neighbor_set]["result"] = res
            seen.add(neighbor_set)

    return bm
