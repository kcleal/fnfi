import uuid
import datetime
import itertools
import os
import multiprocessing
import time
from collections import defaultdict, Counter
from subprocess import call
import click
import networkx as nx
import numpy as np
import pysam
import sys
import difflib
import pandas as pd
import array
import caller
from sklearn.neighbors import KDTree
import data_io
import c_io_funcs
import assembler

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    xrange
except NameError:
    xrange = range


class Scoper:
    """Keeps track of which reads are in scope. Maximum distance depends on the template insert_median"""
    def __init__(self, max_dist, clip_length=20, max_depth=60):
        self.max_dist = max_dist
        self.max_depth = max_depth
        self.scope = []
        self.clip_length = clip_length

    def update(self, current_read):
        if len(self.scope) == 0:
            self.scope.append(current_read)
            return

        elif len(self.scope) > 0 and current_read.rname != self.scope[-1].rname:
            self.scope = [current_read]  # Empty scope on new chromosome
            return

        current_pos = current_read.pos

        # Go through reads and check if they are in scope
        trim_scope_idx = 0
        start = 0
        if len(self.scope) >= self.max_depth:
            start = 1  # Has the effect of dropping first entry from scope
        for scope_index in xrange(start, len(self.scope)):
            if abs(self.scope[scope_index].pos - current_pos) < self.max_dist:
                trim_scope_idx = scope_index
                break

        if trim_scope_idx > 0:
            # assert (abs(self.scope[trim_scope_idx].pos - current_pos) < self.max_dist) != \
            # (abs(self.scope[trim_scope_idx - 1].pos - current_pos) < self.max_dist)
            if trim_scope_idx < len(self.scope):
                self.scope = self.scope[slice(trim_scope_idx, len(self.scope))]

        self.scope.append(current_read)

    @staticmethod
    def overlap(start1, end1, start2, end2):
        return max(0, min(end1, end2) - max(start1, start2))

    def iterate(self):
        if len(self.scope) > 1:

            for i in xrange(len(self.scope)):
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


def align_match(r, t, max_mismatches=2):  # Todo add parameter

    def extent(a):
        start = a.pos
        end = a.pos + len(a.seq)
        if a.cigartuples[0][0] == 4:
            start -= a.cigartuples[0][1]
            end -= a.cigartuples[0][1]
        return start, end

    r_start, r_end = extent(r)
    t_start, t_end = extent(t)

    r_offset = r_start - (min(r_start, t_start))
    t_offset = t_start - (min(r_start, t_start))

    size = max(len(t.seq) + t_offset, len(r.seq) + r_offset)
    arr = np.chararray((2, size), itemsize=1)
    arr[:] = ""

    # Convert string to character array
    r_s = np.array(list(r.seq))
    t_s = np.array(list(t.seq))

    arr[0, r_offset: len(r.seq) + r_offset] = r_s
    arr[1, t_offset: len(t.seq) + t_offset] = t_s

    istart = max(r_offset, t_offset)
    iend = min(len(r.seq) + r_offset, len(t.seq) + t_offset)

    mm = arr[0, istart:iend] != arr[1, istart:iend]
    mismatches = np.sum(mm)
    if mismatches > max_mismatches:
        return 0., 0.

    matches = len(mm) - mismatches
    identity = float(matches) / (mismatches + matches)

    p = 1.0
    # Find probability of mismatches

    if mismatches > 0:
        mm_idxs = np.where(mm)[0]
        r_idxs = [i - r_offset + istart for i in mm_idxs]
        t_idxs = [i - t_offset + istart for i in mm_idxs]

        rq = [r.query_qualities[i] for i in r_idxs]
        tq = [t.query_qualities[i] for i in t_idxs]

        minqual = np.min([rq, tq], axis=0).astype(float)

        p = np.multiply.reduce(10 ** (-minqual / 10))

    return identity, p


def explore_local(starting_nodes, large_component, other_nodes, look_for, upper_bound):
    """
    Search the large component graph for additional nodes (evidence) that doesnt conflict with other nodes or the
    upper bound on copy number. If a conflict is found, no nodes of that type are returned
    :param starting_nodes: set of nodes in current nuclei
    :param large_component: the larger graph to search
    :param other_nodes: conflicting nodes already assigned to other nuclei
    :param look_for: the types of edges to look for i.e. grey or yellow
    :param upper_bound: an upper bound on the expected number of reads
    :return: a set of additional nodes to safely add to the current nuclei
    """
    additional_nodes = {k: set([]) for k in look_for}
    seen = set(starting_nodes)
    if len(starting_nodes) == 0:
        return set([])

    while True:
        nd = starting_nodes.pop()
        seen.add(nd)
        for edge in large_component.edges(nd, data=True):
            if edge[2]['c'] in additional_nodes:
                for n in (0, 1):  # Check both ends of the edge for uniqueness
                    if edge[n] not in seen:
                        if edge[n] in other_nodes:
                            del additional_nodes[edge[2]['c']]  # This evidence conflicts with another source
                            continue
                        starting_nodes.add(edge[n])
                        additional_nodes[edge[2]['c']].add(edge[n])
                        if len(additional_nodes[edge[2]['c']]) > upper_bound:  # Remove this type of evidence if over upper bound
                            del additional_nodes[edge[2]['c']]

        if len(starting_nodes) == 0:
            break
    return set().union(*additional_nodes.values())


def make_nuclei(g, _debug=False):
    """
    Nuclei are primarily clusters of overlapping soft-clipped reads. If none are found in the graph, try using
    supplementary edges, or finally grey edges only.
    Additional nuclei are then created from yellow, or grey edges, only if they do not intersect with already defined
    nuclei. The reason for doing this is to help seperate overlapping events, with some events having soft-clips but
    other events having no soft-clips. If the non-soft-clipped events are not nucleated then they will be dropped.
    :param g:
    :return:
    """
    sub_grp = nx.Graph()
    look_for = ["y", "g"]
    edge_colors = defaultdict(list)
    # t = ('simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69', 83, 46697201)
    for e in g.edges(data=True):
        # if e[0] == t or e[1] == t:
        #     echo("edge", e)
        edge_colors[e[2]['c']].append(e)

    if _debug:
        echo("sub edges", edge_colors)

    edges = []
    if "b" in edge_colors:
        edges = edge_colors["b"]
    elif "y" in edge_colors:
        edges = edge_colors["y"]
        look_for.remove("y")
    elif "g" in edge_colors:
        edges = edge_colors["g"]
    elif "w" in edge_colors:  # Sometimes single read-pair mapping event, so only white edges available
        edges = edge_colors["w"]

    # Look for other unique connected components
    sub_grp.add_edges_from(edges)

    # If no edges were added, try adding single reads, check for spanning reads
    # if len(edges) == 0:  # and len(new_edges) == 0:
    #     for n, dta in g.nodes(data=True):
    #         if dta["s"] == 1:
    #             sub_grp.add_node(n)

    return sub_grp


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


def make_block_model(g, insert_size, insert_stdev, read_length, _debug=False):
    """
    Make a block model representation of the graph. The nodes in the block model correspond to each side of a SV. The
    edge attributes give the different types of support between nodes i.e. number of split reads, soft-clips, pairs.
    The block model is then analysed for SV calling, later in the pipeline.
    The model is constructed from the overlap and paired-end input graph, using a concept of nucleation. Nuclei are
    generated from overlapping soft-clips, or if these are lacking from the 'yellow' or 'grey' edges. Additional edges
    are added to the nuclei if they can be uniquely assigned to that nuclei. For example, a nuclei of nodes connected
    with black edges, nodes with grey edges are added to the nuclei if NONE of these nuclei-specific grey edges
    are connected with a different nuclei. If a node is shared between two nuclei, this implies that kind type of
    evidence is unreliable, such as when a very dense overlapping cluster occurs and grey edges cannot be assigned
    uniquely.
    This approach proceeds by first creating an 'inner' block model consisting of black edges only. Next grey edges
    from the input graph are assigned to each node in the block-model, only if ALL of the edges can be uniquely assigned
    to that node.
    :param:
    :return: Block model, nodes correspond to break sites, edges give connectivity information with other break sites
    """
    # Make the inner block-model
    sub_grp = make_nuclei(g, _debug)

    sub_grp_cc = list(nx.connected_component_subgraphs(sub_grp))

    if _debug:
        echo("sub_grp nodes", len(sub_grp.nodes()))
        echo("sub_grp components", len(sub_grp_cc))

    # Drop any reads that are'nt in a connected component
    intersection = g.subgraph(sub_grp.nodes())

    inner_model = blockmodel(intersection, partitions=[i.nodes() for i in sub_grp_cc])

    # Each node in the inner_model is actually a set of nodes.
    # For each node in the model, explore the input graph for new edges
    # return inner_model

    all_node_sets = set(inner_model.nodes())
    updated = True
    updated_partitons = []
    for node_set in inner_model.nodes():

        other_nodes = set([item for sublist in list(all_node_sets.difference(set([node_set]))) for item in sublist])

        # Expected local nodes
        # Read clouds uncover variation in complex regions of the human genome. Bishara et al 2016.
        # max template size * estimated coverage / read_length; 2 just to increase the upper bound a bit
        local_upper_bound = ((insert_size + (2 * insert_stdev)) * float(2 + len(node_set))) / float(read_length)
        additional_nodes = explore_local(set(node_set), g, other_nodes, ["y", "g"], local_upper_bound)
        if len(additional_nodes) > 0:  # If any new nodes are found
            updated = True
        updated_partitons.append(set(node_set).union(additional_nodes))

    if updated:  # Re-partition the graph
        intersection = g.subgraph([item for sublist in updated_partitons for item in sublist])
        inner_model = blockmodel(intersection, partitions=updated_partitons)

    return inner_model


def block_model_evidence(bm, parent_graph):
    """
    Quantify the level of support over each edge in the block model
    :param bm: input block model
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

            # total_reads = len(node_set) + len(neighbor_set)
            pe_support = len(read_names_a.intersection(read_names_b))

            # Reads connected with black edges give soft-clip support at each side
            sub = parent_graph.subgraph(node_set.union(neighbor_set))
            black_connected = [(j[0], j[1]) for j in [i for i in sub.edges(data=True) if i[2]['c'] == 'b']]
            black_nodes = len(set(item for sublist in black_connected for item in sublist))

            supplementary = len([i[1] for i in sub.nodes() if i[1] & 2048])

            res = {# "total_reads": total_reads,
                   "pe": pe_support,
                   "sc": black_nodes,  # Todo this is overlapping soft-clips not number of soft-clips
                   "supp": supplementary}

            # if ('simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69', 83, 46697201) in node_set:
            #     echo("nodes set", node_set)
            #     echo("neigh", neighbor_set)
            #     echo("union", node_set.union(neighbor_set))
            #     echo(sub.nodes())
            #     ec = []
            #     for i in sub.edges(data=True):
            #         ec.append(i[2]['c'])
            #     echo(Counter(ec))
            #     echo(res)
            #     quit()

            bm[node_set][neighbor_set]["result"] = res
            seen.add(neighbor_set)
    return bm


def get_reads(infile, sub_graph, max_dist, rl, read_buffer):
    """Using infile, collect reads belonging to sub_graph.
    :param infile: the input file
    :param sub_graph: the subgraph of interest, nodes are reads to collect
    :param max_dist: used to define an interval around reads of interest; this interval is then searched using pysam
    :param rl: padding to use, approximately read length
    :returns read dictionary, u nodes, v nodes, edge data"""

    rd = defaultdict(lambda: defaultdict(dict))  # rname: (flag, pos): alignment
    c = 0
    coords = []
    for node, dta in sub_graph.nodes(data=True):
        if node in read_buffer:
            rd[node[0]][(node[1], node[2])] = read_buffer[node]
            c += 1
        else:
            coords.append(dta['p'])

    if len(coords) > 0:  # fetch any reads that were'nt in the buffer
        coords = itertools.groupby(sorted(set(coords), key=lambda x: (x[0], x[1])), key=lambda x: x[0])  # Groupby chrom

        # Make intervals of reads that are near one another
        # https://stackoverflow.com/questions/10016802/python-group-a-list-of-integer-with-nearest-values
        for chrom, crds in coords:
            d = [i[1] for i in crds]
            chrom = infile.get_reference_name(chrom)
            m = [[chrom, d[0] - rl, d[0] + 2 * rl]]

            for x in d[1:]:
                if x - m[-1][-1] < max_dist:
                    m[-1][-1] = x + 2 * rl
                else:
                    m.append([chrom, x - rl, x + 2*rl])

            # Collect reads
            for chrom, start, end in m:
                for aln in infile.fetch(chrom, 1 if start <= 0 else start, end):
                    if (aln.qname, aln.flag, aln.pos) in sub_graph:
                        if aln.flag is not None:
                            rd[aln.qname][(aln.flag, aln.pos)] = aln
                            c += 1
    # assert(c == len(nodes))  #!
    if c != len(sub_graph.nodes()):
        click.echo("WARNING: reads missing from collection - {} vs {}".format(c, len(sub_graph.nodes())), err=True)

    return rd


def construct_graph(infile, max_dist, tree, buf_size=100000, input_windows=(), window_pad=1000):
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

    approx_read_length = []

    G = nx.Graph()
    winds = [(c, s - window_pad if s - window_pad > 0 else 0, e + window_pad) for c, s, e in input_windows]
    for widx, window in enumerate(winds):

        for r in infile.fetch(*window):

            if len(r.cigar) == 0 or r.flag & 4:
                continue

            n1 = (r.qname, r.flag, r.pos)

            if len(approx_read_length) < 100:
                if not r.flag & 2048:
                    if r.seq:
                        approx_read_length.append(len(r.seq))

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
            scope.update(r)  # Add alignment to scope

            ol_include = False
            if tree:
                if data_io.intersecter(tree, infile.get_reference_name(r.rname), r.pos - 150, r.pos + 150) and \
                   data_io.intersecter(tree, infile.get_reference_name(r.rnext), r.pnext - 150, r.pnext + 150):
                    ol_include = True

            # Keeps all singletons that dont overlap --include, even if no 'black' or 'grey' edges.
            if not ol_include:
                G.add_node(n1, {"p": (r.rname, r.pos)})

            gg = 0
            bb = 0
            yy = 0

            # Todo parallelize this loop:
            for t, ol, sup_edge in scope.iterate():  # Other read, overlap current read, supplementary edge
                if not t:
                    break

                n2 = (t.qname, t.flag, t.pos)
                all_flags[t.qname].add((t.flag, t.pos))

                # Four types of edges; black means definite match with overlapping soft-clips. Grey means similar
                # rearrangement with start and end co-ords on reference genome.
                # Yellow means supplementary matches a primary
                # read; these edges need to be checked when both primary reads are available
                # Finally white edge means a read1 to read2 edge
                if ol:  # and not sup_edge:
                    identity, prob_same = align_match(r, t)
                    # if r.qname == "simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69":
                    #     echo("black edge", identity, prob_same, "r", n1, n2)
                    # if t.qname == "simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69":
                    #     echo("black edge", identity, prob_same, "t", n1, n2)
                    if prob_same > 0.01:  # Add a black edge
                        G.add_edge(n1, n2, {"c": "b"})
                        black_edges += 1
                        bb += 1
                        if bb > 6:
                            break
                    continue

                # Make sure both are discordant, mapped to same chromosome
                if (r.flag & 2) or (t.flag & 2) or r.rnext == -1 or r.rnext != t.rnext:
                    continue

                if abs(r.pnext - t.pnext) < max_dist:  # Other ends are within clustering distance
                    # Skip if read is not soley within --include regions
                    if ol_include and tree is not None:
                        continue

                    # Add a grey edge if they are both discordant reads
                    if not n1[1] & 2 and not n2[1] & 2:
                        G.add_edge(n1, n2, {"c": "g"})
                        gg += 1
                        grey_added += 1
                        if gg > 6:
                            break  # Stop the graph getting unnecessarily dense
                        continue

                elif sup_edge:
                    # Add a yellow edge
                    G.add_edge(n1, n2, {"c": "y"})
                    yy += 1
                    if yy > 6:
                        break
                    continue

    # Add read-pair information to the graph, link regions together
    new_edges = []
    for g in nx.connected_component_subgraphs(G, copy=False):

        # Add white edges between read-pairs which are NOT in the subgraph
        # all_flags: rname: (flag, pos)
        nodes_to_check = [(n, all_flags[n]) for n, f, p in g.nodes()]
        for n, flags in nodes_to_check:  # 2 reads, or 3 if supplementary read

            for f1, f2 in itertools.combinations_with_replacement(flags, 2):  # Check all edges within template
                u, v = (n, f1[0], f1[1]), (n, f2[0], f2[1])
                has_v, has_u = g.has_node(v), g.has_node(u)
                # if g.has_node(('simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69', 2145, 46696345)):
                #     echo(has_u, has_v, has_u != has_v, u, v)
                # Only add an edge if either u or v is missing, not if both (or neither) missing
                if has_u != has_v:  # XOR
                    new_edges.append((u, v, {"c": "w"}))

        # if g.has_node(('simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69', 2145, 46696345)):  # todo find out which no white edges added
        #     ec = []
        #     echo(nodes_to_check)
        #     for item in g.edges(data=True):
        #         ec.append(item[2]['c'])
        #     echo("!subg size", len(g.nodes()), Counter(ec))
        #     quit()

    # Regroup based on white edges (link together read-pairs)
    G.add_edges_from(new_edges)
    if len(approx_read_length) > 0:
        approx_read_length = int(np.mean(approx_read_length))
    else:
        approx_read_length = 150
    click.echo("Built cluster graph", err=True)
    return G, read_buffer, approx_read_length


def merge_intervals(intervals):
    """
    >>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [['chr1', 1, 5], ['chr2', 3, 5]]
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))  # by chrom, start, end
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
            continue
        elif higher[0] != merged[-1][0]:  # Dont merge intervals on different chroms
            merged.append(higher)
        else:
            lower = merged[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= lower[2]:
                merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]))
            else:
                merged.append(higher)
    return merged


def knn_cluster(include, infile, max_dist):
    """
    Utilize k nearest neighbour to cluster discordant read-pairs; uses a 2d genome representation
    :param include: regions to limit search to. Only regions with an overlap with 'include' are returned.
    :param infile: input alignment file
    :param max_dist: maximum clustering distance (manhattan distance)
    :return: region U, region V, sug_graph_id, sub_graph_nodes, nodes in U, nodes in V
    """
    # regions = data_io.overlap_regions(args["include"])  # Todo add search and limit options
    # list(u), list(v), sub_graph_id, sug_nodes, links_d[u], links_d[v]
    cluster_map = defaultdict(lambda: array.array('L'))

    # Cluster reads by genome pos, using a 2d plane, pos1 --> pos2
    # read_name --> Add all supplementary, primary only once (use mate pair info), no secondary
    primary_seen = set([])
    count = 0

    # Only need one read of the pair to make the cluster map, so sufficient to just search the 'include' regions
    infile_itr = data_io.get_include_reads(include, infile)

    graph_key_map = defaultdict(list)
    for r in infile_itr:
        f = r.flag
        # fails quality checks, PCR duplicate, unmapped, mate unmapped,  not primary, proper pair
        if f & 1806:
            continue

        # Skip if both positions are non-canonical chromosomes
        if len(infile.get_reference_name(r.rname)) > 5 and len(infile.get_reference_name(r.rnext)) > 5:
            continue
        # Only primary and supplementary left
        name = r.qname
        if name in primary_seen:
            continue

        # Make sure chrom name is in sort order
        cmap_key = tuple(sorted([r.rname, r.rnext]))
        if r.rname == cmap_key[0]:
            pos1 = r.pos
            pos2 = r.pnext
        else:
            pos2 = r.pos
            pos1 = r.pnext

        cluster_map[cmap_key].append(pos1)
        cluster_map[cmap_key].append(pos2)

        # Only save graph node info for one read per pair
        # This value will be identical to the node-name in the main graph later on:
        graph_key_map[cmap_key].append((name, r.flag, r.pos))

        primary_seen.add(name)
        count += 1

    if len(primary_seen) == 0:
        raise ValueError("No reads in input file")

    del primary_seen  # Garbage collect

    sub_graph_id = 0

    # For each chromosome pair. Cluster k-NN on a graph. Nodes are INDEXES into the positions array

    clusters = {}  # (cmap_key, n1, n2): {subg_id, n_nodes, out_edges_n1, out_edges_n2}
    pos_d = defaultdict(tuple)  # The parray positions  (cmap_key, n1, n2): (chr1_pos, chr2_pos)

    for cmap_key, positions in cluster_map.iteritems():

        parray = np.array(positions, dtype=np.uint32).reshape((-1, 2))
        chr1, chr2 = infile.get_reference_name(cmap_key[0]), infile.get_reference_name(cmap_key[1])

        if chr1 == chr2 and len(chr1) < 6:
            click.echo("Calculating NN; links={} {}<-->{}".format(len(parray), chr1, chr2), err=True)

        neigh = KDTree(parray, leaf_size=2, metric="manhattan")

        G = nx.Graph()
        max_nn = 10 if 10 < len(parray) else len(parray)
        min_nn = 3 if 3 < max_nn else max_nn
        numb_neigh = max_nn  # Add edges between n number of neighbours to a point
        for index, (x, y) in enumerate(parray):
            G.add_node(index)

            dists, rng = neigh.query([[x, y]], k=numb_neigh)
            dists = dists[0]
            rng = rng[0]

            if dists[-1] < max_dist:
                numb_neigh -= 1
            else:
                numb_neigh += 1
            if numb_neigh < min_nn:
                numb_neigh = min_nn
            elif numb_neigh > max_nn:
                numb_neigh = max_nn

            for index_2, d2 in zip(rng, dists):
                if index == index_2:  # Skip self
                    continue
                if d2 > max_dist:  # Ensure within right distance
                    continue

                G.add_edge(index, index_2)

        # Need to merge any overlapping intervals, otherwise mutliple calls can be made
        # Only merge intervals locally (within the chromosome-pair, not between all chromosomes)
        cc = 0

        chr1_itv = []
        chr2_itv = []
        links_d = defaultdict(int)  # The number of links (out edges) from each node

        gnodes = graph_key_map[cmap_key]
        for sub in nx.connected_component_subgraphs(G):

            # Get the corresponding positions pos1, pos2
            chr1_pos, chr2_pos = zip(*[parray[i] for i in sub.nodes()])

            # Seed nodes used for indexing the main graph, quick way to find components regions of interest;
            # i.e. use a BFS using these nodes to find components (black-edge subgraphs)
            seed_nodes = [gnodes[i] for i in sub.nodes()]
            # if
            links = len(sub.nodes())

            # Bounding interval for pos1 and pos2
            chr1_start, chr1_end = min(chr1_pos), max(chr1_pos)
            chr2_start, chr2_end = min(chr2_pos), max(chr2_pos)

            # Add interval node
            n1 = (chr1, chr1_start, chr1_end)
            n2 = (chr2, chr2_start, chr2_end)

            # Save the parray information, useful for indexing graph later on
            n1, n2 = sorted([n1, n2], key=lambda x: (infile.get_tid(x[0]), int(x[1])))
            pos_d[(cmap_key, n1, n2)] = seed_nodes

            chr1_itv.append(n1)
            chr2_itv.append(n2)
            links_d[n1] += links
            links_d[n2] += links
            cc += 1

        # Get connected components for this chromosome pair
        ol_G = nx.Graph()
        ol_G.add_edges_from(zip(chr1_itv, chr2_itv))

        for subg in nx.connected_component_subgraphs(ol_G, copy=False):
            sub_graph_id += 1
            sug_nodes = len(subg.nodes())
            for u, v in subg.edges():

                if len(u[0]) <= 5 or len(v[0]) <= 5:  # Note no calls between non-canonical chroms
                    # Might need to re-order (u, v) as (v, u) for cluster map indexing
                    u, v = sorted([u, v], key=lambda x: (infile.get_tid(x[0]), int(x[1])))
                    key = (cmap_key, u, v)
                    clusters[key] = {"gid": sub_graph_id, "gnodes": sug_nodes, "ue": links_d[u], "ve": links_d[v]}

    for k in pos_d.keys():  # Drop unneeded pos_d values
        if k not in clusters:
            del pos_d[k]

    return pos_d, clusters


def analyse_connectivity(f):

    node_counts = []
    for cmapk, n1, n2 in f.keys():
        node_counts.append(n1)
        node_counts.append(n2)

    # Number of times node appears is equal to how many other regions it is connected to
    neighbour_counts = Counter(node_counts)

    for cmapk, n1, n2 in f.keys():
        n = (cmapk, n1, n2)

        u_out = neighbour_counts[n1] - 1
        v_out = neighbour_counts[n2] - 1

        f[n]["n1_connectivity"] = u_out
        f[n]["n2_connectivity"] = v_out

        unodes = f[n]["ue"]
        vnodes = f[n]["ve"]

        # A measure of how likely this edge is the 'real' edge, smaller is better
        f[n]["connectivity"] = (((float(u_out) / unodes) + 1) * (float(v_out) / vnodes) + 1)

    return f


def get_region_depth(cm, inputbam, dest, unique_file_id, pad=10000, delete=True):
    intervals = []
    for (cmk, n1, n2) in cm.keys():
        c1, s1, e1 = n1
        c2, s2, e2 = n2

        s1 = int(s1) - pad
        s2 = int(s2) - pad
        e1 = int(e1) + pad
        e2 = int(e2) + pad

        s1 = (int(s1) / 100) * 100
        if s1 < 0:
            s1 = 0
        s2 = (int(s2) / 100) * 100
        if s2 < 0:
            s2 = 0
        e1 = ((int(e1) / 100) * 100) + 100
        e2 = ((int(e2) / 100) * 100) + 100
        intervals += [[c1, s1, e1], [c2, s2, e2]]

    temp_1 = "{}/regionstempdin.{}.csv".format(dest, unique_file_id)
    temp_2 = "{}/regionstempdout.{}.csv".format(dest, unique_file_id)

    intervals = merge_intervals(intervals)
    with open(temp_1, "w") as depthf:
        for chrom, start, end in intervals:
            starts = np.arange(start, end, 100).astype(str)
            ends = np.arange(start + 100, end + 100, 100).astype(str)
            k = [chrom] * len(starts)
            lines = ["\t".join(item) + "\n" for item in zip(k, starts, ends)]
            depthf.writelines(lines)
    call("samtools bedcov {t1} {bam} > {t2}".format(bam=inputbam, t1=temp_1, t2=temp_2), shell=True)

    depth_d = {}
    with open(temp_2, "r") as depth_i:
        for line in depth_i:
            chrom, s, e, depth = line.strip().split("\t")
            depth_d[(chrom, int(s))] = float(depth) / 100.

    # Remove temp files
    if delete:
        os.remove(temp_1)
        os.remove(temp_2)

    return depth_d


def filter_high_cov(cm, pm, regions_depth, tree, max_cov):

    removed = 0
    bad_keys = set([])
    for kk in cm.keys():
        region_pair = (kk[1], kk[2])
        max_covs = []
        for chrom, s, e in region_pair:
            # Check for intersection with --include bed
            ol = data_io.intersecter(tree, chrom, s, e)
            if ol:
                continue
            max_covs.append(max(caller.calculate_coverage(chrom, s, e, regions_depth)))

        if any([i > max_cov for i in max_covs]):
            removed += 1
            bad_keys.add(kk)
            continue

    for k in bad_keys:
        del cm[k]
        del pm[k]

    click.echo("Dropped {} high coverage regions. Kept {}".format(removed, len(cm)), err=True)
    return cm, pm


def get_regions_windows(args, unique_file_id, infile, max_dist, tree, delete=True):
    """
    :param args:
    :param unique_file_id:
    :param infile:
    :param max_dist:
    :param tree:
    :return:
    """
    position_map, cluster_map = knn_cluster(args["include"], infile, max_dist)

    with open("rmap.txt", "w")as rm:
        for ck, u, v in cluster_map.keys():
            rm.write("{}\t{}\t{}\tu\n".format(*u))
            rm.write("{}\t{}\t{}\tv\n".format(*v))

    target = None  # "simulated_reads.intra.0.2-id7_A_chr21:46699155_B_chr21:46696339-215"
    tk = None
    for k, v in position_map.items():
        for n, f, p in v:
            if n == target:
                echo(k, n)
                tk = k

    cluster_map = analyse_connectivity(cluster_map)
    regions_depth = get_region_depth(cluster_map, args["sv_aligns"], args["dest"], unique_file_id, pad=10000, delete=delete)
    cluster_map, position_map = filter_high_cov(cluster_map, position_map, regions_depth, tree, args["max_cov"])

    for k, v in position_map.items():
        for n, f, p in v:
            if n == target:
                echo("get_region_windows still here", k, (n, f, p))
    return cluster_map, position_map, regions_depth, tk


def get_component_from_seed(G, seed_nodes):

    # Explore all out edges from original seed nodes
    nodes_to_visit = set(seed_nodes)
    found_nodes = set(seed_nodes)
    seen = set([])
    while len(nodes_to_visit) > 0:
        nd = nodes_to_visit.pop()
        if nd in seen:
            continue
        seen.add(nd)
        for edge in G.edges(nd, data=True):
            v = edge[1]
            if v not in seen:
                found_nodes.add(v)
                #if edge[2]['c'] == "white":
                nodes_to_visit.add(v)
    return G.subgraph(found_nodes)


def filter_potential(input_events, tree):
    potential = []

    for i in input_events:
        if "posB" not in i:  # Skip events for which no posB was identified
            continue
        if "contig" not in i or i["contig"] == "":
            i["contig"] = None

        # Remove events for which both ends are in --include but not contig was found
        if (data_io.intersecter(tree, i["chrA"], i["posA"], i["posA"] + 1) and
            data_io.intersecter(tree, i["chrB"], i["posB"], i["posB"] + 1)) and i["contig"] is None:
            continue

        potential.append(i)
    return potential


def merge_events(potential, max_dist, tree, try_rev=False):

    if len(potential) <= 1:
        return potential

    G = nx.Graph()
    seen = set([])  # Tested edges
    printy = False

        #
        #
        # total_matched = shortest - result["editDistance"]
        # # total_matched = sum([blk.size for blk in matching_blocks])
        # identity = total_matched / float(shortest)
        # return identity

    for i in range(len(potential)):
        for j in range(len(potential)):
            ei = potential[i]
            ej = potential[j]

            # if ei["posA"] == 43782665 or ej["posA"] == 43782665 or ei["posB"] == 43782665 or ej["posB"] == 43782665:
            #     printy = True

            if i == j or (i, j) in seen or (j, i) in seen:
                continue

            seen.add((i, j))

            ei = potential[i]
            ej = potential[j]

            # Check if events point to the same loci
            loci_same = False
            if ei["chrA"] == ej["chrA"]:  # Try chrA matches chrA
                if abs(ei["posA"] - ej["posA"]) < max_dist:
                    if ei["chrB"] == ej["chrB"]:
                        if abs(ei["posB"] - ej["posB"]) < max_dist:
                            loci_same = True
            if not loci_same:  # Try chrA matches chrB
                if ei["chrA"] == ej["chrB"]:
                    if abs(ei["posA"] - ej["posB"]) < max_dist:
                        if ei["chrB"] == ej["chrA"]:
                            if abs(ei["posB"] - ej["posA"]) < max_dist:
                                loci_same = True
            if "contig" in ei and "contig" in ej:
                ci, cj = str(ei["contig"]).upper(), str(ej["contig"]).upper()
            else:
                continue

            # if "AATGGTGATTTCCTTTATAGATGTAAATTTTC" in ci and "AATGGTGATTTCCTTTATAGATGTAAATTTTC" in cj:
            #     echo("hi", loci_same)
            #     echo(ci)
            #     echo(cj)
            #     echo(assembler.check_contig_match(ci, cj))
            #
            #     matcher = difflib.SequenceMatcher(None, ci, cj)
            #     matching_blocks = matcher.get_matching_blocks()
            #     echo(matching_blocks)
            #     echo(len(ci), len(cj))
            #     printy = True
            # else:
            #     printy = False

            if loci_same:
                if ci != "None" and len(ci) > 0 and cj != "None" and len(cj) > 0:  # Both have contigs
                    idt = assembler.check_contig_match(ci, cj)
                    if idt == 1:
                        G.add_edge(i, j)
                    elif try_rev:
                        if assembler.check_contig_match(ci, c_io_funcs.reverse_complement(cj, len(cj))) == 1:
                            G.add_edge(i, j)
                if printy:
                    echo("idt", idt)
                # Only merge loci if they are not within --include regions
                elif not (data_io.intersecter(tree, ei["chrA"], ei["posA"], ei["posA"] + 1) and
                          data_io.intersecter(tree, ei["chrB"], ei["posB"], ei["posB"] + 1)):
                    G.add_edge(i, j)

    found = []
    for i in range(len(potential)):  # Add singletons, non-merged
        if not G.has_node(i):
            found.append(potential[i])
            if printy:
                echo("1", potential[i])

    for grp in nx.connected_component_subgraphs(G):
        c = [potential[n] for n in grp.nodes()]
        if printy:
            echo(grp.nodes())
        best = sorted(c, key=lambda x: sum([x["pe"], x["supp"]]), reverse=True)
        w0 = best[0]["pe"] + best[0]["supp"]  # Weighting for base result
        for k in range(1, len(best)):
            for t in ["pe", "supp", "sc", "SP", "block_edge", "joinA", "joinB"]:  # Sum these
                best[0][t] += best[k][t]

            if best[k]["maxASsupp"] > best[0]["maxASsupp"]:
                best[0]["maxASsupp"] = best[k]["maxASsupp"]

            for t in ["DN", "MAPQsupp", "MAPQpri", "DApri", "DAsupp", "DP", "NMpri", "NMsupp"]:  # Average these
                w = best[k]["pe"] + best[k]["supp"]
                denom = w0 + w
                if denom == 0:
                    weighted_av = 0
                else:
                    weighted_av = ((best[0][t] * w0) + (best[k][t] * w)) / denom
                best[0][t] = weighted_av
            w0 = best[0]["pe"] + best[0]["supp"]
        found.append(best[0])
        if printy:
            echo(best[0])
    return found


def get_prelim_events(G, read_buffer, cluster_map, positions_map, infile, args, max_dist, regions, regions_depth, approx_rl, _debug_k):

    # Go through edges in cluster map and make calls on block-model edges/nodes
    # Then merge all events to prevent duplication

    block_edge_events = []
    for c_edge in cluster_map.keys():
        printy = False
        if c_edge == _debug_k:
            echo("here")
            printy = True
        # if c_edge == ((19, 19), ('chr17', 113572, 114777), ('chr17', 115171, 116903)):
        #     printy = True
        #     echo("here2")
        seed_reads = positions_map[c_edge]
        grp = get_component_from_seed(G, seed_reads)

        # Could add a step to check for region-spanning reads

        reads = get_reads(infile, grp, max_dist, approx_rl, read_buffer)
        bm = make_block_model(grp, args["insert_median"], args["insert_stdev"], approx_rl, _debug=printy)

        if printy:
            echo("get_prelim_events")
            echo("seeds", seed_reads)
            echo("grp len", len(grp.nodes()))
            echo("grp nodes:", grp.nodes())
            # quit()
            echo("reads collected", len(reads))
            echo("bm nodes", len(bm.nodes()))
        if len(bm.nodes()) == 0:
            continue

        bm = block_model_evidence(bm, grp)  # Annotate block model with evidence
        potential_events = []
        for event in caller.call_from_block_model(bm, grp, reads, infile, args["clip_length"],
                                                  args["insert_median"], args["insert_stdev"], printy):

            if event:
                if printy:
                    item = event
                    echo("-->", "{}:{}".format(item["chrA"], item["posA"]), "{}:{}".format(item["chrB"], item["posB"]))

                potential_events.append(event)

                if c_edge == _debug_k or printy:
                    item = event
                    if "contig" in item and "posB" in item:
                        echo(event)
                        echo("output event!", "{}:{}".format(item["chrA"], item["posA"]), "{}:{}".format(item["chrB"], item["posB"]), item["contig"])

        potential_events = filter_potential(potential_events, regions)
        potential_events = merge_events(potential_events, max_dist, regions)

        if printy:
            echo("merged block-nodes")
            for item in potential_events:
                echo("--2>", "{}:{}".format(item["chrA"], item["posA"]), "{}:{}".format(item["chrB"], item["posB"]), item["contig"])
        block_edge_events += potential_events

    # Perform a final merge across block nodes

    merged = merge_events(block_edge_events, max_dist, regions, try_rev=True)
    if merged:
        for event in merged:
            if event["posA"] == 46699837 or event["posB"] == 46699837:
                echo("pe", event)
            # Collect coverage information
            event_string = caller.get_raw_cov_information(event, regions, cluster_map[c_edge], regions_depth)
            if event_string:

                yield event_string


def cluster_reads(args):
    t0 = time.time()

    data_io.mk_dest(args["dest"])
    if args["dest"] is None:
        args["dest"] = "."

    kind = args["sv_aligns"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "rs"}

    click.echo("Input file is {}, (.{} format). {} processes".format(args["sv_aligns"], kind, args["procs"]), err=True)
    infile = pysam.AlignmentFile(args["sv_aligns"], opts[kind])

    max_dist = int(args["insert_median"] + (5 * args["insert_stdev"]))  # > distance reads drop out of clustering scope
    click.echo("Maximum clustering distance is {}".format(max_dist), err=True)

    if args["svs_out"] == "-" or args["svs_out"] is None:
        click.echo("SVs output to stdout", err=True)
        outfile = sys.stdout
    else:
        click.echo("SVs output to {}".format(args["svs_out"]), err=True)
        outfile = open(args["svs_out"], "w")

    head = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "cipos95A", "cipos95B",
            "DP", "DApri", "DN", "NMpri", "SP", "DAsup",
            "NMsup", "maxASsup", "contig", "pe", "supp", "sc", "block_edge", "MAPQpri", "MAPQsup", "raw_reads_10kb",
            "kind", "connectivity", "linked"]  # Todo add strandadness

    unique_file_id = str(uuid.uuid4())

    regions = data_io.overlap_regions(args["include"])

    #unique_file_id = "fnfi_tmp"
    cluster_map, positions_map, regions_depth, _debug_k = get_regions_windows(args, unique_file_id, infile, max_dist, regions)

    # Collect all intervals to search
    all_intervals = []
    for cmk, n1, n2 in cluster_map.keys():
        all_intervals.append(n1)
        all_intervals.append(n2)

    # Add all 'include' intervals
    all_intervals += data_io.get_bed_regions(args["include"])

    all_m = merge_intervals(all_intervals)

    G, read_buffer, approx_rl = construct_graph(infile, max_dist, regions, input_windows=all_m, buf_size=args["buffer_size"])

    prelim_ev = get_prelim_events(G, read_buffer, cluster_map, positions_map, infile, args, max_dist, regions,
                                  regions_depth, approx_rl, _debug_k)

    # temp = "{}/callstemp.{}.csv".format(args["dest"], unique_file_id)
    # with open(temp, "r") as call_t:
    #     for event in prelim_ev:
    #         call_t.write(event)
    #
    # os.remove(temp)
    c = 0
    with outfile:
        outfile.write("\t".join(head) + "\n")
        for event in prelim_ev:
            outfile.write(event)
            c += 1

    click.echo("call-events complete, n={}, {} h:m:s".format(c, str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)
