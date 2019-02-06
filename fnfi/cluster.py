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
from annoy import AnnoyIndex
import pandas as pd
import array
import caller
import data_io
import pickle

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


def make_nuclei(g):
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
    for e in g.edges(data=True):
        edge_colors[e[2]['c']].append(e)

    edges = []
    if "b" in edge_colors:
        edges = edge_colors["b"]
    elif "y" in edge_colors:
        edges = edge_colors["y"]
        look_for.remove("y")
    elif "g" in edge_colors:
        edges = edge_colors["g"]
        look_for = []
    else:
        pass  # Single read, no connections

    # Look for other unique connected components
    sub_grp.add_edges_from(edges)

    new_edges = []
    for color in look_for:
        if color in edge_colors:
            query_clusters = nx.Graph()
            query_clusters.add_edges_from(edge_colors[color])
            for cc in nx.connected_component_subgraphs(query_clusters):
                if all(i not in sub_grp for i in cc.nodes()):  # Check for overlap with a different nuclei
                    new_edges += list(cc.edges(data=True))

    sub_grp.add_edges_from(new_edges)

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


def make_block_model(g, insert_size, insert_stdev, read_length):
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
    :param subg:
    :return: Block model, nodes correspond to break sites, edges give connectivity information with other break sites
    """
    # Make the inner block-model
    sub_grp = make_nuclei(g)

    sub_grp_cc = list(nx.connected_component_subgraphs(sub_grp))

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

            total_reads = len(node_set) + len(neighbor_set)
            pe_support = len(read_names_a.intersection(read_names_b))

            # Reads connected with black edges give soft-clip support at each side
            sub = parent_graph.subgraph(node_set.union(neighbor_set))
            black_connected = [(j[0], j[1]) for j in [i for i in sub.edges(data=True) if i[2]['c'] == 'b']]
            black_nodes = len(set(item for sublist in black_connected for item in sublist))
            supplementary = len([i[1] for i in sub.nodes() if i[1] & 2048])

            res = {"total_reads": total_reads,
                   "pe": pe_support,
                   "sc": black_nodes,
                   "sup": supplementary}

            bm[node_set][neighbor_set]["result"] = res
            seen.add(neighbor_set)
    return bm


def get_reads(infile, sub_graph, max_dist, rl, read_buffer):
    """Using infile, collect reads belonging to sub_graph.
    :param infile: the input file
    :param sub_graph: the subgraph of interest, nodes are reads to collect
    :param max_dist: used to define an interval around reads of interest; this interval is then searched using pysam
    :param rl: the read length
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


def construct_graph(args, infile, max_dist, buf_size=10000, input_windows=()):  # todo add option for read buffer length

    regions = data_io.overlap_regions(args["include"])

    # nodes = []
    # edges = []

    all_flags = defaultdict(set)  # Linking clusters together rname: set([(flag, pos)..])

    # Make a buffer of reads to help prevent file lookups later
    read_buffer = dict()  # keys are (rname, flag, pos)
    read_index_buffer = dict()  # keys are int, values are (rname, flag, pos)

    scope = Scoper(max_dist=max_dist)

    count = 0
    buf_del_index = 0

    grey_added = 0
    black_edges = 0
    chroms_seen = set([])

    lastc_time = time.time()
    lastc = None

    G = nx.Graph()
    for window in input_windows:

        for r in infile.fetch(*window):

            # if count % 1000 == 0:
            #     cc = infile.get_reference_name(r.rname)
            #     if cc not in chroms_seen:
            #         chroms_seen.add(cc)
            #         click.echo("Working on {}".format(cc), err=True)
            #         t = time.time()
            #         if lastc:
            #             click.echo("{} took {}s".format(lastc, t - lastc_time), err=True)
            #             lastc_time = t
            #         lastc = cc
            #
            #     if count % 50000 == 0:
            #         if count != 0:
            #             click.echo("{}:{}, alignnmets {}, b-edges {}, g-edges {}".format(cc, r.pos, count, black_edges, grey_added), err=True)

            # Limit to regions
            if regions:
                rname = infile.get_reference_name(r.rname)  # Todo do a fetch, rather than scan the whole bam file
                if not data_io.intersecter(regions, rname, r.pos, r.pos + 1):
                    rnext = infile.get_reference_name(r.rnext)
                    if not data_io.intersecter(regions, rnext, r.pnext, r.pnext + 1):
                        continue

            if len(r.cigar) == 0:
                continue

            if r.flag & 4:  # Unmapped read
                continue

            n1 = (r.qname, r.flag, r.pos)

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

            G.add_node(n1, {"p": (r.rname, r.pos)})
            # nodes.append((n1, {"p": (r.rname, r.pos)}))
            all_flags[r.qname].add((r.flag, r.pos))

            scope.update(r)  # Add alignment to scope

            gg = 0
            bb = 0
            yy = 0
            for t, ol, sup_edge in scope.iterate():  # Other read, overlap current read, supplementary edge

                if not t:
                    break

                n2 = (t.qname, t.flag, t.pos)
                all_flags[t.qname].add((t.flag, t.pos))

                # Four types of edges; black means definite match with overlapping soft-clips. Grey means similar
                # rearrangement with start and end co-ords on reference genome. Yellow means supplementary matches a primary
                # read; these edges need to be checked when both primary reads are available
                # Finally white edge means a read1 to read2 edge
                if ol and not sup_edge:
                    identity, prob_same = align_match(r, t)

                    if prob_same > 0.01:  # Add a black edge
                        # e = (n1, n2, {"c": "b"})
                        G.add_edge(n1, n2, {"c": "b"})
                        # edges.append(e)
                        black_edges += 1
                        bb += 1
                        if bb > 3:
                            break
                    continue

                # elif not ol and not sup_edge:
                # Make sure both are discordant, mapped to same chromosome
                if (r.flag & 2) or (t.flag & 2) or r.rnext == -1 or r.rnext != t.rnext:
                    continue

                if abs(r.pnext - t.pnext) < max_dist:
                    # Add a grey edge
                    # e = (n1, n2, {"c": "g"})
                    G.add_edge(n1, n2, {"c": "g"})
                    gg += 1
                    grey_added += 1
                    if gg > 3:
                        break  # Stop the graph getting unnecessarily dense
                    # edges.append(e)
                    continue

                elif sup_edge:
                    # Add a yellow edge
                    #e = (n1, n2, {"c": "y"})
                    #edges.append(e)
                    G.add_edge(n1, n2, {"c": "y"})
                    yy += 1
                    if yy > 3:
                        break
                    continue

    # Add read-pair information to the graph, link regions together
    #echo("Constructing graph")

    # G.add_nodes_from(nodes)
    # G.add_edges_from(edges)
    new_edges = []
    #echo("Adding read-pair information")
    for g in nx.connected_component_subgraphs(G, copy=True):

        # Add white edges between read-pairs which are NOT in the subgraph
        # all_flags: rname: (flag, pos)
        nodes_to_check = [(n, all_flags[n]) for n, f, p in g.nodes()]
        for n, flags in nodes_to_check:  # 2 reads, or 3 if supplementary read

            for f1, f2 in itertools.combinations_with_replacement(flags, 2):  # Check all edges within template
                u, v = (n, f1[0], f1[1]), (n, f2[0], f2[1])
                has_v, has_u = g.has_node(v), g.has_node(u)

                # Only add an edge if either u or v is missing, not if both (or neither) missing
                if has_u != has_v:  # XOR
                    new_edges.append((u, v, {"c": "w"}))

    # Regroup based on white edges (link together read-pairs)
    G.add_edges_from(new_edges)
    #click.echo("{} alignments processed into overlap graph".format(count - 1), err=True)
    return G, read_buffer


def partition_intervals(intervals):
    """Generates partitions rather.
    >>> partition_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [[('chr1', 1, 4), ('chr1', 2, 5)], [('chr2', 3, 5)]], True
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1], tup[2]))  # by chrom, start
    merged = []
    any_merged = False
    current_i_end = None  # Keep track of the bounds of the current interval, grows to capture nested intervals
    for higher in sorted_by_lower_bound:

        if not merged:
            merged.append([higher])
            current_i_end = higher[2]
            continue
        elif higher[0] != merged[-1][0][0]:  # Dont merge intervals on different chroms
            merged.append([higher])
            current_i_end = higher[2]
        else:
            lower_group = merged[-1]
            lower = lower_group[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= current_i_end: #lower[2]:
                merged[-1] = lower_group + [higher]
                if higher[2] > current_i_end:
                    current_i_end = higher[2]
                any_merged = True
            else:
                merged.append([higher])
                current_i_end = higher[2]

    return merged, any_merged


def merge_intervals(intervals):
    """
    >>> partition_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
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


def merge_intervals_fast(intervals):
    """Find the min and max of a set of intervals in single pass"""
    itr = iter(intervals)
    i_min, i_max = next(itr)
    for lower, higher in itr:
        if lower < i_min:
            i_min = lower
        if higher > i_max:
            i_max = higher
    return i_min, i_max


def make_merge_graph(chr1_itv, chr2_itv, sub_graph_id, links_d):
    ol_G = nx.Graph()
    ol_G.add_edges_from(zip(chr1_itv, chr2_itv))

    # Get partitions for each chrom region
    # Some intervals overlap, merge these into single intervals, any1==True when some intervals were merged
    parts, any1 = partition_intervals(chr1_itv + chr2_itv)

    # if any1:
    #     bloc_model_unconnected = blockmodel(ol_G, parts)
    #
    #
    #     for bm in nx.connected_component_subgraphs(bloc_model_unconnected, copy=False):
    #         bmnodes = len(bm.nodes())
    #         an = []
    #         for si in bm.nodes():
    #             an += list(si)
    #         echo("blocksub", set(an))
    #         if ('chr9', 7371526, 7371776) in an:
    #             echo(bm.nodes())
    #             echo("edges", bm.edges())
    #             quit()
    #
    #         sub_graph_id += 1
    #         for u, v in bm.edges():
    #             u = list(u)
    #
    #             v = list(v)
    #             if len(u) > 1:
    #                 n_links_u = sum([links_d[i] for i in u])
    #                 u_intervals = merge_intervals_fast([i[1:] for i in u])  # Create a final merge of block node
    #                 u = [u[0][0], u_intervals[0], u_intervals[1]]
    #             else:
    #                 n_links_u = links_d[u[0]]
    #                 u = list(u[0])
    #             if len(v) > 1:
    #                 n_links_v = sum([links_d[i] for i in v])
    #                 v_intervals = merge_intervals_fast([i[1:] for i in v])
    #                 v = [v[0][0], v_intervals[0], v_intervals[1]]
    #             else:
    #                 n_links_v = links_d[v[0]]
    #                 v = list(v[0])
    #             yield u, v, sub_graph_id, bmnodes, n_links_u, n_links_v
    #
    # else:
    for subg in nx.connected_component_subgraphs(ol_G):
        sub_graph_id += 1
        sug_nodes = len(subg.nodes())
        for u, v in subg.edges():
            yield list(u), list(v), sub_graph_id, sug_nodes, links_d[u], links_d[v]


def kd_cluster(args, infile, max_dist):

    click.echo("Finding reads", err=True)
    # regions = data_io.overlap_regions(args["include"])  # Todo?

    cluster_map = defaultdict(lambda: array.array('L'))

    # Cluster reads by genome pos, using a 2d plane, pos1 --> pos2
    # read_name --> Add all supplementary, primary only once (use mate pair info), no secondary
    primary_seen = set([])
    count = 0
    for r in infile:
        f = r.flag
        if f & 1548:  # fails quality checks, PCR duplicate, unmapped, mate unmapped
            continue
        if f & 256 or f & 2048:  # Marked as a not primary or supplementary
            continue
        if f & 2:  # Skip reads in proper pair
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

        if not f & 2048:
            primary_seen.add(name)
        count += 1

    del primary_seen  # Garbage collect

    sub_graph_id = 0
    for cmap_key, positions in cluster_map.iteritems():
        # For each chromosome pair
        parray = np.array(positions, dtype=np.uint32).reshape((-1, 2))
        chr1, chr2 = infile.get_reference_name(cmap_key[0]), infile.get_reference_name(cmap_key[1])
        if chr1 == chr2:
            click.echo("Calculating NN; links={} {}<-->{}".format(len(parray), chr1, chr2), err=True)

        # Use fast nearest neighbour lookup
        neigh = AnnoyIndex(2, metric="manhattan")
        for idx, item in enumerate(parray):
            neigh.add_item(idx, item)
        neigh.build(10)

        G = nx.Graph()

        numb_neigh = 6  # Add edges between n number of neighbours to a point
        for index, (x, y) in enumerate(parray):
            G.add_node(index)

            rng, dists = neigh.get_nns_by_vector([x, y], numb_neigh, include_distances=True)

            if dists[-1] < max_dist:
                numb_neigh -= 1
            else:
                numb_neigh += 1
            if numb_neigh < 3:
                numb_neigh = 3
            elif numb_neigh > 6:
                numb_neigh = 6

            for index_2, d2 in zip(rng, dists):
                if index == index_2:  # Skip self
                    continue
                if d2 > max_dist:  # Ensure within right distance
                    continue

                G.add_edge(index, index_2)

        # Need to merge any overlapping intervals, otherwise mutliple calls can be made
        # Only merge intervals locally
        cc = 0

        chr1_itv = []
        chr2_itv = []
        links_d = defaultdict(int)  # The number of links between nodes

        for sub in nx.connected_component_subgraphs(G):

            chr1_pos, chr2_pos = zip(*[parray[i] for i in sub.nodes()])  # Get the corresponding positions pos1, pos2
            links = len(sub.nodes())

            chr1_start, chr1_end = min(chr1_pos), max(chr1_pos)  # Bounding interval for pos1 and pos2
            chr2_start, chr2_end = min(chr2_pos), max(chr2_pos)

            cc += 1
            n1 = (chr1, chr1_start, chr1_end)  # Add interval node
            n2 = (chr2, chr2_start, chr2_end)
            chr1_itv.append(n1)
            chr2_itv.append(n2)
            links_d[n1] += links
            links_d[n2] += links

        # Merge all overlapping intervals for this chromosome pair
        for u, v, sub_graph_id, n_nodes, n_links_u, n_links_v in make_merge_graph(chr1_itv, chr2_itv, sub_graph_id, links_d):

            yield u, v, sub_graph_id, n_nodes, n_links_u, n_links_v


def region_connectivity(tempin):

    df = pd.read_csv(tempin, sep="\t", header=None, index_col=None, low_memory=False)

    # Only regions between pairs of chromosomes will be merged into subgraphs
    # Need to perform a genome wide merge of intervals, across all chromosome pairs

    df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "component", "n_nodes", "u_nodes", "v_nodes"]
    df = df[(df["chr1"].str.len() <= 5) & (df["chr2"].str.len() <= 5)]
    df = df.sort_values(["chr1", "start1"])

    links_d = {}
    itv1 = zip(df["chr1"], df["start1"], df["end1"])
    links_d.update(dict(zip(itv1, df["u_nodes"])))

    itv2 = zip(df["chr2"], df["start2"], df["end2"])
    links_d.update(dict(zip(itv2, df["v_nodes"])))

    for u, v, sub_graph_id, n_nodes, n_links_u, n_links_v in make_merge_graph(itv1, itv2, 0, links_d):
        yield u, v, sub_graph_id, n_nodes, n_links_u, n_links_v


def analyse_connectivity(f, tmps3):

    df = pd.read_csv(f, sep="\t", header=None, index_col=None, low_memory=False)

    # Only regions between pairs of chromosomes will be merged into subgraphs
    # Need to perform a genome wide merge of intervals, across all chromosome pairs

    df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "component", "n_nodes", "u_nodes", "v_nodes"]

    n1, n2 = zip(df["chr1"], df["start1"], df["end1"]), zip(df["chr2"], df["start2"], df["end2"])
    neighbour_counts = Counter(n1 + n2)

    con = []
    df["n1_connectivity"] = [neighbour_counts[i] for i in n1]
    df["n2_connectivity"] = [neighbour_counts[i] for i in n2]
    for u, v, unodes, vnodes in zip(n1, n2, df["u_nodes"], df["v_nodes"]):

        u_out = neighbour_counts[u]
        v_out = neighbour_counts[v]

        con.append(((float(u_out) / unodes) * (float(v_out) / vnodes)) / 2.)
    df["connectivity"] = con
    df = df[(df["u_nodes"] > 3) & (df["v_nodes"] > 3)]

    # Re-order by chromosome sorted order
    res = []
    for idx, r in df.iterrows():

        dr = dict(r)
        if r["chr1"] == r["chr2"]:

            if r["start1"] > r["start2"]:  # Swap sides
                start1, end1 = r["start1"], r["end1"]
                start2, end2 = r["start2"], r["end2"]
                unodes = r["u_nodes"]
                vnodes = r["v_nodes"]

                dr["start1"] = start2
                dr["end1"] = end2
                dr["start2"] = start1
                dr["end2"] = end1
                dr["u_nodes"] = vnodes
                dr["v_nodes"] = unodes

            start1, end1 = dr["start1"], dr["end1"] + 150
            start2, end2 = dr["start2"], dr["end2"] + 150

            # Check if intervals overlap
            if max(start1, start2) <= min(end1, end2):
                low_start = min(start1, start2)
                high_end = max(end1, end2)
                mid = (high_end - low_start) / 2
                start1 = low_start
                end1 = low_start + mid
                start2 = end1
                end2 = high_end
                dr["start1"] = start2
                dr["end1"] = end2
                dr["start2"] = start1
                dr["end2"] = end1
            res.append(dr)
        else:
            res.append(dr)
    df = pd.DataFrame.from_records(res)

    df = df.sort_values(["chr1", "start1", "chr2", "start2"])[['chr1','start1','end1','chr2','start2','end2','component','connectivity','n1_connectivity','n2_connectivity','n_nodes','u_nodes','v_nodes']]
    # Lenient filtering
    predict = [True if (i < 0.04 and j > 3 and k > 3 and l < 4 and m < 4) else False for i, j, k, l, m in zip(df["connectivity"], df["u_nodes"], df["v_nodes"], df["n1_connectivity"], df["n2_connectivity"])]

    df = df[predict]
    df.to_csv(tmps3, sep="\t", index=None)
    # return df


def get_region_depth(f, inputbam, temp_1, temp_2):

    df = pd.read_csv(f, sep="\t", index_col=None)
    intervals = []
    for c1, s1, e1, c2, s2, e2 in zip(df["chr1"], df["start1"], df["end1"], df["chr2"], df["start2"], df["end2"]):
        s1 = (int(s1) / 100) * 100
        s2 = (int(s2) / 100) * 100
        e1 = ((int(e1) / 100) * 100) + 100
        e2 = ((int(e2) / 100) * 100) + 100
        intervals += [[c1, s1, e1], [c2, s2, e2]]

    intervals = merge_intervals(intervals)
    with open(temp_1, "w") as depthf:
        for chrom, start, end in intervals:
            starts = np.arange(start, end, 100).astype(str)
            ends = np.arange(start + 100, end + 100, 100).astype(str)
            k = [chrom] * len(starts)
            lines = ["\t".join(item) + "\n" for item in zip(k, starts, ends)]
            depthf.writelines(lines)

    t0 = time.time()
    call("samtools bedcov {t1} {bam} > {t2}".format(bam=inputbam, t1=temp_1, t2=temp_2), shell=True)
    click.echo("Region coverage: {}s".format(int(time.time() - t0)), err=True)

    depth_d = {}
    with open(temp_2, "r") as depth_i:
        for line in depth_i:
            chrom, s, e, depth = line.strip().split("\t")
            depth_d[(chrom, int(s))] = float(depth) / 100.
    return depth_d


def filter_high_cov(region_windows, regions_depth, max_cov=150):

    new = []
    removed = 0
    for region_pair in region_windows:

        max_covs = []
        for chrom, s, e in region_pair:

            # Round start and end to get dict key
            start = (int(s) / 100) * 100
            end = ((int(e) / 100) * 100) + 100
            max_covs.append(max([regions_depth[(chrom, i)] for i in range(start, end, 100)]))

        if any([i > max_cov for i in max_covs]):
            removed += 1
            continue
        new.append(region_pair)
    click.echo("Dropped {} high coverage regions. Kept {}".format(removed, len(new)), err=True)
    return new


def get_regions_windows(args, unique_file_id, infile, max_dist):
    t0 = time.time()
    temps = ["{}/regionstemp{}.{}.csv".format(args["dest"], i, unique_file_id) for i in range(3)]

    with open(temps[0], "w") as tempout:
        for u, v, sub_id, n_nodes, nu, nv in kd_cluster(args, infile, max_dist):
            tempout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*tuple(u + v + [sub_id, n_nodes, nu, nv])))

    with open(temps[1], "w") as tempout:
        for u, v, sub_id, n_nodes, nu, nv in region_connectivity(tempin=temps[0]):
            tempout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*tuple(u + v + [sub_id, n_nodes, nu, nv])))

    analyse_connectivity(temps[1], temps[2])

    depth_temp_in = "{}/regionstempdin.{}.csv".format(args["dest"], unique_file_id)
    depth_temp_out = "{}/regionstempdout.{}.csv".format(args["dest"], unique_file_id)
    regions_depth = get_region_depth(temps[2], args["sv_aligns"], depth_temp_in, depth_temp_out)

    df = pd.read_csv(temps[2], sep="\t", index_col=None)

    region_windows = zip(zip(df["chr1"], df["start1"], df["end1"]), zip(df["chr2"], df["start2"], df["end2"]))

    region_windows = filter_high_cov(region_windows, regions_depth)

    # Remove temp files
    for f in temps + [depth_temp_in, depth_temp_out]:
        os.remove(f)
    click.echo("Regions prepared in {} h:m:s\n".format(str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)
    return region_windows


def worker(argument_set):
    args, infile_params, max_dist, region_pair = argument_set
    infile = pysam.AlignmentFile(*infile_params)
    G, read_buffer = construct_graph(args, infile, max_dist, input_windows=region_pair)

    potential_events = []
    for grp in nx.connected_component_subgraphs(G):  # Get large components, possibly multiple events

        reads = get_reads(infile, grp, max_dist, int(args['read_length']), read_buffer)

        bm = make_block_model(grp, args["insert_median"], args["insert_stdev"], args["read_length"])

        if len(bm.nodes()) == 0:
            continue
        bm = block_model_evidence(bm, grp)  # Annotate block model with evidence

        for event in caller.call_from_block_model(bm, grp, reads, infile, args["clip_length"],
                                                  args["insert_median"], args["insert_stdev"]):
            if event:
                potential_events.append(event)

    return potential_events


def check_pkl(val):
    try:
        pickle.dumps(val)
        return True
    except:
        return False


def cluster_reads(args):
    t0 = time.time()

    data_io.mk_dest(args["dest"])
    if args["dest"] is None:
        args["dest"] = "."

    kind = args["sv_aligns"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "rs"}

    click.echo("Input file is {}, (.{} format)".format(args["sv_aligns"], kind), err=True)
    infile = pysam.AlignmentFile(args["sv_aligns"], opts[kind])
    infile_params = (args["sv_aligns"], opts[kind])

    max_dist = int(args["insert_median"] + (5 * args["insert_stdev"]))  # > distance reads drop out of clustering scope
    click.echo("Maximum clustering distance is {}".format(max_dist), err=True)

    if args["svs_out"] == "-" or args["svs_out"] is None:
        click.echo("SVs output to stdout", err=True)
        outfile = sys.stdout
    else:
        click.echo("SVs output to {}".format(args["svs_out"]), err=True)
        outfile = open(args["svs_out"], "w")

    unique_file_id = str(uuid.uuid4())
    click.echo("Tempfile name {}".format(unique_file_id), err=True)

    region_windows = get_regions_windows(args, unique_file_id, infile, max_dist)

    regions = data_io.overlap_regions(args["include"])

    click.echo("Calling SVs with processors={}...".format(args["procs"]), err=True)

    args2 = {k: v for k, v in args.items() if check_pkl(v)}  # Picklable args

    with outfile:
        head = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "total_reads", "cipos95A", "cipos95B",
         "DP", "DApri", "DN", "NMpri", "SP", "EVsup", "DAsup",
         "NMsup", "maxASsup", "contig", "pe", "sup", "sc", "block_edge", "MAPQpri", "MAPQsup", "raw_reads_10kb", "kind"]

        outfile.write("\t".join(head) + "\n")

        jobs = []
        reg_counter = 0
        c = 0
        for region_pair in region_windows:

            jobs.append((args2, infile_params, max_dist, region_pair))

            if len(jobs) == 20 or reg_counter == len(region_windows) - 1:  # Todo add chunk size parameter
                # Process jobs
                if args["procs"] > 1:
                    pool = multiprocessing.Pool(args["procs"])
                    target = pool.map_async(worker, jobs)  # Asynchronous return
                    pool.close()
                    pool.join()
                    result = target.get()
                else:
                    result = map(worker, jobs)

                for potential_events in result:
                    if potential_events:
                        if potential_events:
                            # Keep one best event per region-pair (default mode)
                            best = sorted(potential_events, key=lambda x: sum([x["pe"], x["sup"]]))[-1]

                            # Collect coverage information
                            event_string = caller.get_raw_cov_information(best, infile, regions)
                            if event_string:
                                outfile.write(event_string)
                                c += 1
                jobs = []

            reg_counter += 1

        assert len(jobs) == 0

        click.echo("{} SV called in {} h:m:s\n".format(c, str(datetime.timedelta(seconds=int(time.time() - t0)))),
                   err=True)
