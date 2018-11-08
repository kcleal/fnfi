import datetime
import itertools
import multiprocessing
import quicksect
import sys
import time
from collections import defaultdict, deque
from threading import Thread

import click
import networkx as nx
import numpy as np
import pysam
from pybedtools import BedTool

import caller
import data_io

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class Scoper:
    """Keeps track of which reads are in scope. Maximum distance depends on the template insert_median"""
    def __init__(self, max_dist, clip_length=20):
        self.max_dist = max_dist
        self.scope = deque([])
        self.clip_length = clip_length

    def update(self, current_read):

        if len(self.scope) == 0:
            self.scope.append(current_read)
            return

        elif len(self.scope) > 0 and current_read.rname != self.scope[-1].rname:
            self.scope = deque([current_read])  # Empty scope on new chromosome
            return

        current_pos = current_read.pos

        # Go through reads and check if they are in scope
        #out_of_scope = []
        while True:
            if len(self.scope) > 0:
                p = self.scope[0].pos
                if current_pos - p < self.max_dist:
                    break
                out = self.scope.popleft()
                #out_of_scope.append((out.qname, out.flag, out.pos))
                continue
            break
        self.scope.append(current_read)
        #return out_of_scope

    @staticmethod
    def overlap(start1, end1, start2, end2):
        return max(0, min(end1, end2) - max(start1, start2))

    def iterate(self):
        if len(self.scope) > 1:
            for i, t in enumerate(self.scope):
                if i == len(self.scope) - 1:
                    break  # Don't pair with self
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
                    if t.cigar[-1][0] == 4: # and t.cigar[-1][1] > self.clip_length:
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


def get_regions(reg_path):
    if not reg_path:
        return None
    else:
        f = BedTool(reg_path)
        regions = defaultdict(lambda: quicksect.IntervalTree())
        for item in f:
            regions[item.chrom].add(item.start, item.end)
        return regions


def align_match(r, t):

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


def make_nuclei(g, reads):
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
    part = list(map(frozenset, partitions))

    # Check for overlapping node partitions
    u = set()
    for p1, p2 in zip(part[:-1], part[1:]):
        u.update(p1)
        # if not u.isdisjoint(p2):  # Python 2.6 required
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


def make_block_model(g, insert_size, insert_stdev, read_length, reads):
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
    sub_grp = make_nuclei(g, reads)

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
                   "supp": supplementary}

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


def construct_graph(args, infile, max_dist, buf_size=1e6):  # todo add option for read buffer length
    click.echo("Constructing graph", err=True)
    regions = get_regions(args["include"])

    nodes = []
    edges = []

    all_flags = defaultdict(lambda: defaultdict(list))  # Linking clusters together rname: (flag, pos): (chrom, pos)

    # Make a buffer of reads to help prevent file lookups later
    read_buffer = dict()  # keys are (rname, flag, pos)
    read_index_buffer = dict()  # keys are int, values are (rname, flag, pos)

    scope = Scoper(max_dist=max_dist)

    # template = pysam.AlignmentFile(args["sv_bam"])
    # tt = pysam.AlignmentFile("{}.test.bam".format(args["sv_bam"][:-4]), "wb", template=template)
    count = 0
    buf_del_index = 0
    # unwritten = {}
    for r in infile:

        # written = False
        if count % 50000 == 0:
            if count != 0:
                click.echo("SV alignmnets processed {}".format(count), err=True)

        # Limit to regions
        if regions:
            rname = infile.get_reference_name(r.rname)
            if rname not in regions or len(regions[rname].search(r.pos, r.pos + 1)) == 0:
                rnext = infile.get_reference_name(r.rnext)
                if rnext not in regions or len(regions[rnext].search(r.pnext, r.pnext + 1)) == 0:
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

        nodes.append((n1, {"p": (r.rname, r.pos)}))
        all_flags[r.qname][(r.flag, r.pos)].append((r.rname, r.pos, r.flag))

        out_of_scope = scope.update(r)  # Add alignment to scope

        # if out_of_scope is not None and len(out_of_scope) > 0:
        #     for kk in out_of_scope:
        #         if kk in unwritten:
        #             del unwritten[kk]

        grey_added = 0

        for t, ol, sup_edge in scope.iterate():  # Other read, overlap current read, supplementary edge

            if not t:
                break

            n2 = (t.qname, t.flag, t.pos)
            all_flags[t.qname][(t.flag, t.pos)].append((t.rname, t.pos))

            # Four types of edges; black means definite match with overlapping soft-clips. Grey means similar
            # rearrangement with start and end co-ords on reference genome. Yellow means supplementary matches a primary
            # read; these edges need to be checked when both primary reads are available, change to black edge or delete
            # Finally white edge means a read1 to read2 edge
            if ol and not sup_edge:
                identity, prob_same = align_match(r, t)

                if prob_same > 0.01:  # Add a black edge
                    e = (n1, n2, {"c": "b"})  # "p": [(r.rname, r.pos), (t.rname, t.pos)],
                    # if not written:
                    #     tt.write(r)
                    #     written = True
                    #     if (t.qname, t.flag, t.pos) in unwritten:
                    #         tt.write(unwritten[(t.qname, t.flag, t.pos)])
                    #         del unwritten[(t.qname, t.flag, t.pos)]
                    edges.append(e)
                    continue

            # elif not ol and not sup_edge:
                # Make sure both are discordant, mapped to same chromosome
            if (r.flag & 2) or (t.flag & 2) or r.rnext == -1 or r.rnext != t.rnext:
                continue

            if abs(r.pnext - t.pnext) < max_dist:
                # Add a grey edge
                e = (n1, n2, {"c": "g"})
                grey_added += 1
                # if not written:
                #     tt.write(r)
                #     written = True
                #     if (t.qname, t.flag, t.pos) in unwritten:
                #         tt.write(unwritten[(t.qname, t.flag, t.pos)])
                #         del unwritten[(t.qname, t.flag, t.pos)]
                if grey_added > 60:
                    break  # Stop the graph getting unnecessarily dense
                edges.append(e)
                continue

            elif sup_edge:
                # Add a yellow edge
                e = (n1, n2, {"c": "y"})
                edges.append(e)
                continue

        # if not written:
        #     unwritten[(r.qname, r.flag, r.pos)] = r

    # template.close()
    # tt.close()
    # from subprocess import call
    # call("samtools sort -o {v}.test.srt.bam {v}.test.bam".format(v=args["sv_bam"][:-4]), shell=True)
    # call("samtools index {}.test.srt.bam".format(args["sv_bam"][:-4]), shell=True)

    # Add read-pair information to the graph, link regions together
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    new_edges = []

    for g in nx.connected_component_subgraphs(G, copy=True):

        # Add white edges between read-pairs which are NOT in the subgraph
        nodes_to_check = [(n, all_flags[n].keys()) for n, f, p in g.nodes()]

        for n, flags in nodes_to_check:  # 2 reads, or 3 if supplementary read

            for f1, f2 in itertools.combinations_with_replacement(flags, 2):  # Check all edges within template
                u, v = (n, f1[0], f1[1]), (n, f2[0], f2[1])
                has_v, has_u = g.has_node(v), g.has_node(u)

                # Only add an edge if either u or v is missing, not if both (or neither) missing
                if has_u != has_v:  # XOR
                    new_edges.append((u, v, {"c": "w"}))

    # Regroup based on white edges (link together read-pairs)
    G.add_edges_from(new_edges)
    click.echo("Read {} alignments into overlap graph".format(count - 1), err=True)
    return G, read_buffer


def setup_multi(args, worker, outf):
    cpus = args["procs"] if args["procs"] != 0 else multiprocessing.cpu_count()
    click.echo("fufi align runnning {} cpus".format(cpus), err=True)

    the_queue = multiprocessing.JoinableQueue(maxsize=cpus + 2)
    out_queue = multiprocessing.Queue()

    the_pool = multiprocessing.Pool(args["procs"] if args["procs"] != 0 else multiprocessing.cpu_count(),
                                    worker, (the_queue, out_queue,))

    def writer_thread(q, outsam):
        while True:
            aln = q.get()
            if aln == "Done":
                break
            elif aln == "Job failed":
                click.echo("job failed", err=True)
            elif len(aln) > 1:
                outsam.write(aln)

        click.echo("Writing done", err=True)

    writer = Thread(target=writer_thread, args=(out_queue, outf,))
    writer.setDaemon(True)
    writer.start()
    return the_queue, out_queue, writer


def worker(queue, out_queue):
    while True:
        job = queue.get(True)

        if job == "Done":
            queue.task_done()
            out_queue.put("Done")
            break

        else:
            big_string = ""
            for data_tuple in job:
                read_template = data_io.make_template(*data_tuple)
                process_template(read_template)
                if read_template['passed']:
                    outstring = data_io.to_output(read_template)
                    if outstring:
                        big_string += outstring
                    else:
                        pass  # No mappings

            if len(big_string) > 0:
                out_queue.put(big_string)
            else:
                click.echo("WARNING: no output from job.", err=True)
        queue.task_done()


def cluster_reads(args):
    t0 = time.time()
    infile = pysam.AlignmentFile(args["sv_bam"], "rb")
    max_dist = int(args["insert_median"] + (5 * args["insert_stdev"]))  # > distance reads drop out of clustering scope
    click.echo("Maximum clustering distance is {}".format(max_dist), err=True)

    G, read_buffer = construct_graph(args, infile, max_dist)

    if args["output"] == "-" or args["output"] is None:
        click.echo("Writing events to stdout", err=True)
        outfile = sys.stdout
    else:
        click.echo("Writing events to {}".format(args["output"]), err=True)
        outfile = open(args["output"], "w")

    regions = data_io.overlap_regions(args["include"])

    # Write events to output
    with outfile:
        head = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "total_reads", "cipos95A", "cipos95B",
         "DP", "DApri", "DN", "NMpri", "SP", "EVsup", "DAsup",
         "NMsup", "maxASsup", "contig", "pe", "sup", "sc", "block_edge", "MAPQpri", "MAPQsup", "raw_reads_10kb", "kind"]

        outfile.write("\t".join(head) + "\n")

        # if args["procs"] != 1:
        #     # Setup multiprocessing threads
        #     the_queue, out_queue, writer = setup_multi(args, worker, outfile)

        c = 0
        for grp in nx.connected_component_subgraphs(G):  # Get large components, possibly multiple events

            reads = get_reads(infile, grp, max_dist, int(args['read_length']), read_buffer)

            if args["procs"] == 1:

                bm = make_block_model(grp, args["insert_median"], args["insert_stdev"], args["read_length"], reads)
                if len(bm.nodes()) == 0:
                    continue
                bm = block_model_evidence(bm, grp)  # Annotate block model with evidence
                for event in caller.call_from_block_model(bm, grp, reads, infile, args["clip_length"],
                                                          args["insert_median"], args["insert_stdev"]):
                    if event:
                        # Collect coverage information
                        event_string = caller.get_raw_cov_information(event, infile, regions)
                        if event_string:
                            outfile.write(event_string)
                            c += 1
            else:
                pass
                # Send off for multiprocessing here. Need to find a way of pickling external objects
                # reads can be converted to a dict like object
                # infile needs the get_reference_name function
                # regions is a quicksetc object
                # dta = (infile, grp, args["insert_median"], args["insert_stdev"], args["read_length"], reads,
                #        args["clip_length"], regions)
                # the_queue.put(dta)

        # if args["procs"] != 0:
        #     the_queue.join()  # Wait for jobs to finish
        #     the_queue.put("Done")  # Send message to stop workers
        #     writer.join()  # Wait for writer to closing


    click.echo("cluster completed in {} h:m:s\n".format(str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)
