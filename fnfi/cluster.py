from __future__ import absolute_import
import datetime
import itertools
import os
import multiprocessing
import time
import numpy as np
import random
from collections import defaultdict
import click
import networkx as nx
import pysam
import sys
import pickle
from . import graph_funcs, coverage, assembler, data_io, caller


def echo(*args):
    click.echo(args, err=True)


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
            if edge[2]['c'] in additional_nodes:  # If color in look_for list

                for n in (0, 1):  # Check both ends of the edge for uniqueness
                    if edge[n] not in seen:
                        if edge[n] in other_nodes:
                            del additional_nodes[edge[2]['c']]  # This evidence conflicts with another source
                            continue

                        starting_nodes.add(edge[n])
                        additional_nodes[edge[2]['c']].add(edge[n])

                        # Remove this type of evidence if over upper bound
                        if len(additional_nodes[edge[2]['c']]) > upper_bound:
                            del additional_nodes[edge[2]['c']]

        if len(starting_nodes) == 0:
            break
    return set().union(*additional_nodes.values())


def make_nuclei(g, _debug=[]):
    """
    Nuclei are primarily clusters of overlapping soft-clipped reads. If none are found in the graph, try using
    supplementary edges, or finally grey edges only.
    Additional nuclei are then created from yellow, or grey edges, only if they do not intersect with already defined
    nuclei. The reason for doing this is to help seperate overlapping events, with some events having soft-clips but
    other events having no soft-clips. If the non-soft-clipped events are not nucleated then they will be dropped.
    :param g:
    :param _debug: if True additional info is printed to stderr
    :return:
    """
    sub_grp = nx.Graph()
    look_for = ["y", "g"]
    edge_colors = defaultdict(list)
    for e in g.edges(data=True):
        edge_colors[e[2]['c']].append(e)

    edges = []
    if "b" in edge_colors:
        edges += edge_colors["b"]

    if "g" in edge_colors:
        edges += edge_colors["g"]

    if len(edges) == 0:
        if "y" in edge_colors:
            edges = edge_colors["y"]
            look_for.remove("y")

        elif "w" in edge_colors:  # Sometimes single read-pair mapping event, so only white edges available
            edges = edge_colors["w"]

    # Look for other unique connected components
    sub_grp.add_edges_from(edges)

    if _debug:
        echo("sub edge colors for nuclei", {k: len(v) for k, v in edge_colors.items()})
        echo("looking for", look_for)
        for item in _debug:
            if sub_grp.has_edge(*item):
                echo("subg has edge", item, sub_grp[item[0]][item[1]]["c"])
    return sub_grp


def make_block_model(g, insert_size, insert_stdev, read_length, _debug=[]):
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
    show_details = []
    _debug2 = []
    if _debug:

        for item in _debug:
            if g.has_node(item):
                item_name = item[0]
                echo("target node in g", item)
                echo("target has edges", g.edges(item, data=True))
                for u, v, d in g.edges(item, data=True):
                    if d["c"] == "w":
                        _debug2.append(u)
                        _debug2.append(v)
                for u, v, dta in g.edges(data=True):
                    if u[0] == item_name or v[0] == item_name:
                        # echo(u, v, dta)
                        show_details.append((u, v))
                echo("subg components", len(list(nx.connected_component_subgraphs(g))))

    # Make the inner block-model

    sub_grp = make_nuclei(g, show_details)

    sub_grp_cc = list(nx.connected_component_subgraphs(sub_grp))

    if _debug2:
        t = False
        for item in _debug2:
            if sub_grp.has_node(item):
                item_name = item[0]
                echo("target node in nuclei sub_grp", item)
                t = True
        if t:
            echo("showing edges with targets in nuclie sub_grp, with components", len(sub_grp_cc))

            for idx, comp in enumerate(sub_grp_cc):
                echo("component", idx)
                for u, v, dta in comp.edges(data=True):
                    if u[0] == item_name or v[0] == item_name:
                        echo(u, v, dta)

    # Drop any reads that are'nt in a connected component
    intersection = g.subgraph(sub_grp.nodes())

    inner_model = graph_funcs.blockmodel(intersection, partitions=[i.nodes() for i in sub_grp_cc])
    # return inner_model
    # Each node in the inner_model is actually a set of nodes.
    # For each node in the inner_model, explore the input graph for new edges
    # return inner_model
    # return inner_model

    all_node_sets = set(inner_model.nodes())

    # Try and find other nodes for each partition i.e. other supplementary or mate-pairs
    updated_partitons = []
    for node_set in inner_model.nodes():

        other_nodes = set([item for sublist in list(all_node_sets.difference(node_set)) for item in sublist])

        # Expected local nodes
        # Read clouds uncover variation in complex regions of the human genome. Bishara et al 2016.
        # max template size * estimated coverage / read_length; 2 just to increase the upper bound a bit
        local_upper_bound = ((insert_size + (2 * insert_stdev)) * float(2 + len(node_set))) / float(read_length)

        additional_nodes = explore_local(set(node_set), g, other_nodes, ["y", "g"], local_upper_bound*2)

        # Make sure additional nodes dont appear in other partitions. Can occasionally happen
        # if len(updated_partitons) > 0:
        #     for item in updated_partitons:
        #         additional_nodes = additional_nodes.difference(item)

        if len(additional_nodes) > 0:  # If any new nodes are found
            updated_partitons.append(set(node_set).union(additional_nodes))

    if len(updated_partitons) > 0:  # Re-partition the graph
        intersection = g.subgraph([item for sublist in updated_partitons for item in sublist])
        inner_model = graph_funcs.blockmodel(intersection, partitions=updated_partitons)

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
            # seen.add(neighbor_set)
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
            try:
                coords.append(dta['p'])
            except KeyError:
                echo("Warning: 'p' key not found on node")
                pass

    if len(coords) > 0:  # fetch any reads that were'nt in the buffer
        coords = itertools.groupby(sorted(set(coords), key=lambda x: (x[0], x[1])), key=lambda x: x[0])  # Groupby chrom

        # Make intervals of reads that are near one another
        # https://stackoverflow.com/questions/10016802/python-group-a-list-of-integer-with-nearest-values
        for chrom, crds in coords:
            d = [i[1] for i in crds]
            chrom = infile.get_reference_name(chrom)
            m = [[chrom, d[0] - (4 * rl), d[0] + (4 * rl)]]

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
    # Todo fix this
    # if c != len(sub_graph.nodes()):
    #     click.echo("WARNING: reads missing from collection - {} vs {}".format(c, len(sub_graph.nodes())), err=True)

    return rd


def filter_potential(input_events, tree):
    potential = []

    for i in input_events:
        if "posB" not in i:  # Skip events for which no posB was identified
            continue
        if "contig" not in i or i["contig"] == "":
            i["contig"] = None

        # Remove events for which both ends are in --include but not contig was found
        posA_intersects = data_io.intersecter(tree, i["chrA"], i["posA"], i["posA"] + 1)
        posB_intersects = data_io.intersecter(tree, i["chrB"], i["posB"], i["posB"] + 1)
        if (posA_intersects and posB_intersects) and i["contig"] is None:
            if i["NP"] == 0:
                continue

        # Remove events for which neither end is in --include (if --include provided)
        if tree:
            if not posA_intersects and not posB_intersects:
                continue

        potential.append(i)
    return potential


def merge_events(potential, max_dist, tree, seen, try_rev=False, pick_best=False):
    # Merging is not efficient currently.
    if len(potential) <= 1:
        return potential, seen

    G = nx.Graph()

    id_to_event_index = {}  # Mapping of event_id to index
    for idx in range(len(potential)):

        ei = potential[idx]
        i_id = ei["event_id"]
        id_to_event_index[i_id] = idx

        for jdx in range(len(potential)):

            ej = potential[jdx]
            j_id = ej["event_id"]
            if i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen:
                continue

            seen.add((i_id, j_id))

            # Check if events point to the same loci
            loci_similar = False
            loci_same = False
            if ei["chrA"] == ej["chrA"]:  # Try chrA matches chrA

                dist1 = abs(ei["posA"] - ej["posA"])
                if dist1 < max_dist:
                    if ei["chrB"] == ej["chrB"]:

                        dist2 = abs(ei["posB"] - ej["posB"])
                        if dist2 < max_dist:
                            loci_similar = True
                        if dist1 < 5 and dist2 < 5:
                            loci_same = True

            if not loci_similar:  # Try chrA matches chrB
                if ei["chrA"] == ej["chrB"]:
                    dist1 = abs(ei["posA"] - ej["posB"])
                    if dist1 < max_dist:
                        if ei["chrB"] == ej["chrA"]:
                            dist2 = abs(ei["posB"] - ej["posA"])
                            if dist2 < max_dist:
                                loci_similar = True
                            if dist1 < 5 and dist2 < 5:
                                loci_same = True

            if "contig" in ei and "contig" in ej:
                ci, cj, cj_rev = str(ei["contig"]), str(ej["contig"]), str(ej["contig_rev"])
                ci = "" if ci == "None" else ci
                cj = "" if cj == "None" else cj
                cj_rev = "" if cj_rev == "None" else cj_rev

                ci2, cj2, cj2_rev = str(ei["contig2"]), str(ej["contig2"]), str(ej["contig2_rev"])
                ci2 = "" if ci2 == "None" else ci2
                cj2 = "" if cj2 == "None" else cj2
                cj2_rev = "" if cj2_rev == "None" else cj2_rev

                any_ci = len(ci) > 0 or len(ci2) > 2
                any_cj = len(cj) > 0 or len(cj2) > 2
            else:
                continue

            if loci_similar:
                if any_ci and any_cj and loci_same:  # Try and match contigs

                    for cont1, cont2, cont2_rev in ((ci, cj, cj_rev), (ci2, cj2, cj2_rev),
                                                    (ci, cj2, cj2_rev), (ci2, cj, cj_rev)):

                        if len(cont1) > 0 and len(cont2) > 0:

                            # Each breakpoint can have a different assembly, only check for match if contigs overlap
                            idt = assembler.check_contig_match(ci, cj, diffs=15, return_int=True)
                            if idt == 1:
                                G.add_edge(i_id, j_id)
                                break

                            elif try_rev and cont2_rev:
                                if assembler.check_contig_match(cont1, cont2_rev, diffs=15, return_int=True) == 1:
                                    G.add_edge(i_id, j_id)
                                    break

                # No contigs to match, merge anyway
                elif not any_ci and not any_cj and loci_same:
                    G.add_edge(i_id, j_id)

                # Only merge loci if they are not both within --include regions. chrA:posA only needs checking
                elif not (data_io.intersecter(tree, ei["chrA"], ei["posA"], ei["posA"] + 1) and
                          data_io.intersecter(tree, ei["chrB"], ei["posB"], ei["posB"] + 1)):
                    G.add_edge(i_id, j_id)

    found = []
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item["event_id"]):
            found.append(item)

    for grp in nx.connected_component_subgraphs(G):

        c = [potential[id_to_event_index[n]] for n in grp.nodes()]

        best = sorted(c, key=lambda x: sum([x["pe"], x["supp"]]), reverse=True)
        w0 = best[0]["pe"] + best[0]["supp"]  # Weighting for base result

        if not pick_best:
            for k in range(1, len(best)):

                # Sum these
                for t in ["pe", "supp", "sc", "NP", "block_edge", "joinA", "joinB"]:
                    best[0][t] += best[k][t]

                if best[k]["maxASsupp"] > best[0]["maxASsupp"]:
                    best[0]["maxASsupp"] = best[k]["maxASsupp"]

                # Average these
                for t in ["DN", "MAPQsupp", "MAPQpri", "DApri", "DAsupp", "DP", "NMpri", "NMsupp"]:
                    w = best[k]["pe"] + best[k]["supp"]
                    denom = w0 + w
                    if denom == 0:
                        weighted_av = 0
                    else:
                        weighted_av = ((best[0][t] * w0) + (best[k][t] * w)) / denom
                    best[0][t] = weighted_av
                w0 = best[0]["pe"] + best[0]["supp"]
        found.append(best[0])

    return found, seen


def get_preliminary_events(G, read_buffer, infile, args, max_dist, regions, regions_depth,
                      approx_rl, _debug_k):

    # Go through edges in cluster map and make calls on block-model edges/nodes
    # Then merge all events to prevent duplication

    block_edge_events = []
    event_id = 0
    edges_merge_tested = set([])
    # _debug_k = {('HISEQ2500-10:539:CAV68ANXX:7:2216:13688:97541', 97, 26348614)}
    if _debug_k:
        echo("_debugk for get_prelim_events is ", _debug_k)
    t0 = time.time()
    for grp in nx.connected_component_subgraphs(G, copy=False):

        reads = get_reads(infile, grp, max_dist, approx_rl, read_buffer)

        bm = make_block_model(grp, args["insert_median"], args["insert_stdev"], approx_rl, _debug=_debug_k)

        for bmn in bm.nodes():

            if any(i in bmn for i in _debug_k):
                echo("read in bmn", len(bmn), len(bm.nodes()), len(bm.edges()))

        if len(bm.nodes()) == 0:
            continue

        # Annotate block model with evidence
        bm = graph_funcs.block_model_evidence(bm, grp)

        if _debug_k:
            for nds in bm.nodes():
                for item in _debug_k:
                    if item in nds:
                        echo(item, "in block model")

        potential_events = []
        for event in caller.call_from_block_model(bm, grp, reads, infile, args["clip_length"],
                                                  args["insert_median"], args["insert_stdev"], _debug_k):
            if event:
                event["event_id"] = event_id
                potential_events.append(event)
                event_id += 1

        potential_events = filter_potential(potential_events, regions)

        tested_edges = set([])

        potential_events, tested_edges = merge_events(potential_events, max_dist, regions, tested_edges,
                                                      try_rev=False, pick_best=True)

        edges_merge_tested = edges_merge_tested.union(tested_edges)
        block_edge_events += potential_events
    click.echo("Processed chunks {}s".format(round(time.time() - t0, 1)), err=True)
    # Perform a final merge across block nodes that haven't already been tested
    # Pick best=True prevents adding up of pe/supp, instead the best result is chosen
    t0 = time.time()
    merged, _ = merge_events(block_edge_events, max_dist, regions, edges_merge_tested, try_rev=True, pick_best=True)
    preliminaries = []
    if merged:
        for event in merged:
            # Collect coverage information
            event_dict = caller.get_raw_coverage_information(event, regions, regions_depth)
            if event_dict:
                preliminaries.append(event_dict)
    click.echo("Merged raw {}s".format(round(time.time() - t0, 1)), err=True)
    return preliminaries


def cluster_reads(args):
    t0 = time.time()
    np.random.seed(1)
    random.seed(1)
    try:
        model = pickle.load(open(args["model"], "rb"))
        click.echo("Model loaded from {}".format(args["model"]), err=True)
    except:
        model = None
        click.echo("No model loaded", err=True)

    data_io.mk_dest(args["dest"])
    if args["dest"] is None:
        args["dest"] = "."

    kind = args["sv_aligns"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "rs"}

    click.echo("Input file is {}, (.{} format). Processes={}".format(args["sv_aligns"], kind, args["procs"]), err=True)
    infile = pysam.AlignmentFile(args["sv_aligns"], opts[kind])
    sample_name = os.path.splitext(os.path.basename(args["sv_aligns"]))[0]

    if "insert_median" not in args and "I" in args:
        im, istd = map(float, args["I"].split(","))
        args["insert_median"] = im
        args["insert_stdev"] = istd

    max_dist = 2 * (int(args["insert_median"] + (5 * args["insert_stdev"])))  # reads drop from clustering scope
    click.echo("Maximum clustering distance is {}".format(max_dist), err=True)

    if args["svs_out"] == "-" or args["svs_out"] is None:
        click.echo("SVs output to stdout", err=True)
        outfile = sys.stdout
    else:
        click.echo("SVs output to {}".format(args["svs_out"]), err=True)
        outfile = open(args["svs_out"], "w")

    _debug_k = []
    regions = data_io.overlap_regions(args["include"])

    regions_merged, approx_rl, regions_depth = coverage.get_low_coverage_regions(infile, args["max_cov"], regions, args["include"])

    G, read_buffer = graph_funcs.construct_graph(infile, max_dist, regions, input_windows=regions_merged,
                                                buf_size=args["buffer_size"])

    preliminaries = get_preliminary_events(G, read_buffer, infile, args, max_dist, regions, regions_depth,
                      approx_rl, _debug_k)

    click.echo("Preliminaries found in {}s".format(round(time.time() - t0, 1)), err=True)

    classified_events_df = caller.calculate_prob_from_model(preliminaries, model)

    # Out order
    k = ["chrA", "posA", "chrB", "posB", "sample", "id", "kind", "svtype", "join_type", "cipos95A", "cipos95B",
         "DP", "DN", "DApri", "DAsupp",  "NMpri", "NMsupp", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "pe", "supp", "sc", "block_edge",
         "raw_reads_10kb",
          "linked", "contigA", "contigB", "mark", "mark_seq", "mark_ed", "templated_ins_info",
         "templated_ins_len", "Prob"]

    c = 0
    if classified_events_df is not None and len(classified_events_df) > 0:
        c = len(classified_events_df)
        classified_events_df["sample"] = [sample_name]*len(classified_events_df)
        classified_events_df["id"] = range(len(classified_events_df))
        classified_events_df = classified_events_df.rename(columns={"contig": "contigA", "contig2": "contigB"})
        classified_events_df[k].to_csv(outfile, index=False)

    click.echo("call-events {} complete, n={}, {} h:m:s".format(args["sv_aligns"],
                                                                c,
                                                                str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)
