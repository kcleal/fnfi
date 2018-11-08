from collections import Counter, defaultdict

import click
import networkx as nx
import numpy as np

import assembler
import data_io


def echo(*args):
    click.echo(args, err=True)


def get_tuple(j):
    # Need (3 or 5 join, soft-clip length, chromosome, break point position)
    return [(5 if j["left_clips"] > j["right_clips"] else 3,
            j["read_names"],
            j["bamrname"],
            j["ref_start"] if j["left_clips"] > j["right_clips"] else j["ref_end"])]


def guess_break_point(read, bam, insert_size, insert_stdev):

    # If the read is non-discordant and no clip is present, skip
    if read.flag & 2 and not any(i[0] == 4 or i[0] == 5 for i in read.cigartuples):
        return []

    # Sometimes a clip may be present, use this as a break point if available
    cl = []  # Get the biggest clipped portion of the read-pair to find the break-point
    if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:  # Left clip
        cl.append((5, read.qname, bam.get_reference_name(read.rname), read.pos))
    if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
        cl.append((3, read.qname, bam.get_reference_name(read.rname), read.reference_end))
    if len(cl) > 0:
        return sorted(cl, key=lambda x: x[1])[-1]
    else:
        # Breakpoint position is beyond the end of the last read
        if read.flag & 16:  # Read on reverse strand guess to the left
            p = read.pos - (insert_size / 2) + insert_stdev
            t = 5
        else:  # Guess right
            p = read.reference_end + (insert_size / 2) - insert_stdev
            t = 3
        return t, read.qname, bam.get_reference_name(read.rname), int(p)


def process_node_set(node_set, all_reads, bam, insert_size, insert_stdev):

    break_points = {}
    seen = set([])
    # for u, v, d in edges:
    for n, f, p in node_set:  # node name, flag
        if n in seen:
            continue
        seen.add(n)  # Seen template
        for (f, p) in all_reads[n].keys():
            aln = all_reads[n][(f, p)]
            guessed = guess_break_point(aln, bam, insert_size, insert_stdev)
            if len(guessed) > 0:
                break_points[(n, f, p)] = guessed

    return break_points


def pre_process_breakpoints(break_points_dict):

    # If number chroms > 2, reduce
    chroms = Counter([i[2] for i in break_points_dict.values()])
    if len(chroms) > 2:
        c = [i[0] for i in sorted(chroms.items(), key=lambda x: x[1], reverse=True)][:2]
        chroms = {k: v for k, v in chroms.items() if k in c}

    # If minimum chrom counts < 0.1, reduce. Drop low coverage chromosomes
    if len(chroms) == 2:
        ci = list(chroms.items())
        total = float(sum(chroms.values()))
        if ci[0][1] / total < 0.05:
            del chroms[ci[0][1]]  # Keep item 1
        elif ci[1][1] / total < 0.05:
            del chroms[ci[1][1]]  # Keep item 0

    return {k: v for k, v in break_points_dict.items() if v[2] in chroms}


def cluster_by_distance(bpt, t):
    """
    Naively merge breakpoints into clusters based on a distance threshold.
    :param bpt: key-value pair of (read_name, flag) - breakpoint info
    :param t: clustering distance threshold
    :return: a list of clusters, each cluster is a dict with items of key=(read_name, flag), value=breakpoint info
    """
    clst = []
    current = []
    for node_name, i in sorted(bpt, key=lambda x: (x[1][2], x[1][3])):  # Sort by chrom and pos
        if len(current) == 0:
            current.append((node_name, i))
            continue
        chrom, pos = i[2], i[3]
        last_chrom, last_pos = current[-1][1][2], current[-1][1][3]
        if chrom == last_chrom and abs(last_pos - pos) < t:
            current.append((node_name, i))
        else:
            clst.append(current)
            current = [(node_name, i)]

    if len(current) > 0:
        clst.append(current)

    if len(clst) > 2:
        clst = sorted(clst, key=lambda x: len(x))[-2:]  # Choose largest 2

    assert len(clst) <= 2
    return list(map(dict, clst))


def separate_mixed(break_points_dict, thresh=500):

    break_points = list(break_points_dict.items())
    c1, c2 = {}, {}
    if len(break_points) == 1:
        c1 = break_points_dict

    elif len(break_points) == 2:
        c1, c2 = dict([break_points[0]]), dict([break_points[1]])

    elif len(break_points) > 2:
        clst = cluster_by_distance(break_points, t=thresh)
        if len(clst) == 2:
            c1, c2 = clst

        else:  # Try seperate into smaller clusters
            clst = cluster_by_distance(break_points, t=25)
            if len(clst) == 2:
                c1, c2 = clst
            else:
                c1, c2 = clst[0], {}  # Assume one cluster

            # X = [[i[0], i[3]] for i in break_points]
            # k_means = KMeans(init='k-means++', n_clusters=2, n_init=5, max_iter=20)
            # labels = k_means.fit_predict(X)
            # echo(labels)
            # grp = {0: [], 1: []}
            # for l, bp in zip(labels, break_points):
            #     grp[l].append(bp)
            # c1, c2 = grp[0], grp[1]

    return c1, c2


def call_break_points(c1, c2):
    """
    Makes a call from a list of break points. Can take a list of break points, or one merged cluster of breaks.
    Outliers are dropped. Breakpoints are clustered using kmeans into sets
    :param break_points: A 4 tuple (3 or 5 join, set([(read_name, flag)..]), chromosome, break point position)
    :param thresh: the distance threshold to determine if clustered
    :return: Info dict containing a summary of the call
    """

    info = {}
    count = 0
    contributing_reads = set([])
    for grp in (c1, c2):
        if len(grp) == 0:
            continue  # When no c2 is found
        for i in grp:
            contributing_reads = contributing_reads.union(i[1])

        chrom = Counter([i[2] for i in grp]).most_common()[0][0]
        grp = [i for i in grp if i[2] == chrom]
        sc_side = Counter([i[0] for i in grp]).most_common()[0][0]
        bp = [i[3] for i in grp]
        mean_pos = int(np.mean(bp))
        mean_95 = abs(int(np.percentile(bp, [97.5])) - mean_pos)
        if count == 0:
            side = "A"
        else:
            side = "B"
        info["chr" + side] = chrom
        info["pos" + side] = mean_pos
        info["cipos95" + side] = mean_95
        info["join" + side] = sc_side
        count += 1

    if "chrB" in info and info["chrA"] == info["chrB"]:
        if info["joinA"] == info["joinB"]:
            info["svtype"] = "INV"
            info["join_type"] = str(info["joinA"]) + "to" + str(info["joinB"])
        else:
            if info["posA"] <= info["posB"]:
                x = "joinA"
                y = "joinB"
            else:
                x = "joinB"
                y = "joinA"
            if info[x] == 3 and info[y] == 5:
                info["svtype"] = "DEL"
                info["join_type"] = "3to5"
            elif info[x] == 5 and info[y] == 3:
                info["svtype"] = "DUP"
                info["join_type"] = "5to3"

    elif "chrB" in info:
        info["svtype"] = "TRA"
        info["join_type"] = str(info["joinA"]) + "to" + str(info["joinB"])
    else:
        info["svtype"] = "BND"
        if "joinA" in info:
            info["join_type"] = str(info["joinA"]) + "to?"
        else:
            info["join_type"] = "?to?"

    return info, contributing_reads


def max_kmer(reads, k=27):
    kmers = defaultdict(int)
    for s in (i.seq for i in reads if not i.flag & 2048):  # Skip supplementary
        if s:
            for j in [s[i:i+k] for i in range(len(s) - k)]:
                kmers[j] += 1
    if len(kmers) > 0:
        return max(kmers.values())
    return 0


def merge_intervals(intervals):
    # thanks https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def count_soft_clip_stacks(reads):
    blocks = []
    for r in reads:  # (i for i in reads if not i.flag & 2048):  # Skip supp
        cigar = r.cigartuples
        if cigar:
            if (cigar[0][0] == 4 and cigar[0][1] > 20) or (cigar[-1][0] == 4 and cigar[-1][1] > 20):
                blocks.append((r.pos, r.reference_end))
    mi = merge_intervals(blocks)
    return len(mi)


def score_reads(read_names, all_reads):

    if len(read_names) == 0:
        return {}

    collected = []
    data = defaultdict(list)
    for name, flag, pos in read_names:
        # for flag in all_reads[name]:  # Could score all alignments from a read-template?
        read = all_reads[name][(flag, pos)]
        collected.append(read)
        if not read.flag & 2048:
            idf = "pri"
            _ = data["SP"]
            if read.has_tag("SP") and int(read.has_tag("SP")) == 1:
                data["SP"].append(1)
            data["MAPQpri"].append(read.mapq)
        else:
            idf = "sup"
            _ = data["EVsup"]
            if read.has_tag("EV"):
                data["EVsup"].append(read.get_tag("EV"))
            _ = data["maxASsup"]
            if read.has_tag("AS"):
                data["maxASsup"].append(read.get_tag("AS"))
            data["MAPQsup"].append(read.mapq)

        for item in ["DA", "NM"]:
            _ = data[item]  # Make sure item exists in default dict, otherwise problems arise with header
            if read.has_tag(item):
                data[item + idf].append(float(read.get_tag(item)))
            else:
                echo("Missing tag {}".format(item))
        if idf == "pri":
            for item in ["DP", "DN"]:
                _ = data[item]  # Make sure item exists in default dict, otherwise problems arise with header
                if read.has_tag(item):
                    data[item].append(float(read.get_tag(item)))
                else:
                    echo("Missing tag {}".format(item))

    averaged = {}
    for k, v in data.items():
        if k == "EVsup":
            if len(v) > 0:
                averaged[k] = min(v)
            else:
                averaged[k] = None
        elif k == "SP":
            averaged[k] = len(v)
        elif k == "maxASsup":
            averaged[k] = max(v)
        elif len(v) > 0:
            averaged[k] = sum(v) / float(len(v))
        else:
            averaged[k] = 0
    #averaged["max_k_coverage"] = max_kmer(collected)
    #averaged["soft_clip_stacks"] = count_soft_clip_stacks(collected)

    return averaged


def breaks_from_one_side(node_set, reads, bam, insert_size, insert_stdev):
    tup = {}
    for nn, nf, np in node_set:
        tup[(nn, nf, np)] = guess_break_point(reads[nn][(nf, np)], bam, insert_size, insert_stdev)
    return pre_process_breakpoints(tup).values()  # Don't return whole dict


def single(bm, reads, bam, insert_size, insert_stdev):
    # Todo single should use an assembly if available
    bmnodes = list(bm.nodes())
    assert len(bmnodes) == 1

    break_points = process_node_set(bmnodes[0], reads, bam, insert_size, insert_stdev)  # Dict, keyed by node
    if break_points is None or len(break_points) == 0:
        return

    break_points = pre_process_breakpoints(break_points)
    dict_a, dict_b = separate_mixed(break_points)
    info, contrib_reads = call_break_points(dict_a.values(), dict_b.values())

    info["total_reads"] = len(contrib_reads)
    both = list(dict_a.keys()) + list(dict_b.keys())

    info['pe'] = len([1 for k, v in Counter([i[0] for i in both]).items() if v > 1])
    info['sup'] = len([1 for i in both if i[1] & 2048])
    info['sc'] = len([1 for name, flag, pos in both if "S" in reads[name][(flag, pos)].cigarstring])
    info["block_edge"] = 0

    info.update(score_reads(bmnodes[0], reads))

    return info


def one_edge(bm, reads, bam, assemblies, clip_length, insert_size, insert_stdev):

    assert len(bm.nodes()) == 2
    ns = list(bm.nodes())
    as1 = assemblies[ns[0]]
    as2 = assemblies[ns[1]]
    # roi = "simulated_reads.0.10-id277_A_chr21:46699632_B_chr17:12568030-36717"
    # if roi in reads:
    #     echo("as1", as1)
    #     echo("as2", as2)

    if as1 is None or len(as1) == 0:
        tuple_a = breaks_from_one_side(ns[0], reads, bam, insert_size, insert_stdev)
    else:
        tuple_a = get_tuple(as1)  # Tuple of breakpoint information

    if as2 is None or len(as2) == 0:
        tuple_b = breaks_from_one_side(ns[1], reads, bam, insert_size, insert_stdev)
    else:
        tuple_b = get_tuple(as2)

    if as1 is not None and len(as1) > 0 and as2 is not None and len(as2) > 0:
        as1, as2 = assembler.link_pair_of_assemblies(as1, as2, clip_length)

    info, contrib_reads = call_break_points(tuple_a, tuple_b)

    # if roi in reads:
    #     echo("info", info)
    #     echo("r", "result" in bm[ns[0]][ns[1]])

    if "result" not in bm[ns[0]][ns[1]]:
        return None
    info.update(bm[ns[0]][ns[1]]["result"])
    info.update(score_reads(ns[0].union(ns[1]), reads))
    info["block_edge"] = 1

    if as1 is not None and "contig" in as1:
        info["contig"] = as1["contig"]
    elif as2 is not None and "contig" in as2:
        info["contig"] = as2["contig"]
    else:
        info["contig"] = ""

    # if roi in reads:
    #     echo(info)

    return info


def multi(bm, reads, bam, assemblies, clip_length, insert_size, insert_stdev):
    # Score edges
    edge_scores = []
    for u, v, data in bm.edges(data=True):
        score = data["weight"]
        if "result" in data:
            score = sum(data["result"].values())
        edge_scores.append((u, v, score))

    edge_scores = sorted(edge_scores, key=lambda k: k[2], reverse=True)
    seen = set([])
    for u, v, scr in edge_scores:
        if u in seen or v in seen:
            continue
        # seen.add(u)
        # seen.add(v)
        sub = bm.subgraph([u, v])
        yield one_edge(sub, reads, bam, assemblies, clip_length, insert_size, insert_stdev)


def call_from_block_model(bm_graph, parent_graph, reads, bam, clip_length, insert_size, insert_stdev):

    # Block model is not guaranteed to be connected
    for bm in nx.connected_component_subgraphs(bm_graph):

        # Try and assemble nodes in the block model
        assemblies = {}
        for node_set in bm.nodes():

            sub = parent_graph.subgraph(node_set)
            assemblies[node_set] = assembler.base_assemble(sub, reads, bam)

            # if "simulated_reads.0.10-id277_A_chr21:46699632_B_chr17:12568030-36717" in reads:
            #     echo("assem", assemblies)
            #     echo(sub.nodes())

            # continue
            # if any(i[2]["c"] == "b" or i[2]["c"] == "y" for i in sub.edges(data=True)):
            #     # assemble reads if any overlapping soft-clips
            #     assemblies[node_set] = assembler.base_assemble(sub, reads, bam)
            #
            # else:
            #     assemblies[node_set] = {}
            #     if "simulated_reads.0.10-id270_A_chr21:46697027_B_chr17:6458938-35841" in reads:
            #         click.echo("hi", err=True)
            #         quit()

        if len(bm.edges) > 1:
            # Break apart connected
            for event in multi(bm, reads, bam, assemblies, clip_length, insert_size, insert_stdev):
                yield event

        if len(bm.nodes()) == 1:
            # Single isolated node
            yield single(bm, reads, bam, insert_size, insert_stdev)

        if len(bm.edges) == 1:
            # Easy case
            yield one_edge(bm, reads, bam, assemblies, clip_length, insert_size, insert_stdev)


def call_to_string(call_info):
    # Tab delimited string, in a bedpe style
    k = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "total_reads", "cipos95A", "cipos95B",
         "DP", "DApri", "DN", "NMpri", "SP", "EVsup", "DAsup",
         "NMsup", "maxASsup", "contig", "pe", "supp", "sc", "block_edge", "MAPQpri", "MAPQsup", "raw_reads_10kb", "kind"]

    return "\t".join([str(call_info[ky]) if ky in call_info else str(None) for ky in k]) + "\n"


def get_raw_cov_information(r, raw_bam, regions):

    # Check if side A in regions
    ar = False
    if data_io.intersecter(regions, r["chrA"], r["posA"], r["posA"] + 1):
        ar = True

    if "chrB" not in r:  # todo fix this
        return None
    br = False
    if data_io.intersecter(regions, r["chrB"], r["posB"], r["posB"] + 1):
        br = True

    # Put non-region first
    kind = None
    if (br and not ar) or (not br and ar):
        kind = "hemi-regional"
        if not br and ar:
            chrA, posA, cipos95A = r["chrA"], r["posA"], r["cipos95A"]
            r["chrA"] = r["chrB"]
            r["posA"] = r["posB"]
            r["cipos95A"] = r["cipos95B"]
            r["chrB"] = chrA
            r["posB"] = posA
            r["cipos95B"] = cipos95A

    if not ar and not br:
        kind = "extra-regional"
    if ar and br:
        if r["chrA"] == r["chrB"]:
            # rA = regions[r["chrA"]].find_overlap(quicksect.Interval(r["posA"], r["posA"]+1))[0]
            # rB = regions[r["chrB"]].find_overlap(quicksect.Interval(r["posB"], r["posB"] + 1))[0]
            rA = list(regions[r["chrA"]].find_overlap(r["posA"], r["posA"] + 1))[0]
            rB = list(regions[r["chrB"]].find_overlap(r["posB"], r["posB"] + 1))[0]

            if rA[0] == rB[0] and rA[1] == rB[1]:
                kind = "intra_regional"
            else:
                kind = "inter-regional"
        else:
            kind = "inter-regional"

    reads_10kb = None
    if kind == "hemi-regional":

        # Get read coverage in raw bam file +/- 10kb around A-side
        start = r["posA"] - 10000
        reads_10kb = len(list(raw_bam.fetch(r["chrA"], 0 if start < 0 else start, r["posA"] + 10000)))

    r["kind"] = kind
    r["raw_reads_10kb"] = reads_10kb

    return call_to_string(r)


