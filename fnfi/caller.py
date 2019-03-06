from collections import Counter, defaultdict
import click
import networkx as nx
import numpy as np
import assembler
import data_io


def echo(*args):
    click.echo(args, err=True)


def get_tuple(j):
    # Get breakpoint info from contig. Need (3 or 5 join, read names, chromosome, break point position, soft-clipped)
    return [(5 if j["left_clips"] > j["right_clips"] else 3,
            j["read_names"],
            j["bamrname"],
            j["ref_start"] if j["left_clips"] > j["right_clips"] else j["ref_end"],
            j["contig"][0].islower() or j["contig"][-1].islower())]


def guess_break_point(read, bam):

    # If the read is non-discordant and no clip is present, skip

    if read.flag & 2 and not any(i[0] == 4 or i[0] == 5 for i in read.cigartuples):
        return []

    # Sometimes a clip may be present, use this as a break point if available
    left = 0
    right = 0
    if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:  # Left soft or hard-clip
        left = read.cigartuples[0][1]
    if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
        right = read.cigartuples[-1][0]
    if left > 0 or right > 0:
        if left > right:
            return 5, read.qname, bam.get_reference_name(read.rname), read.pos, True
        else:
            return 3, read.qname, bam.get_reference_name(read.rname), read.reference_end, True
    else:
        # Breakpoint position is beyond the end of the last read
        if read.flag & 16:  # Read on reverse strand guess to the left
            p = read.pos
            t = 5
        else:  # Guess right
            p = read.reference_end
            t = 3
        return t, read.qname, bam.get_reference_name(read.rname), int(p), False


def process_node_set(node_set, all_reads, bam):

    break_points = {}
    seen = set([])
    for n, f, p in node_set:  # node name, flag, position
        if n in seen:
            continue
        seen.add(n)  # Seen template
        for (f, p) in all_reads[n].keys():
            aln = all_reads[n][(f, p)]
            guessed = guess_break_point(aln, bam)
            if len(guessed) > 0:
                break_points[(n, f, p)] = guessed

    return break_points


def pre_process_breakpoints(break_points_dict):

    # If number chroms > 2, reduce
    vals = break_points_dict.values()
    if not vals:
        return {}
    chroms = Counter([i[2] for i in vals if len(i) > 0])
    if len(chroms) == 0:
        return {}

    if len(chroms) > 2:
        c = [i[0] for i in sorted(chroms.items(), key=lambda x: x[1], reverse=True)][:2]
        chroms = {k: v for k, v in chroms.items() if k in c}

    # If minimum chrom counts < 0.1, reduce. Drop low coverage chromosomes
    if len(chroms) == 2:
        ci = list(chroms.items())
        total = float(sum(chroms.values()))
        if ci[0][1] / total < 0.05:  # Todo add parameter to list
            del chroms[ci[0][0]]  # Keep item 1
        elif ci[1][1] / total < 0.05:
            del chroms[ci[1][0]]  # Keep item 0

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

    #if len(clst) > 2:
    #    clst = sorted(clst, key=lambda x: len(x))[-2:]  # Choose largest 2

    #assert len(clst) <= 2
    return [dict(i) for i in clst]


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

        else:  # Try separate into smaller clusters
            clst = cluster_by_distance(break_points, t=25)
            if len(clst) == 2:
                c1, c2 = clst
            else:

                # Try separate using only supplementary
                supps = set([])
                for name, flg, pos in break_points_dict.keys():
                    if flg & 2048:
                        supps.add((name, bool(flg & 64)))
                        # if name == 'HISEQ2500-10:539:CAV68ANXX:7:2107:1342:6170':
                        #     echo((name, bool(flg & 64)))
                # Get only reads that are split i.e. primary + supplementary pairs
                break_points_dict = {k: v for k, v in break_points_dict.items() if (k[0], bool(k[1] & 64)) in supps}
                break_points = break_points_dict.items()
                # for item in break_points:
                #     if item[0] == 'HISEQ2500-10:539:CAV68ANXX:7:2107:1342:6170':
                #         echo("breakpoints", break_points)
                #         break

                if len(break_points) == 1:
                    c1 = break_points_dict

                elif len(break_points) == 2:
                    c1, c2 = dict([break_points[0]]), dict([break_points[1]])

                elif len(break_points) > 2:
                    clst = cluster_by_distance(break_points, t=25)
                    if len(clst) == 2:
                        c1, c2 = clst
                    elif len(clst) > 2:
                        c1, c2 = sorted(clst, key=lambda x: len(x))[-2:]  # Choose largest 2
                    else:  # Couldn't separate
                        c1, c2 = clst[0], {}  # Assume one cluster

                else:  # Couldn't separate by supplementary
                    c1, c2 = clst[0], {}  # Assume one cluster

    return c1, c2


def call_break_points(c1, c2):
    """
    Makes a call from a list of break points. Can take a list of break points, or one merged cluster of breaks.
    Outliers are dropped. Breakpoints are clustered using kmeans into sets
    :param c1: A 5 tuple (3 or 5 join, set([(read_name, flag)..]), chromosome, break point position,
                          soft-clipped)
    :param c2: Same as c1
    :return: Info dict containing a summary of the call
    """

    info = {}
    count = 0
    contributing_reads = set([])
    for grp in (c1, c2):

        if len(grp) == 0:
            continue  # When no c2 is found

        grp2 = [i for i in grp if i[4]]
        if len(grp2) > 0:
            grp = grp2

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
    # Todo other scores need adding, DAsup, NMsup. Bug in MAPQsup? sometimes events with supp have no MAPQsup
    if len(read_names) == 0:
        return {}

    data = {k: [] for k in ('DP', 'DApri', 'DN', 'NMpri', 'NP', 'DAsupp', 'NMsupp', 'maxASsupp', 'MAPQpri', 'MAPQsupp')}

    # roi = "HISEQ2500-10:541:CATW5ANXX:6:2113:11913:39334"
    # p = False
    for name, flag, pos in read_names:

        # if name == roi:
        #     p = True
        #     echo(str(all_reads[name].keys()))
        #     echo(str(all_reads[name][(flag, pos)]))

        read = all_reads[name][(flag, pos)]
        if not read.flag & 2048:
            idf = "pri"
            for item in ["DP", "DN"]:
                if read.has_tag(item):
                    data[item].append(float(read.get_tag(item)))
            data["MAPQpri"].append(read.mapq)
            if read.has_tag("NP") and float(read.get_tag("NP")) == 1:
                data["NP"].append(1)
        else:
            idf = "supp"
            if read.has_tag("AS"):
                data["maxASsupp"].append(read.get_tag("AS"))
            data["MAPQsupp"].append(read.mapq)

        for gk in ["DA", "MAPQ", "NM"]:
            key = gk + idf
            if read.has_tag(gk):
                data[key].append(float(read.get_tag(gk)))
    # if p:
    #     echo("data collected", data)
    averaged = {}
    for k, v in data.items():
        if k == "NP":
            averaged[k] = len(v)
        elif k == "maxASsupp":
            if len(v) > 0:
                averaged[k] = max(v)
            else:
                averaged[k] = 0
        elif len(v) > 0:
            averaged[k] = sum(v) / float(len(v))
        else:
            averaged[k] = 0
    # averaged["max_k_coverage"] = max_kmer(collected)
    # averaged["soft_clip_stacks"] = count_soft_clip_stacks(collected)

    return averaged


def breaks_from_one_side(node_set, reads, bam):
    tup = {}
    for nn, nf, np in node_set:
        tup[(nn, nf, np)] = guess_break_point(reads[nn][(nf, np)], bam)
    return pre_process_breakpoints(tup).values()  # Don't return whole dict


def single(parent_graph, bm, reads, bam, insert_size, insert_stdev, _debug_k):

    bmnodes = list(bm.nodes())

    assert len(bmnodes) == 1
    if _debug_k:
        if any(item in bmnodes[0] for item in _debug_k):
            echo("_debug_k in single", bmnodes[0])

    break_points = process_node_set(bmnodes[0], reads, bam)  # Dict, keyed by node

    if break_points is None or len(break_points) == 0:
        return

    break_points = pre_process_breakpoints(break_points)

    if not break_points:
        return

    dict_a, dict_b = separate_mixed(break_points, thresh=insert_size + insert_stdev)

    dict_a_subg = parent_graph.subgraph(dict_a.keys())
    dict_b_subg = parent_graph.subgraph(dict_b.keys())

    assembl1 = assembler.base_assemble(dict_a_subg, reads, bam)
    assembl2 = assembler.base_assemble(dict_b_subg, reads, bam)

    roi = "simulated_reads.0.10-id120_A_chr17:114709_B_chr1:3108649-16543"
    if roi in reads:
        echo("roi in reads (one_node)", reads)
        echo("assembly", assembl1)
        echo("dict_a", dict_a)
        echo("dict_b", dict_b)

    info, contrib_reads = call_break_points(dict_a.values(), dict_b.values())
    info["linked"] = 0  # Only 1 region so no linkage defined

    info["total_reads"] = len(contrib_reads)
    both = list(dict_a.keys()) + list(dict_b.keys())

    info['pe'] = len([1 for k, v in Counter([i[0] for i in both]).items() if v > 1])
    info['supp'] = len([1 for i in both if i[1] & 2048])
    info['sc'] = len([1 for name, flag, pos in both if "S" in reads[name][(flag, pos)].cigarstring])
    info["block_edge"] = 0

    info["contig"] = None
    info["contig_rev"] = None
    info["contig2"] = None
    info["contig2_rev"] = None
    if assembl1:
        if "contig" in assembl1 and (assembl1["left_clips"] > 0 or assembl1["right_clips"] > 0):
            info["contig"] = assembl1["contig"]
            info["contig_rev"] = assembl1["contig_rev"]
    if assembl2:
        if "contig" in assembl2 and (assembl2["left_clips"] > 0 or assembl2["right_clips"] > 0):
            info["contig2"] = assembl2["contig"]
            info["contig2_rev"] = assembl2["contig_rev"]
    #         if ass["left_clips"] > ass["right_clips"]:
    #             info["posA"] = ass["ref_start"]
    #         else:
    #             info["posA"] = ass["ref_end"]
    #         info["chrA"] = ass["bamrname"]

    info.update(score_reads(bmnodes[0], reads))

    if roi in reads:
        echo("info", info)

    return info


def one_edge(parent_graph, bm, reads, bam, clip_length, from_func=None):

    assert len(bm.nodes()) == 2
    ns = list(bm.nodes())

    # assembler.base_assemble(dict_a_subg, reads, bam)
    sub = parent_graph.subgraph(ns[0])
    sub2 = parent_graph.subgraph(ns[1])
    #     assemblies[node_set] = assembler.base_assemble(sub, reads, bam)

    as1 = assembler.base_assemble(sub, reads, bam)  #  assemblies[ns[0]]
    as2 = assembler.base_assemble(sub2, reads, bam)

    if as1 is None or len(as1) == 0:
        tuple_a = breaks_from_one_side(ns[0], reads, bam)
    else:
        tuple_a = get_tuple(as1)  # Tuple of breakpoint information

    if not tuple_a:
        return None

    if as2 is None or len(as2) == 0:
        tuple_b = breaks_from_one_side(ns[1], reads, bam)
    else:
        tuple_b = get_tuple(as2)
    if not tuple_b:
        return None

    info, contrib_reads = call_break_points(tuple_a, tuple_b)

    info["linked"] = 0
    if as1 is not None and len(as1) > 0 and as2 is not None and len(as2) > 0:
        as1, as2 = assembler.link_pair_of_assemblies(as1, as2, clip_length)
        if as1["linked"] == 1:
            info["linked"] = 1

    if "result" not in bm[ns[0]][ns[1]]:
        return None
    info.update(bm[ns[0]][ns[1]]["result"])
    info.update(score_reads(ns[0].union(ns[1]), reads))

    info["block_edge"] = 1
    info["contig"] = None
    info["contig_rev"] = None
    info["contig2"] = None
    info["contig2_rev"] = None

    if as1 is not None and "contig" in as1:
        info["contig"] = as1["contig"]
        info["contig_rev"] = as1["contig_rev"]

    if as2 is not None and "contig" in as2:
        info["contig2"] = as2["contig"]
        info["contig2_rev"] = as2["contig_rev"]

    roi = "simulated_reads.0.10-id120_A_chr17:114709_B_chr1:3108649-16543"
    if roi in reads:
        echo("roi in reads (one_edge)", from_func, len(reads), reads.keys())
        echo("as1", as1)
        echo("as2", as2)
        echo("info", info)

    return info


def multi(parent_graph, bm, reads, bam, clip_length):

    for u, v in bm.edges():
        rd = reads_from_bm_nodeset([u, v], reads)
        sub = bm.subgraph([u, v])
        yield one_edge(parent_graph, sub, rd, bam, clip_length, from_func="multi")


def reads_from_bm_nodeset(bm_nodeset, reads):
    qnames = set([])
    for nd in bm_nodeset:
        for item in nd:
            qnames.add(item[0])
    return {nm: reads[nm] for nm in qnames}


def call_from_block_model(bm_graph, parent_graph, reads, bam, clip_length, insert_size, insert_stdev, _debug_k=False):

    # Block model is not guaranteed to be connected
    for bm in nx.connected_component_subgraphs(bm_graph):

        rds = reads_from_bm_nodeset(bm.nodes(), reads)

        if len(bm.edges()) > 1:
            # Break apart connected
            for event in multi(parent_graph, bm, rds, bam, clip_length):
                yield event

        elif len(bm.nodes()) == 1:
            # Single isolated node
            yield single(parent_graph, bm, rds, bam, insert_size, insert_stdev, _debug_k)

        elif len(bm.edges()) == 1:
            # Easy case
            yield one_edge(parent_graph, bm, rds, bam, clip_length, from_func="bm")


def call_to_string(call_info):
    # Tab delimited string, in a bedpe style
    k = ["chrA", "posA", "chrB", "posB", "svtype", "join_type", "cipos95A", "cipos95B",
         "DP", "DApri", "DN", "NMpri", "NP", "DAsupp",
         "NMsupp", "maxASsupp", "contig", "contig2", "pe", "supp", "sc", "block_edge", "MAPQpri", "MAPQsupp", "raw_reads_10kb",
         "kind", "connectivity", "linked"]

    return "\t".join([str(call_info[ky]) if ky in call_info else str(None) for ky in k]) + "\n"


def calculate_coverage(chrom, start, end, region_depths):
    # Round start and end to get dict key
    start = (int(start) / 100) * 100
    if start < 0:
        start = 0
    end = ((int(end) / 100) * 100) + 100
    return [region_depths[(chrom, i)] for i in range(start, end, 100) if (chrom, i) in region_depths]


def get_raw_cov_information(r, regions, window_info, regions_depth):

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

    if not ar and not br:
        kind = "extra-regional"
        # Skip if regions have been provided; almost always false positives
        if regions is not None:  # Todo this throws away all genomic stuff! Too harsh
            return None

    switch = False
    if (br and not ar) or (not br and ar):
        kind = "hemi-regional"
        if not br and ar:
            switch = True

    if ar and br:

        if r["chrA"] == r["chrB"]:
            rA = list(regions[r["chrA"]].find_overlap(r["posA"], r["posA"] + 1))[0]
            rB = list(regions[r["chrB"]].find_overlap(r["posB"], r["posB"] + 1))[0]

            if rA[0] == rB[0] and rA[1] == rB[1]:
                kind = "intra_regional"
                # Put posA first
                if r["posA"] > r["posB"]:
                    switch = True

            else:
                kind = "inter-regional"
                if r["chrA"] != sorted([r["chrA"], r["chrB"]])[0]:
                    switch = True
        else:
            kind = "inter-regional"

    if switch:
        chrA, posA, cipos95A, contig2 = r["chrA"], r["posA"], r["cipos95A"], r["contig2"]
        r["chrA"] = r["chrB"]
        r["posA"] = r["posB"]
        r["cipos95A"] = r["cipos95B"]
        r["chrB"] = chrA
        r["posB"] = posA
        r["cipos95B"] = cipos95A
        r["contig2"] = r["contig"]
        r["contig"] = contig2

    reads_10kb = 0
    if kind == "hemi-regional":
        reads_10kb = sum(calculate_coverage(r["chrA"], r["posA"] - 10000, r["posA"] + 10000, regions_depth))

    r["kind"] = kind
    r["raw_reads_10kb"] = reads_10kb
    r["connectivity"] = window_info["connectivity"]
    return call_to_string(r)


