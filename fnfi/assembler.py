"""
A basic assembler. Takes an overlap graph and merges reads in-place in a pileup style. Different soft-clipped regions
are then overlapped of 'linked'.
"""

import networkx as nx
import numpy as np
from collections import defaultdict
import click
from skbio.alignment import StripedSmithWaterman
from c_io_funcs import reverse_complement
from graph_funcs import find_cycle, dag_longest_path

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
                o += range(1, length)

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


def base_assemble(g, reads, bam, iid=0):
    """
    Assembles reads that have overlaps. Uses alignment positions to determine contig construction
    :param g: The overlap graph
    :param reads: Dict of read_name: flag: alignment
    :param bam: Original bam for header access
    :param iid: Unique ID for the event
    :return: Returns None if no soft-clipped portion of the cluster was assembled, otherwise a result dict is returned
    """
    # Note supplementary are included in assembly; helps link regions
    # Get reads of interest

    rd = [reads[n[0]][(n[1], n[2])] for n in g.nodes() if n[0] in reads]

    #
    # rnames = set([r.qname for r in rd])
    # roi = "simulated_reads.intra.0.2-id2_A_chr21:46696380_B_chr21:46697291-69"
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

    path = None
    try:
        path = dag_longest_path(G, weight="weight")
    except:
        cy = find_cycle(G, orientation="original")
        echo("Cycle in graph", cy)
        for r in rd:
            echo(bam.get_reference_name(r.rname))
            echo(str(r).split("\t"))

    if not path:
        return None  # No assembly

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
           "ref_end": matches[-1],
           "read_names": set(read_names),
           "contig": bases,
           "contig_rev": reverse_complement(bases, len(bases)),  # .upper()
           "id": iid}

    return res


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


def check_contig_match(a, b, diffs=8, ol_length=21):

    query = StripedSmithWaterman(str(a), suppress_sequences=True)
    alignment = query(str(b))

    qs, qe = alignment.query_begin, alignment.query_end
    als, ale = alignment.target_begin, alignment.target_end_optimal

    # Find the length of any unaligned overhangs
    extent_left = min((qs, als))
    extent_right = min((len(a) - qe, len(b) - ale))
    total_overhangs = extent_left + extent_right
    aln_s = alignment.optimal_alignment_score
    expected = (qe - qs) * 2  # +2 is the score for a match

    if expected < 2 * ol_length:  # Match score * Minimum clip length
        return 0
    diff = expected - aln_s + total_overhangs  # old diff thresh = 8
    # diff = (expected - aln_s - total_overhangs) / expected
    if diff > diffs:  # e.g. 2 mis-matches + 2 unaligned overhanging bits
        return 0
    else:
        return 1


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

        if i is None:
            a["linked"] = 0
            return a

        seqs.append((i["contig"], sc_support, i["bamrname"], i["ref_start"], i["ref_end"] + 1, left_clipped,
                     negative_strand, i["contig_rev"]))

    ainfo, binfo = seqs
    aseq, bseq, aseq_rev, bseq_rev = ainfo[0], binfo[0], ainfo[7], binfo[7]
    if ainfo[2] == binfo[2]:  # If seqs are on the same chrom
        if ainfo[5] == binfo[5]:  # If clips are on same side, rev comp one of them
            bseq = bseq_rev  # rev_comp(bseq)
    else:
        if ainfo[6]:  # negative strand
            aseq = aseq_rev  # rev_comp(aseq)
        if binfo[6]:
            bseq = bseq_rev  # rev_comp(bseq)

    if check_contig_match(aseq, bseq) == 1:  # .upper()
        a["linked"] = 1
    else:
        a["linked"] = 0

    return a, b
