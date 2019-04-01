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


def check_contig_match(a, b, diffs=8, ol_length=21, supress_seq=True, return_int=False):

    query = StripedSmithWaterman(str(a), suppress_sequences=supress_seq)
    alignment = query(str(b))
    # echo(alignment)
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
        if return_int:
            return 1
        return (qs, qe, als, ale,
                alignment.cigar,
                alignment.aligned_query_sequence,
                alignment.aligned_target_sequence)


def get_upper_start_end(a):
    a_start, a_end = -1, 0
    for idx, l in enumerate(a):
        if l.isupper():
            if a_start == -1:
                a_start = idx
            if idx > a_end:
                a_end = idx
    return a_start, a_end + 1


def get_mark_result(res, insertion, a, b, sa, sb, a_start, a_end, b_start, b_end, b_rev=False):
    if insertion < 0:
        res["mark"] = "microh"
    else:
        res["mark"] = "ins"
    edit_dis = len([1 for i, j in zip(sa, sb) if i != j])
    res["mark_seq"] = sa
    res["mark_ed"] = edit_dis

    if insertion > 0:
        # Look for templated insertion
        a_align = a[a_start:a_end].replace("-", "")  # Remove any deletion markers
        b_align = b[b_start:b_end].replace("-", "")

        a_query = StripedSmithWaterman(sa)
        a_alignment = a_query(a_align)
        b_alignment = a_query(b_align)

        if a_alignment.optimal_alignment_score >= b_alignment.optimal_alignment_score:
            cont = "A"
            aln = a_alignment
        else:
            cont = "B"
            aln = b_alignment

        aqs = aln.aligned_query_sequence
        tqs = aln.aligned_target_sequence
        if tqs:
            edit_dis = len([1 for i, j in zip(aqs, tqs) if i.upper() != j])
            if b_rev and cont == "B":

                v = "contB_rev:pos={}:ed={}:align={}".format(aln.target_begin, edit_dis, aqs)
            else:
                v = "cont{}:pos={}:ed={}:align={}".format(cont, aln.target_begin, edit_dis, aqs)
            l = aln.target_end_optimal - aln.target_begin + 1

            res["templated_ins_info"] = v
            res["templated_ins_len"] = l

    return res


def get_microh_or_ins(aln_idx):
    qs, qe, als, ale, q_cigar, q_aln, t_aln = aln_idx
    a = q_aln  # Use actual alignment sequence - keep deletions and insertions in place
    b = t_aln

    a_start, a_end = get_upper_start_end(a)
    b_start, b_end = get_upper_start_end(b)

    # Check for overlap of gap
    res = {"mark": "blunt", "mark_seq": "", "mark_ed": "", "templated_ins_info": "", "templated_ins_len": ""}
    if a_start >= b_start:
        insertion = a_start - b_end
        if insertion != 0:
            v = slice(*sorted([a_start, b_end]))
            sa = a[v]
            sb = b[v]
            res = get_mark_result(res, insertion, a, b, sa, sb, a_start, a_end, b_start, b_end)

    else:
        insertion = b_start - a_end
        if insertion != 0:
            v = slice(*sorted([b_start, a_end]))
            sa = a[v]
            sb = b[v]
            res = get_mark_result(res, insertion, a, b, sa, sb, a_start, a_end, b_start, b_end)

    return res


def link_pair_of_assemblies(a, b, clip_length):

    # Safest way is to try forward and reverse complement
    m = check_contig_match(a["contig"], b["contig"], supress_seq=False)

    if m != 0:
        a["linked"] = 1
        h = get_microh_or_ins(m)

    else:

        m = check_contig_match(a["contig"], b["contig_rev"], supress_seq=False)
        if m != 0:
            a["linked"] = 1
            h = get_microh_or_ins(m)
        else:
            a["linked"] = 0
            h = {"mark": "None", "mark_seq": "", "mark_ed": "", "templated_ins_info": "",
                 "templated_ins_len": ""}
    a.update(h)
    return a, b


if __name__ == "__main__":

    a = "cgcccgcccaggtctgacctcagaagaactctgctccgccttcgcaatacccccgaagtctgtgcagagaagaacgcagctccgccctggcgatgctccctAACCCTAACCCTAACCCTAACCCTAACCCTTCCTCAGCCTCTCAACCTGCTTGGGTTACAGGTATGAGCCCGGGTGCCTAGCCAAACATTCCATTTTATATGTATATGCTAGGAAT"
    b = "ctgtgcagagaagaacgcagctccgccctggcgatgctccctAACCCTAACCCTAACCCTAACCCTAACCCTTCCTCAGCCTCTCAACCTGCTTGGGTTACAGGTATGAGCCCGGGTGCCTAGCCAAACATTCCATTTTATATGTATATGCTAGGAATGAATAATCT"
    c = "attcctagcatatacatataaaatggaatgtttggctaggcacccgggctcatacctgtaacccaagcaggttgagaggctgaggaagggttagggttagggttagggttagggttaggGAGCATCGCCAGGGCGGAGCTGCGTTCTTCTCTGCACAGACTTCGGGGGTATTGCGAAGGCGGAGCAGAGTTCTTCTGAGGTCAGACCTGGGCGGGCG"
    print check_contig_match(a, b, supress_seq=False)
    print get_microh_or_ins(check_contig_match(a, reverse_complement(c, len(c)), diffs=8, supress_seq=False))
    pass