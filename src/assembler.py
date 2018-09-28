"""
A basic assembler. Takes an overlap graph and merges reads in-place in a pileup style. Different soft-clipped regions
are then overlapped of 'linked'.
"""

import networkx as nx
import itertools
import difflib
import heapq
import numpy as np
from collections import deque
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
    rd = [reads[n[0]][n[1]] for n in g.nodes()]

    G = nx.DiGraph()
    for r in rd:
        if r.seq is None:
            continue
        seq_pos = deque(list(zip(r.seq, r.get_reference_positions(full_length=True), r.query_qualities)))
        c_chrom = r.rname
        c_pos = r.pos
        strand = -1 if r.flag & 16 else 1
        pred = -1  # Start of sequence node
        v_first = None
        v_last = None
        for idx, (opp, length) in enumerate(r.cigartuples):

            if opp == 4:
                if idx == 0:  # Left clip
                    for j in range(length)[::-1]:
                        base, pos, qual = seq_pos.popleft()
                        offset = j + 1
                        u = pred
                        v = (base, c_chrom, c_pos, offset)
                        if not v_first:
                            v_first = v
                        update_edge(u, v, qual, G, "softclip_l", strand)
                        pred = v
                else:  # right clip
                    for j in range(length):
                        base, pos, qual = seq_pos.popleft()
                        offset = j + 1
                        u = pred
                        v = (base, c_chrom, c_pos, offset)
                        if j == length - 1:
                            v_last = v
                        update_edge(u, v, qual, G, "softclip_r", strand)
                        pred = v
            elif opp == 0 or opp == 7 or opp == 8:  # All match, match (=), mis-match (X)
                for j in range(length):
                    base, c_pos, qual = seq_pos.popleft()  # New c_pos defined
                    offset = 0
                    u = pred
                    v = (base, c_chrom, c_pos, offset)
                    if idx == 0:
                        v_first = v
                    elif idx == len(r.cigartuples) - 1 and j == length - 1:
                        v_last = v
                    update_edge(u, v, qual, G, "match", strand)
                    pred = v
            elif opp == 1:  # Insertion
                for j in range(length):
                    base, pos, qual = seq_pos.popleft()
                    offset = j + 1
                    u = pred
                    v = (base, c_chrom, c_pos, offset)
                    update_edge(u, v, qual, G, "insertion", strand)
                    pred = v
            elif opp == 5 or opp == 2:  # Hard clip or deletion
                continue
            elif opp == 3:
                break  # N

        # Add a start and end tag for each read, so reads contributing to contig sequence can be determined
        if v_first:
            G.node[v_first]["rid"].append((r.qname, r.flag))
        if v_last:
            G.node[v_last]["rid"].append((r.qname, r.flag))

    G.remove_node(-1)

    path = nx.algorithms.dag.dag_longest_path(G, weight="w")
    bases = []
    quals = []
    weights = []
    chroms = None
    sc_support_l, sc_support_r = 0, 0
    break_qual_l, break_qual_r = 0, 0
    strand_l, strand_r = 0, 0
    clipped = False
    ref_start, ref_end = 1e12, -1
    read_list = []

    for i in range(len(path)):
        ni = G.node[path[i]]
        weights.append(ni["w"])
        base, chrom, pos, offset = path[i]
        if ni["kind"] == "match":
            bases.append(base.upper())
        else:
            bases.append(base.lower())
        quals.append(254 if ni["w"] > 254 else ni["w"])  # Max qual is 254

        if ni["kind"] == "softclip_l":
            clipped = True
            if ni["n"] > sc_support_l:
                sc_support_l = ni["n"]
                break_qual_l = ni["w"] / float(ni["n"])
                strand_l = ni["strand"]
        if ni["kind"] == "softclip_r":
            clipped = True
            if ni["n"] > sc_support_r:
                sc_support_r = ni["n"]
                break_qual_r = ni["w"] / float(ni["n"])
                strand_r = ni["strand"]

        if pos < ref_start:
            ref_start = pos
        if pos > ref_end:
            ref_end = pos
        if not chroms:
            chroms = bam.get_reference_name(chrom)

        read_list += G.node[path[i]]["rid"]

    if clipped:
        res = {"bamrname": chroms,
                "qual_l": break_qual_l,
                "qual_r": break_qual_r,
                "base_quals": quals,
                "left_clips": sc_support_l,
                "right_clips": sc_support_r,
                "strand_l": strand_l,
                "strand_r": strand_r,
                "ref_start": ref_start,
                "ref_end": ref_end + 1,
                "read_names": set(read_list),
                "contig": "".join(bases),
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
        qual = int(i["qual_r"])
        if i["left_clips"] > i["right_clips"]:
            left_clipped = True
            negative_strand = True if i["strand_l"] < 0 else False
            sc_support = i["left_clips"]
            qual = int(i["qual_l"])
        seqs.append((i["contig"], sc_support, i["bamrname"], i["ref_start"], i["ref_end"] + 1, left_clipped,
                     negative_strand, qual))

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
