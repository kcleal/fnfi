"""
Find the best path through a collection of alignments
Works with paried-end reads rather or single contigs. Borrows the pairing heuristic from bwa.
"""

import sys
import numpy as np
import os
import click
import align_path_c


def read_orientations(t, r1_len, r2_len):
    """
    # This code prooves that the orientation of the read dowesnt matter, only the order of the reads:
    # Permute all orientations, [AB|CD] + AB|DC, BA|DC, BA|CD
    ab, cd = [t[t[:, 7] == 1], t[t[:, 7] == 2]]

    # Make AB reversed
    ba = ab.copy()
    s, e = ba[:, 2].copy(), ba[:, 3].copy()
    ba[:, 2] = r1_l - e
    ba[:, 3] = r1_l - s
    ba = np.flipud(ba)

    # Make CD reversed
    dc = cd.copy()
    dc[:, 2:4] -= r1_l
    s, e = dc[:, 2].copy(), dc[:, 3].copy()
    dc[:, 2] = r2_l - e
    dc[:, 3] = r2_l - s
    dc[:, 2:4] += r1_l
    dc = np.flipud(dc)

    # Add combos
    ori.append(np.concatenate([ab, dc]))
    ori.append(np.concatenate([ba, dc]))
    ori.append(np.concatenate([ba, cd]))

    These 3 plus the original, produce two sets of identical paths; therefore only the order of reads is important
    """
    yield t

    read1_arr, read2_arr = [t[t[:, 7] == 1], t[t[:, 7] == 2]]
    read2_arr[:, 2:4] -= r1_len
    read1_arr[:, 2:4] += r2_len

    yield np.concatenate([read2_arr, read1_arr])


def process(rt):
    """
    Assumes that the reads are ordered read1 then read2, in the FR direction
    :param rt: Read_template object, contains all parameters within the pairing_params array
    """
    r1_len = rt['read1_length']
    r2_len = rt['read2_length']
    if not rt["paired_end"]:
        single_end = True
        contig_l = r1_len
    else:
        if r2_len is None and r1_len:
            single_end = True
            contig_l = r1_len
        elif r1_len is None and r2_len:
            single_end = True
            contig_l = r2_len
        elif r1_len is None and r2_len is None:
            return False
        else:
            single_end = False
            contig_l = r1_len + r2_len

    mu, sigma = rt['isize']

    pp = map(float, rt["pairing_params"])
    max_insertion = pp[0]
    min_aln = pp[1]
    max_homology = pp[2]
    ins_cost = pp[3]
    hom_cost = pp[4]  # 2
    inter_cost = pp[5]  # 10
    U = pp[6]
    match_score = rt["match_score"]

    args = [contig_l, mu, sigma, max_insertion, min_aln, max_homology, ins_cost,
            hom_cost, inter_cost, U, match_score]

    table = rt['data'][:, range(8)]

    if not single_end:

        # If it is unclear which read comes first this function can be used to generate both orientations:
        both_ways = []
        for r in read_orientations(table, r1_len, r2_len):
            a_res = align_path_c.optimal_path(r, *args)
            if len(a_res) > 0:
                if a_res[1] == a_res[4] and len(a_res[0]) == 2 and a_res[1] - a_res[2] > U:
                    # Cant do better. Normal pairing with 2 good alignments
                    return a_res
            both_ways.append(a_res)

        if len(both_ways) == 0:
            return False
        path, length, second_best, dis_to_normal, norm_pairings = sorted(both_ways, key=lambda x: x[1])[-1]  # Best

    else:
        path, length, second_best, dis_to_normal, norm_pairings = align_path_c.optimal_path(table, *args)

    if int(length) < int(second_best):
        sys.stderr.write("WARNING: primary path < secondary path\n")

    # if rt["name"] == "chr22-163734":
    #     click.echo(path, err=True)
    #     for idx, row in enumerate(rt["data"].astype(int)):
    #         click.echo((idx, list(row)), err=True)
    #     for idx, row in enumerate(rt["inputdata"]):
    #         click.echo(row, err=True)
    #     click.echo(rt["data"].astype(int).tolist(), err=True)

    return path, length, second_best, dis_to_normal, norm_pairings


if __name__ == "__main__":

    rt = {"pairing_params": (100, 17, 100, 1, 3, 2, 9), "match_score": 1, "read1_length": 100, "read2_length": 100,
          "isize": (500, 50), "name": None, "paired_end": True}

    array = np.array
    # [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, num_mis-matches]
    rt["data"] = [[0, 20462937, 0, 66, 61, 0, -1, 1, 0], [0, 20459467, 0, 35, 35, 1, 1, 1, 0], [0, 20462492, 100, 200, 95, 2, 1, 2, 0]]
    rt["data"] = [[0, 50629342, 0, 100, 96, 0, 1, 1, 0], [0, 50639946, 0, 100, 96, 1, 1, 1, 0], [0, 50630825, 144, 200, 56, 2, -1, 2, 0], [0, 50629781, 160, 200, 40, 3, -1, 2, 0]]

    # rt["data"] = [i for i in rt["data"] if i[5] in[1, 93, 102, 39, 7, 37]]
    rt["data"] = np.array(rt["data"]).astype(float)

    print(rt["data"].astype(int))
    print(process(rt))


