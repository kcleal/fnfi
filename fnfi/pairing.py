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

    # Todo second best can be negative?
    return path, length, second_best, dis_to_normal, norm_pairings


if __name__ == "__main__":
    array = np.array
    rt = {'paired_end': 1, 'bias': 1.0, 'fq_read2_seq': 0, 'isize': (210.0, 175.0), 'read2_q': '1111111111111111111111111111111111111111111111111111111111111111111111111111', 'max_d': 1085.0, 'read2_seq': 'TATTTTGAGGTTTCTTAAGGTGTTTAGCTGCAGATTTCCAATCAATAGCTGAATATTTAGATTTGAAACCTAGTAA', 'chrom_ids': {'chr7': 7, 'chr5': 8, 'chr4': 4, 'chr3': 2, 'chr2': 0, 'chr9': 1, 'chrX': 5, 'chr12': 6, 'chr10': 3}, 'read2_length': 76, 'passed': 0, 'replace_hard': 0, 'read2_reverse': 0, 'inputdata': [['65', 'chr2', '137487647', '54', '88M60S', '=', '137488063', '417', 'CAGAAATACAAACTACCATCAGAGAATAGTACAAACACCTCTACGCAAATAAAATAGAAAATCTAGAAGAAATATTTTGAGGTTTCTTAAGGTGTTTAGCTGCAGATTTCCAATCAATAGCTGAATATTTAGATTTGAAACCTAGTAA', '1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111', 'NM:i:0', 'MD:Z:88', 'AS:i:88', 'XS:i:68', 'SA:Z:chr2,137488063,+,72S76M,54,0;'], ['2113', 'chr2', '137488063', '54', '72H76M', '=', '137488063', '0', 'TATTTTGAGGTTTCTTAAGGTGTTTAGCTGCAGATTTCCAATCAATAGCTGAATATTTAGATTTGAAACCTAGTAA', '1111111111111111111111111111111111111111111111111111111111111111111111111111', 'NM:i:0', 'MD:Z:76', 'AS:i:76', 'XS:i:0', 'SA:Z:chr2,137487647,+,88M60S,54,0;'], ['321', 'chr9', '29453905', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:28C44', 'AS:i:68'], ['321', 'chr9', '26252968', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:28C44', 'AS:i:68'], ['321', 'chr3', '65512490', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:28C44', 'AS:i:68'], ['321', 'chr10', '26006565', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:53C19', 'AS:i:68'], ['321', 'chr9', '26587569', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:53C19', 'AS:i:68'], ['321', 'chr4', '16560401', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:53C19', 'AS:i:68'], ['337', 'chrX', '121697471', '0', '75H73M', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:19G53', 'AS:i:68'], ['337', 'chr9', '28114669', '0', '75H73M', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:44G28', 'AS:i:68'], ['337', 'chr4', '186113154', '0', '75H73M', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:19G53', 'AS:i:68'], ['321', 'chr12', '109000145', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:53C19', 'AS:i:68'], ['321', 'chr9', '26446877', '0', '73M75H', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:28C44', 'AS:i:68'], ['337', 'chr7', '30417234', '0', '75H73M', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:19G53', 'AS:i:68'], ['337', 'chrX', '13135077', '0', '75H73M', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:19G53', 'AS:i:68'], ['337', 'chr5', '99831714', '0', '75H73M', 'chr2', '137488063', '0', '*', '*', 'NM:i:1', 'MD:Z:19G53', 'AS:i:68'], ['129', 'chr2', '137488063', '60', '76M', '=', '137487647', '-417', 'TATTTTGAGGTTTCTTAAGGTGTTTAGCTGCAGATTTCCAATCAATAGCTGAATATTTAGATTTGAAACCTAGTAA', '1111111111111111111111111111111111111111111111111111111111111111111111111111', 'NM:i:0', 'MD:Z:76', 'AS:i:76', 'XS:i:0']], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'read1_q': '1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111', 'read1_length': 148, 'data': array([[ 0.00000000e+00,  1.37487647e+08,  0.00000000e+00,
         8.80000000e+01,  8.80000000e+01,  0.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  2.94539050e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  2.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  2.62529680e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  3.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 2.00000000e+00,  6.55124900e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  4.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 3.00000000e+00,  2.60065650e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  5.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  2.65875690e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  6.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 4.00000000e+00,  1.65604010e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  7.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 5.00000000e+00,  1.21697471e+08,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  8.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  2.81146690e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  9.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 4.00000000e+00,  1.86113154e+08,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  1.00000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 6.00000000e+00,  1.09000145e+08,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  1.10000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  2.64468770e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  1.20000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 7.00000000e+00,  3.04172340e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  1.30000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 5.00000000e+00,  1.31350770e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  1.40000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 8.00000000e+00,  9.98317140e+07,  0.00000000e+00,
         7.30000000e+01,  6.80000000e+01,  1.50000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.37488063e+08,  7.20000000e+01,
         1.48000000e+02,  7.60000000e+01,  1.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.37488063e+08,  1.48000000e+02,
         2.24000000e+02,  7.60000000e+01,  1.60000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00]]), 'name': 'HWI-D00360:5:H814YADXX:1:1202:6573:92356',
          'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 'CAGAAATACAAACTACCATCAGAGAATAGTACAAACACCTCTACGCAAATAAAATAGAAAATCTAGAAGAAATATTTTGAGGTTTCTTAAGGTGTTTAGCTGCAGATTTCCAATCAATAGCTGAATATTTAGATTTGAAACCTAGTAA', 'last_seen_chrom': 'chr2', 'score_mat': {},
          'pairing_params': (100.0, 17.0, 100.0, .5, 2.2, 14.0, 9.0)}

    print(process(rt))


