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

    # if rt["name"] == "simulated_reads.0.10-id217_A_chr21:46697360_B_chr6:157282148-28985":
    #     click.echo(path, err=True)
    #     for idx, row in enumerate(rt["data"].astype(int)):
    #         click.echo((idx, list(row)), err=True)
    #     click.echo(rt["data"].astype(int).tolist(), err=True)

    return path, length, second_best, dis_to_normal, norm_pairings


if __name__ == "__main__":

    rt = {"pairing_params": (100, 17, 100, 1, 3, 2, 9), "match_score": 1, "read1_length": 125, "read2_length": 125,
          "isize": (300, 150), "name": None}

    array = np.array

    rt["data"] = [[1, 46697206, 0, 125, 141, 1, 1, 1, 0], [0, 248943309, 0, 125, 123, 0, 1, 1, 0], [2, 50805563, 0, 125, 123, 2, 1, 1, 0], [3, 113605310, 0, 125, 113, 3, -1, 1, 0], [4, 59802900, 0, 87, 82, 4, -1, 1, 0], [2, 40775679, 0, 87, 82, 5, -1, 1, 0], [5, 43980445, 0, 87, 82, 6, -1, 1, 0], [0, 204171422, 0, 81, 81, 7, -1, 1, 0], [6, 98606459, 0, 81, 81, 8, -1, 1, 0], [7, 94343709, 0, 81, 81, 9, -1, 1, 0], [8, 9134289, 0, 80, 80, 10, -1, 1, 0], [9, 2167926, 0, 81, 76, 11, -1, 1, 0], [3, 187604049, 0, 84, 69, 15, -1, 1, 0], [6, 101163814, 0, 88, 68, 16, -1, 1, 0], [10, 18096396, 0, 83, 68, 17, -1, 1, 0], [5, 64586379, 0, 88, 68, 19, -1, 1, 0], [7, 73759684, 0, 89, 64, 23, -1, 1, 0], [5, 142667073, 0, 88, 63, 24, -1, 1, 0], [13, 45671280, 0, 90, 60, 26, -1, 1, 0], [15, 373402, 0, 90, 60, 27, -1, 1, 0], [0, 40454613, 0, 90, 50, 31, -1, 1, 0], [13, 1201895, 2, 86, 67, 20, -1, 1, 0], [13, 8609030, 3, 94, 53, 30, -1, 1, 0], [8, 22647563, 6, 93, 66, 22, -1, 1, 0], [11, 187313656, 9, 93, 75, 13, -1, 1, 0], [12, 122814974, 9, 96, 68, 18, -1, 1, 0], [14, 50877744, 37, 139, 41, 36, 1, 1, 0], [11, 39676689, 38, 163, 42, 35, 1, 1, 0], [10, 152311907, 39, 164, 76, 12, 1, 1, 0], [10, 125458533, 39, 164, 71, 14, 1, 1, 0], [11, 101529556, 39, 164, 66, 21, 1, 1, 0], [10, 186980476, 39, 164, 46, 34, 1, 1, 0], [14, 70721343, 40, 165, 60, 25, 1, 1, 0], [16, 50737539, 41, 154, 57, 29, 1, 1, 0], [0, 213916185, 41, 166, 49, 32, 1, 1, 0], [12, 154777980, 41, 166, 49, 33, 1, 1, 0], [7, 38079164, 42, 167, 58, 28, 1, 1, 0], [16, 248081, 84, 209, 34, 37, -1, 1, 0], [17, 207874, 84, 206, 31, 38, -1, 1, 0], [18, 157282298, 125, 250, 120, 39, -1, 2, 0], [5, 18811230, 125, 237, 82, 41, -1, 2, 0], [8, 90164704, 125, 237, 74, 53, -1, 2, 0], [10, 125390819, 125, 223, 72, 55, -1, 2, 0], [18, 161490568, 125, 237, 72, 56, -1, 2, 0], [7, 96000479, 125, 223, 70, 60, -1, 2, 0], [16, 21838577, 125, 223, 70, 62, -1, 2, 0], [16, 7128391, 125, 213, 70, 63, -1, 2, 0], [12, 150478530, 125, 227, 70, 65, -1, 2, 0], [13, 45135143, 125, 237, 69, 66, -1, 2, 0], [10, 128017277, 125, 238, 68, 68, -1, 2, 0], [6, 86259045, 125, 227, 67, 69, -1, 2, 0], [10, 138550540, 125, 237, 67, 70, -1, 2, 0], [4, 63004434, 125, 243, 65, 71, -1, 2, 0], [7, 95184760, 125, 250, 63, 73, -1, 2, 0], [16, 58101697, 125, 250, 62, 76, -1, 2, 0], [14, 83084687, 125, 249, 61, 79, -1, 2, 0], [3, 84943598, 125, 221, 61, 80, 1, 2, 0], [23, 75905083, 125, 240, 47, 87, -1, 2, 0], [0, 8142172, 125, 240, 42, 91, -1, 2, 0], [0, 211533467, 125, 228, 41, 94, 1, 2, 0], [23, 104771997, 125, 240, 32, 101, -1, 2, 0], [25, 54275176, 125, 167, 28, 103, 1, 2, 0], [0, 145454989, 125, 156, 26, 104, 1, 2, 0], [19, 67803368, 126, 251, 74, 52, 1, 2, 0], [9, 3312166, 126, 245, 68, 67, -1, 2, 0], [23, 87274832, 126, 229, 62, 78, 1, 2, 0], [4, 110688626, 126, 227, 33, 100, 1, 2, 0], [2, 30540737, 127, 252, 56, 83, 1, 2, 0], [6, 129916366, 127, 199, 42, 93, 1, 2, 0], [23, 133403648, 127, 182, 39, 96, 1, 2, 0], [7, 14599999, 127, 176, 34, 98, 1, 2, 0], [7, 22011814, 128, 177, 36, 97, 1, 2, 0], [16, 42264889, 129, 241, 78, 42, -1, 2, 0], [3, 172180973, 129, 244, 58, 81, -1, 2, 0], [0, 51028083, 130, 255, 47, 85, 1, 2, 0], [18, 107283755, 131, 256, 64, 72, 1, 2, 0], [3, 60447349, 131, 234, 62, 77, 1, 2, 0], [24, 84295617, 131, 234, 42, 92, 1, 2, 0], [19, 4912827, 132, 257, 85, 40, 1, 2, 0], [11, 155876992, 135, 238, 63, 74, 1, 2, 0], [22, 41615, 135, 238, 63, 75, 1, 2, 0], [11, 24058675, 137, 262, 51, 84, 1, 2, 0], [2, 41921621, 138, 263, 77, 43, 1, 2, 0], [12, 124858366, 138, 263, 77, 44, 1, 2, 0], [0, 212470930, 138, 241, 75, 46, 1, 2, 0], [3, 229960838, 138, 263, 75, 50, 1, 2, 0], [17, 28968903, 138, 263, 75, 51, 1, 2, 0], [3, 235293681, 138, 260, 72, 54, 1, 2, 0], [12, 135068333, 138, 263, 71, 58, 1, 2, 0], [21, 57130108, 138, 263, 70, 61, 1, 2, 0], [24, 77976321, 139, 262, 47, 86, -1, 2, 0], [14, 30664946, 147, 259, 75, 47, -1, 2, 0], [20, 2950483, 147, 259, 75, 48, -1, 2, 0], [19, 8007483, 147, 259, 75, 49, -1, 2, 0], [17, 119545508, 147, 260, 71, 59, -1, 2, 0], [3, 33474502, 147, 269, 70, 64, -1, 2, 0], [21, 22859021, 147, 271, 45, 88, -1, 2, 0], [11, 71159563, 147, 269, 45, 89, -1, 2, 0], [11, 64836968, 147, 271, 44, 90, -1, 2, 0], [8, 94128350, 147, 272, 41, 95, -1, 2, 0], [5, 48635056, 148, 273, 76, 45, 1, 2, 0], [5, 136846542, 151, 270, 57, 82, -1, 2, 0], [0, 183153560, 151, 275, 33, 99, -1, 2, 0], [3, 169663932, 152, 277, 71, 57, 1, 2, 0], [13, 28412837, 201, 326, 29, 102, -1, 2, 0]]

    rt["data"] = [i for i in rt["data"] if i[5] in[1, 93, 102, 39, 7, 37]]
    rt["data"] = np.array(rt["data"]).astype(float)

    print(rt["data"].astype(int))
    print(process(rt))


