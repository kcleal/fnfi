"""
Find the best path through a collection of alignments
Works with paried-end reads rather or single contigs. Borrows the pairing heuristic from bwa.
"""
import array
import math
import time
import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.stats import norm

from align import align_path_c


def bwa_pair_score(pos1, pos2, strand1, strand2, mu, sigma, match_score, u=9):
    """
    Calculate the bwa pairing cost
    :param pos1: Position 1
    :param pos2: Position 2
    :param strand1: +1 or -1
    :param strand2: as above
    :param mu: insert size mean
    :param sigma: insert size std
    :param match_score: score gained from a matched base
    :param u: constant parameter, worst cost possible
    :return: The bwa pairing cost (float), wher the read is FR (bool)
    """

    # Check if read pair is FR type
    if not (pos1 < pos2 and (strand1 == 1 and strand2 == -1)) \
            and not (pos2 < pos1 and (strand2 == 1 and strand1 == -1)):
        return u, False

    # prob is probability of observing an insert size larger than d assuming a normal distribution
    d = abs(pos1 - pos2)
    if d > 2500:
        prob = 1e-9
    else:
        prob = (1 - norm.cdf(d, mu, sigma)) + 1e-9  # Add 1e-9 to prevent 0 and math domain error
    c = min([-match_score * math.log(prob, 4), u])
    proper_pair = True if c < 8 else False
    # bwa is [-match_score * math.log(prob, 4), 9]
    return c, proper_pair


def remove_low_qual(a):
    new = []
    for arr in a:
        m = arr[:, 8].max() - 10  # Max biased scrores
        b = arr[arr[:, 8] > m]
        new.append(b)
    print(new[0].shape[0], new[1].shape[0])
    return new


def distance_next_best_aln(d):
    # Sort the array by non-biased scores
    if len(d) > 1:
        srt_ls = sorted(zip(d[:, 4], range(len(d))), reverse=True)
        ls_diffs = [(srt_ls[i][1], srt_ls[0][0] - srt_ls[1][0]) if i == 0 else
                    (srt_ls[i][1], srt_ls[i][0] - srt_ls[0][0]) for i in range(len(srt_ls))]
    else:
        ls_diffs = [[0, d[0, 4]]]

    for i, dis in ls_diffs:
        d[i, 10] = dis
    return d


def optimal_path(segments, contig_l, mu, sigma, max_insertion=100, min_aln=20, max_homology=100, ins_cost=0,
                 hom_cost=3, inter_cost=1, U=9, match_score=1):
    """
    The scoring has become quite complicated.
     = Last alignment score - jump cost

    where each jump has a cost:
    jump cost = S + microhomology_length + insertion_length        # S = 10 for intra, 20 for inter

    The last score takes into account mismatched and gaps

    :param segments:  The input data array with candidate alignments as rows
    :param mu: Mean insert size
    :param sigma: Insert size standard deviation
    :param contig_l:   Length of the contig
    :param max_insertion: Arbitrary high number, stop looking foralignments of a big insertion is found.
    :param min_aln: The minimum sequence length allowed is not part of a microhomology overlap
    :param max_homology: Arbitrary high number
    :return: array([chosen row indexs]), path score
    """
    # Start at first node then start another loop running backwards through the preceeding nodes.
    # Choose the best score out of the in edges.
    # Use a special start and end node to score the first and last alignments

    pred = array.array('i', [0] * len(segments))  # Keep track of the predecessor for the node
    node_scores = array.array('f', [0] * len(segments))  # Keeps track of the cumulative score at each node

    # adj_dd = defaultdict(dict)  # {i: {j: {'weight'}}, ...}  A representation of the final graph

    # Deal with first score
    node_scores[0] = segments[0][4] - (segments[0][2] * ins_cost)
    pred[0] = -1

    # Segment data is [chrom, pos, query_start, query_end, score, row_index, strand, read]
    # start from segment two because the first has been scored
    for i in range(1, len(segments)):

        chr1, pos1, start1, end1, sc1, row_index1, strand1, r1 = segments[i]  # sc1 is the score
        p = -1  # -1 means there is no presecessor
        best_score = sc1 - (start1 * ins_cost)  # All preceeding alignments skipped!

        # Walking backwards means the search may be terminated at some point
        for j in range(i-1, -1, -1):

            chr2, pos2, start2, end2, sc2, row_index2, strand2, r2 = segments[j]

            if r1 == r2:
                if start1 != start2:
                    pass #print start1, start2, end1, end2

            # Allow alignments with minimum sequence and max overlap
            if start1 > end2 - max_homology and end1 > end2 + min_aln and start1 - start2 > min_aln:

                if start1 > end2 and start1 - end2 > max_insertion:
                    continue

                # Define insertion and microhomology cost
                micro_h = end2 - start1
                if micro_h < 0:
                    ins = abs(micro_h)
                    micro_h = 0
                else:
                    ins = 0

                # Define jump cost
                if chr1 == chr2:
                    if r1 == r2:  # If the jump occurs on the same read, means there is an SV
                        jump_cost = U

                    else:
                        jump_cost, FR = bwa_pair_score(pos1, pos2, strand1, strand2, mu, sigma, match_score)

                else:
                    jump_cost = inter_cost + U

                # Calculate node score, last_score = node_scores[j]
                S = - (micro_h * hom_cost) - (ins * ins_cost) - jump_cost + sc1
                current_score = node_scores[j] + S

                # adj_dd[j][i] = {'weight': 1, 'cost': S}  # Construct the actual graph using this, for debugging

                if current_score > best_score:
                    best_score = current_score
                    p = j

        node_scores[i] = best_score
        pred[i] = p

    #
    # Update the score for jumping to the end of the sequence
    max_s = -1e6
    end_i = 0
    for i in range(len(segments)):
        node_scores[i] -= (ins_cost * (contig_l - segments[i][2]))
        if node_scores[i] > max_s:
            max_s = node_scores[i]
            end_i = i
    path_score = node_scores[end_i]

    # Get path backwards from max score
    indexes = [end_i]  # [node_scores.index(max(node_scores))]
    while True:
        # Use len(indexes) - 1 to get around cython wraparound constraint (cant used indexes[-1] to get last item)
        next_i = pred[indexes[len(indexes) - 1]]
        if next_i == -1:
            break
        indexes.append(next_i)

    # convert_to_graph(adj_dd)

    # Convert indexes, to row indexes
    found = segments[indexes[::-1], 5].astype(int)

    return found, path_score


def remove_extraneous(arr, diff=10):
    """
    Not currently used, as it is quicker not to do so.
    Remove extraneous mappings, i.e. short mappings which fall within a larger mapping, or alignments with a large
    overlap but one has a low score.
    :param arr: alignments in psl style format for a single read, sorted by query_start, query_end
    :param diff: if an overlapping alignment has a score difference < x, throw away
    :return: a subset which has the shorter mappings removed
    """
    import intervaltree
    tree = intervaltree.IntervalTree()
    for idx, i in enumerate(arr):
        tree.addi(i[2], i[3], (i[4], idx))  # Start, end, score

    bad_rows = set([])  # Indexes
    for i in arr:
        overlaps = tree.search(i[2], i[3], strict=True)
        # Find which intervals to remove; interior intervals with low score
        score = i[4]
        bad_rows |= set([j[2][1] for j in overlaps if j[2][0] < (score - diff)])

    return np.delete(arr, list(bad_rows), axis=0)


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
    ori = [t]
    read1_arr, read2_arr = [t[t[:, 7] == 1], t[t[:, 7] == 2]]
    read2_arr[:, 2:4] -= r1_len
    read1_arr[:, 2:4] += r2_len

    r = np.concatenate([read2_arr, read1_arr])

    ori.append(r)

    return ori


def process(rt):
    """
    Assumes that the reads are ordered read1 then read2, in the FR direction
    :param rt: Read_template object
    :param mu: mean insert size
    :param sigma: standard deviation of insert size
    :return: path=The original row_idxs (not array indexes)
    """
    r1_len = rt['read1_length']
    r2_len = rt['read2_length']

    if r2_len is None:
        single_end = True
        contig_l = r1_len
    else:
        single_end = False
        contig_l = r1_len + r2_len

    mu, sigma = rt['isize']

    max_insertion = 100
    min_aln = 17
    max_homology = 100
    ins_cost = 1
    hom_cost = 2
    inter_cost = 10
    U = 9
    match_score = 1

    args = [contig_l, mu, sigma, max_insertion, min_aln, max_homology, ins_cost,
            hom_cost, inter_cost, U, match_score]

    table = rt['data'][:, range(8)]

    if not single_end:
        # If it is unclear which read comes first this function can be used to generate both orientations:
        orientations = read_orientations(table, r1_len, r2_len)
        both_ways = [align_path_c.optimal_path(i, *args) for i in orientations]

        # Remove any bad paths
        both_ways = [i for i in both_ways if len(i) > 0]
        if len(both_ways) == 0:
            return False
        path, length, second_best, dis_to_normal, norm_pairings = sorted(both_ways, key=lambda x: x[1])[-1]  # Best

    else:
        path, length, second_best, dis_to_normal, norm_pairings = align_path_c.optimal_path(table, *args)

    if int(length) < int(second_best):
        sys.stderr.write("WARNING: primary path < secondary path\n")

    return path, length, second_best, dis_to_normal, norm_pairings


if __name__ == "__main__":

    # Components are
    # chromID, pos, qstart, qend, alignScore, readID
    # Read1 and Read2 are concatenated in both orientations as ordering is unknown.

    # Use mean template and std to find a p-value for observed template lengths for bwa pairing heuristic
    mean_template = 400  # These should be empirically found
    template_std = 150

    # Data is [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, num_matches, biased]
    # Data should be sorted by the query_start column
    # Path is 2,4,6
    # Todo make a test out of this
    read1 = [[0, 1100, 0, 75, 75, 0, -1,1,40],
             #[0, 1200, 0, 16, 16, 1, -1,1, 16],
             [0, 1100, 0, 100, 100, 1, -1,1, 80],
             #[4, 256700, 40, 65, 20, 3, 1,1, 20]
             ]

    read2 = [[0, 1000, 100, 200, 100, 2, 1,2, 100],
             #[13, 5500, 100, 200, 100, 5, -1,2, 100],
             #[0, 800, 200, 250, 50, 6, -1,2, 50]
             ]


    class Template():
        def __init__(self):
            self.data = np.array(read1 + read2).astype(float)
            self.read1_length = 100
            self.read2_length = 100

    t0 = time.time()
    print(process(Template(), mu=mean_template, sigma=template_std))
    print(time.time() - t0)