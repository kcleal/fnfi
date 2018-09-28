#!python
#cython: language_level=2, boundscheck=False, wraparound=False
#distutils: language=c++
#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
"""
http://stackoverflow.com/questions/7403966/most-efficient-way-to-build-a-1-d-array-list-vector-of-unknown-length-using-cython
"""

import array
from cpython cimport array
import numpy as np
cimport numpy as np
from libc.math cimport exp, log, sqrt
# Import vector templates from the STL
from libcpp.vector cimport vector

#DTYPE = np.int
#ctypedef np.int_t DTYPE_t
DTYPE = np.float
ctypedef np.float_t DTYPE_t


cdef erfcc(float x):
    """Complementary error function."""
    cdef float z, t, r, p1, p2, p3, p4, p5, p6, p7, p8, p9
    z = abs(x)
    t = (1. / (1. + 0.5*z))
    p1 = -.82215223+t*.17087277
    p2 = 1.48851587+t*p1
    p3 = -1.13520398+t*p2
    p4 = .27886807+t*p3
    p5 = -.18628806+t*p4
    p6 = .09678418+t*p5
    p7 = .37409196+t*p6
    p8 = 1.00002368+t*p7
    p9 = 1.26551223+t*p8
    r = t * exp(-z*z-p9)
    if x >= 0.:
        return r
    else:
        return 2. - r


cdef normcdf(float x, float mu, float sigma):
    cdef float t, y
    t = x-mu
    y = 0.5*erfcc(-t/(sigma*sqrt(2.0)))
    if y > 1.0:
        y = 1.0
    return y


cdef bwa_pair_score(float pos1, float pos2, float strand1, float strand2, float mu, float sigma, float match_score,
                   float u=9.):
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
    cdef float d, prob, c
    cdef int proper_pair = 0

    # Check if read pair is FR type
    if not (pos1 < pos2 and (strand1 == 1. and strand2 == -1.)) \
            and not (pos2 < pos1 and (strand2 == 1. and strand1 == -1.)):
        return u, proper_pair

    # prob is probability of observing an insert size larger than d assuming a normal distribution
    d = abs(pos1 - pos2)
    if d > 2500:
        prob = 1e-9
    else:
        prob = (1 - normcdf(d, mu, sigma)) + 1e-9  # Add 1e-9 to prevent 0 and math domain error

    # Log4 is calculated as log(x)/log(base)
    c = min([-match_score * (log(prob)/log(4)), u])
    if c < 8:
        proper_pair = 1
    # bwa is [-match_score * math.log(prob, 4), 9]
    return c, proper_pair


def optimal_path(
                 np.ndarray[DTYPE_t, ndim=2] segments,
                 float contig_length,
                 float mu,
                 float sigma,
                 float max_insertion=100,
                 float min_aln=20,
                 float max_homology=100,
                 float ins_cost=1,
                 float hom_cost=3,
                 float inter_cost=2,
                 float U=9,
                 float match_score=1):
    """
    The scoring has become quite complicated.
     = Last alignment score - jump cost

    where each jump has a cost:
    jump cost = S + microhomology_length + insertion_length        # S = 10 for intra, 20 for inter

    The last score takes into account mismatched and gaps

    :param mapped:  The input dataframe with candidate alignments
    :param contig_length:   Length of the contig
    :param max_insertion: Arbitrary high number
    :param min_aln: The minimum sequence which is not part of a microhomology overlap
    :param max_homology: Arbitrary high number
    :return: A dataframe with the inferred alignments
    """

    # Start at first node then start another loop running backwards through the preceeding nodes.
    # Choose the best score out of the in edges.
    # Use a special start and end node to score the first and last alignments

    cdef np.ndarray[np.int_t, ndim=1] pred = np.zeros(segments.shape[0], dtype=np.int)
    cdef np.ndarray[np.float_t, ndim=1] node_scores = np.zeros(segments.shape[0], dtype=np.float)
    # Next best node score, for finding the secondary path
    cdef np.ndarray[np.float_t, ndim=1] nb_node_scores = np.zeros(segments.shape[0], dtype=np.float)

    normal_jumps = set([])  # Keep track of which alignments form 'normal' pairs between read-pairs. Alas native python

    cdef int i, j, p, FR
    cdef float chr1, pos1, start1, end1, score1, row_index1, strand1, r1,\
               chr2, pos2, start2, end2, score2, row_index2, strand2, r2, \
               micro_h, ins, best_score, next_best_score, best_normal_orientation, current_score, total_cost,\
               S, sc, max_s, path_score, cst, jump_cost, normal_score

    # Deal with first score
    for i in range(segments.shape[0]):
        node_scores[i] = segments[i, 4] - (segments[i, 2] * ins_cost)
    pred[0] = -1

    nb_node_scores.fill(-1e6)  # Must set to large negative, otherwise a value of zero can imply a path to that node

    best_score = 0  # Declare here in case only 1 alignment
    next_best_score = 0
    best_normal_orientation = 0  # Keep track of the best normal pairing score, F first R second

    # start from segment two because the first has been scored
    for i in range(1, segments.shape[0]):
        chr1 = segments[i, 0]
        pos1 = segments[i, 1]
        start1 = segments[i, 2]
        end1 = segments[i, 3]
        score1 = segments[i, 4]
        row_index1 = segments[i, 5]
        strand1 = segments[i, 6]
        r1 = segments[i, 7]

        p = -1  # -1 means there is no presecessor
        best_score = score1 - (start1 * ins_cost)  # Implies all preceding alignments skipped!
        next_best_score = - (start1 * ins_cost)  # Worst case

        # Walking backwards mean the search may be terminated at some point
        for j in range(i-1, -1, -1):

            chr2 = segments[j, 0]
            pos2 = segments[j, 1]
            start2 = segments[j, 2]
            end2 = segments[j, 3]
            score2 = segments[j, 4]
            row_index2 = segments[j, 5]
            strand2 = segments[j, 6]
            r2 = segments[j, 7]

            # Allow alignments with minimum sequence and max overlap
            if start1 > end2 - max_homology and end1 > end2 + min_aln and \
                                    start1 - start2 > min_aln:

                if start1 > end2 and start1 - end2 > max_insertion:
                    continue

                # Microhomology and insertion lengths between alignments on same read only
                micro_h = 0
                ins = 0
                if r1 == r2:
                    micro_h = end2 - start1
                    if micro_h < 0:
                        ins = abs(micro_h)
                        micro_h = 0

                # Define jump cost
                FR = 0

                if chr1 == chr2:
                    if r1 == r2:  # If the jump occurs on the same read, means there is an SV
                        jump_cost = U

                    else:
                        jump_cost, FR = bwa_pair_score(pos1, pos2, strand1, strand2, mu, sigma, match_score)
                        if FR:
                            normal_jumps.add((i, j))
                else:
                    jump_cost = inter_cost + U

                # Calculate score, last_score = node_scores[j]
                total_cost = (micro_h * hom_cost) + (ins * ins_cost) + jump_cost
                S = score1 - total_cost  # Score 1 is 'ahead' in the sort order from sc2

                current_score = node_scores[j] + S

                if current_score >= best_score:
                    next_best_score = best_score
                    best_score = current_score
                    p = j

                elif current_score > next_best_score:
                    next_best_score = current_score

                if FR and abs(pos1 - pos2) < 2000:  # Intra, normal pair

                    normal_score = score1 + score2 + jump_cost
                    if normal_score > best_normal_orientation:
                        best_normal_orientation = normal_score

        node_scores[i] = best_score
        pred[i] = p
        nb_node_scores[i] = next_best_score

    # Update the score for jumping to the end of the sequence
    # Basically calculates what the best and secondary scores are for the 'virtual' end node
    cdef float node_to_end_cost
    cdef float secondary

    path_score = -(contig_length * ins_cost)  # Worst case
    secondary = float(path_score)
    end_i = -1
    for i in range(segments.shape[0]):
        cst = (ins_cost * (contig_length - segments[i, 3]))

        node_to_end_cost = node_scores[i] - cst
        if node_to_end_cost > path_score:
            path_score = node_to_end_cost
            end_i = i
        elif node_to_end_cost > secondary:
            secondary = node_to_end_cost

    # Need to check if any branches from main path have a higher secondary than the 'virtual' end node
    # This can happen if the secondary path joins the main path within the graph, rather than at the end node.

    cdef float dis_2_end = 0  # Can be negative if the last node has an insertion before the end (soft-clip)
    cdef float s, potential_secondary

    # Use STL vector instead of list for possible speedup. Doesnt compile using pyximport in pairing script
    cdef vector[int] v
    cdef int normal_pairings = 0
    cdef int last_i
    cdef int next_i
    if end_i != -1:

        v.push_back(end_i)
        while True:
            last_i = v.back()
            next_i = pred[last_i]
            if next_i == -1:
                break
            v.push_back(next_i)

            if (last_i, next_i) in normal_jumps:
                normal_pairings += 1

            s = nb_node_scores[next_i]
            dis_2_end = path_score - node_scores[next_i]  # Distance of main path to the end
            potential_secondary = dis_2_end + s

            if potential_secondary > secondary:
                secondary = potential_secondary

    cdef np.ndarray[np.int_t, ndim=1] a = np.empty(len(v), dtype=np.int)
    for i in range(v.size()):
        a[v.size() - 1 - i] = v[i]  # Virtual reversal of the indexes array

    return segments[a, 5], path_score, secondary, best_normal_orientation, normal_pairings
