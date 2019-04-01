#!python
#cython: language_level=2, boundscheck=False, wraparound=False

cimport cpython.array

cdef mult(float a, float b):
    return a * b


cpdef float alignments_match(str rseq, str tseq, int r_left_clip, int t_left_clip, int r_pos, int t_pos,
                     char [:] r_quals, char [:] t_quals, int max_mismatches=2):

    cdef int len_r_seq, len_t_seq, r_start, t_start, end, mis, ri, rj, i
    cdef float p

    len_r_seq = len(rseq)
    len_t_seq = len(tseq)

    if r_left_clip > 0:
        r_pos -= r_left_clip

    if t_left_clip > 0:
        t_pos -= t_left_clip

    r_start = 0
    t_start = 0

    if r_pos < t_pos:
        r_start = t_pos - r_pos
        end = min(len_t_seq, len_r_seq - r_start)
    elif t_pos < r_pos:
        t_start = r_pos - t_pos
        end = min(len_r_seq, len_t_seq - t_start)
    else:
        end = min(len_r_seq, len_t_seq)

    mis = 0
    p = 1.
    for i in range(end):
        ri = r_start + i
        ti = t_start + i

        if rseq[ri] != tseq[ti]:
            mis += 1
            if mis > max_mismatches:
                return 0.

            minqual = float(min(r_quals[ri], t_quals[ti]))
            p *= 10 ** (-minqual / 10)

    return p