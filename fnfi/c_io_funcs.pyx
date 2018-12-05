#!python
#cython: language_level=2, boundscheck=False, wraparound=False
#distutils: language=c++


import numpy as np
cimport numpy as np
import click

DTYPE = np.float
ctypedef np.float_t DTYPE_t


from libc.stdlib cimport malloc
import re

cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0',
                       '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0',  'a',  'a' ]


def reverse_complement(str seq, int seq_len):
    """https://bioinformatics.stackexchange.com/questions/3583/\
    what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho/3595#3595"""

    cdef char *seq_dest = <char *>malloc(seq_len + 1)
    seq_dest[seq_len] = '\0'

    cdef bytes py_bytes = seq.encode('UTF-8')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('UTF-8')


def get_start_end(str cigar):
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    cdef int end = 0
    cdef int start = 0
    cdef int i

    for i in range(0, len(c)-1, 2):
        if i == 0 and (c[i+1] == "S" or c[i+1] == "H"):
            start += int(c[i])
            end += int(c[i])
        elif c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += int(c[i])
    return start, end


def sam_to_array(template):

    data, overlaps = zip(*template["inputdata"])
    template["inputdata"] = [[i[1], i[2], i[3]] + i[4].strip().split("\t") for i in data]

    # [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, num_mis-matches]
    cdef np.ndarray[np.float_t, ndim=2] arr = np.zeros((len(data), 9), dtype=np.float)

    chrom_ids = {}

    if template["inputdata"][0][1] == "*":
        template["inputdata"][0][1] = template["last_seen_chrom"]

    cdef int cc = 0
    cdef int idx, pos, flag, seq_len, query_start, query_end, new_start, new_end
    cdef str chromname, cigar, k, t, v
    cdef float bias = template["bias"]

    cdef int read1_strand_set, read2_strand_set, current_l
    read1_set = 0  # Occationally multiple primaries, take the longest
    read2_set = 0

    for idx in range(len(template["inputdata"])):

        l = template["inputdata"][idx]

        if l[1] != "*":
            chromname = l[1]
            template["last_seen_chrom"] = chromname
        else:
            l[1] = template["last_seen_chrom"]
            chromname = l[1]
        if chromname not in chrom_ids:
            chrom_ids[chromname] = cc
            cc += 1
        arr[idx, 0] = chrom_ids[chromname]  # l.rname  # chrom name

        pos = int(l[2])  # Add hard clips
        arr[idx, 1] = pos  # l.pos
        arr[idx, 5] = idx
        flag = int(l[0])
        arr[idx, 6] = -1 if flag & 16 else 1  # Flag

        # Sometimes both first and second are not flagged. Assume first
        if not flag & 64 and not flag & 128:
            arr[idx, 7] = 1
        else:
            arr[idx, 7] = 1 if flag & 64 else 2

        tags = [i.split(":") for i in l[11:]]
        seq_len = len(l[8])

        for k, t, v in tags:
            if k == "NM":
                arr[idx, 8] = float(v)
            elif k == "AS":
                if overlaps[idx]:
                    arr[idx, 4] = float(v) * bias
                else:
                    arr[idx, 4] = float(v)

        cigar = l[4]
        if not cigar:
            query_start = 0  # Unmapped read, no cigar
            query_end = 0
        else:
            query_start, query_end = get_start_end(cigar)

            # If current alignment it not primary, and on different strand from primary, count from other end
            if flag & 256 or flag & 2048:
                if flag & 64 and template["read1_reverse"] != bool(flag & 16):
                        # Different strand to primary, count from end
                        new_end = template["read1_length"] - query_start
                        new_start = template["read1_length"] - query_end
                        query_start = new_start
                        query_end = new_end

                elif flag & 128 and (template["read2_reverse"] != bool(flag & 16)):
                    new_end = template["read2_length"] - query_start
                    new_start = template["read2_length"] - query_end
                    query_start = new_start
                    query_end = new_end

        arr[idx, 2] = query_start
        arr[idx, 3] = query_end  # query_start + query_end

        current_l = len(l[9])

        if template["paired_end"]:
            if flag & 64 and read1_set < current_l and len(l[8]) > 1:  # First in pair
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]
                template["read1_length"] = seq_len
                read1_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template["read1_reverse"] = 1 if flag & 16 else 0

            elif flag & 128 and read2_set < current_l and len(l[8]) > 1:  # First in pair
                template["read2_seq"] = l[8]
                template["read2_q"] = l[9]
                template["read2_length"] = seq_len
                read2_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template["read2_reverse"] = 1 if flag & 16 else 0

        else:
            if template["read1_seq"] == 0 and not (flag & 256) and (len(l[8]) > 1) and read1_set < current_l:
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]
                template["read1_length"] = len(l[8])
                read1_set = current_l

    # Save any input fastq information
    fq1, fq2 = template["inputfq"]
    if fq1 == 1:
        template["fq_read1_seq"] = fq1[1]
        template["fq_read1_q"] = fq1[2]
        template["read1_length"] = len(fq1[1])
    if fq2 == 1:
        template["fq_read2_seq"] = fq2[1]
        template["fq_read2_q"] = fq2[2]
        template["read2_length"] = len(fq2[1])

    cdef int j
    if template["paired_end"]:  # Increment the contig position of read 2

        for j in range(len(arr)):
            if arr[j, 7] == 2:  # Second in pair
                arr[j, 2] += template['read1_length']
                arr[j, 3] += template['read1_length']
            if arr[j, 3] > template["read1_length"] + template["read2_length"]:
                click.echo((template["read1_length"], template["read2_length"]), err=True)
                click.echo(arr[j, 3], err=True)
                raise ValueError

    template['data'] = np.array(sorted(arr, key=lambda x: (x[2], -x[4]))).astype(float)
    template['chrom_ids'] = chrom_ids

    # if template["name"] == u'HISEQ1:9:H8962ADXX:1:1101:1144:90311':
    #     click.echo(template["read1_seq"], err=True)
    #     quit()
    del template["inputfq"]


def choose_supplementary(template):
    template['ri'] = dict(zip(template['data'][:, 5], range(len(template['data']))))  # Map of row_index and array index
    actual_rows = [template['ri'][j] for j in template['rows']]

    cdef np.ndarray[np.float_t, ndim=2] d = template['data'][actual_rows, :]

    cdef float read1_max = 0.
    cdef float read2_max = 0.
    cdef int i = 0
    cdef int read1_alns = 0
    cdef int read2_alns = 0

    for i in range(len(d)):
        if d[i, 7] == 1 and d[i, 4] > read1_max:
            read1_max = d[i, 4]
            read1_alns += 1
        elif d[i, 7] == 2 and d[i, 4] > read2_max:
            read2_max = d[i, 4]
            read2_alns += 1

    template["splitters"] = [read1_alns > 1, read2_alns > 1]
    template['score_mat']["splitter"] = [read1_alns - 1, read2_alns - 1]

    ids_to_name = {v: k for k, v in template["chrom_ids"].items()}

    locs = []
    cdef float m
    for i in range(len(d)):

        loc = "{}-{}-{}-{}".format(ids_to_name[int(d[i, 0])], int(d[i, 1]), int(d[i, 6]), int(d[i, 7]) )
        locs.append(loc)

        if loc not in template['score_mat']:
                template['score_mat'][loc] = []
        # Values are popped when setting supplementary; prevents bug where read contains two identical aligns
        if d[i, 7] == 1:
            m = read1_max
        else:
            m = read2_max

        if np.round(d[i, 4], 3) == m:  # Primary, next best s
            template['score_mat'][loc] += [True, 0]
        else:
            template['score_mat'][loc] += [False, 0]
    template['locs'] = locs


def score_alignments(template, ri, np.ndarray[np.int64_t, ndim=1]  template_rows, np.ndarray[DTYPE_t, ndim=2] template_data):
    # Scans all alignments for each query, slow for long reads but ok for short read data

    all_xs = []
    cdef int i, actual_row, item, idx
    cdef float xs = 0
    cdef float size = 0
    cdef float qstart = 0
    cdef float qend = 0
    cdef float isize = 0
    cdef float istart = 0
    cdef float iend = 0
    cdef float iscore = 0
    cdef float ol = 0

    idx = 0
    for item in template_rows:
        actual_row = ri[item]
        qstart, qend = template_data[actual_row, [2, 3]]
        size = qend - qstart
        xs = 0
        for i in range(len(template_data)):
            if i == actual_row:
                continue
            istart, iend, iscore = template_data[i, [2, 3, 4]]
            isize = (iend - istart) + 1e-6
            # Check for overlap
            ol = max(0, min(qend, iend) - max(qstart, istart))
            if ol and (ol / size > .85) and (ol / isize > .85):  # check for 85% reciprocal overlap with current
                if iscore > xs:
                    xs = iscore

        template["score_mat"][template["locs"][idx]][1] = xs
        idx += 1


def add_scores(template, np.ndarray[np.float_t, ndim=1] rows, float path_score, float second_best, float dis_to_normal, int norm_pairings):

    # The rows correspond to the indexes of the input array, not the original ordering of the data
    template['rows'] = rows.astype(int)  # list(map(int, rows))

    # To get the original rows use, col 5 is the row index: mapped to actual index
    template['score_mat']["dis_to_next_path"] = path_score - second_best
    template['score_mat']["dis_to_normal"] = dis_to_normal
    template['score_mat']["path_score"] = path_score
    template['score_mat']['normal_pairings'] = norm_pairings


