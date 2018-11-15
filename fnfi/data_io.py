"""
Iterate over the output from last and send data object back to the pair_reads_from_LAST script
"""

import re
from collections import defaultdict
import os
import click
import ncls
import numpy as np
import samclips


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory")


def make_template(rows, args, max_d, last_seen_chrom, fq):
    # Make a picklable data object for multiprocessing
    return {"isize": (args["insert_median"], args["insert_stdev"]),
            "max_d": max_d,
            "match_score": args["match_score"],
            "pairing_params": (args["max_insertion"], args["min_aln"], args["max_overlap"], args["ins_cost"],
                               args["ol_cost"], args["inter_cost"], args["u"]),
            "paired_end": args["paired"] == "True",
            "inputdata": rows,
            "bias": args["bias"],
            "read1_length": None,
            "read2_length": None,
            "score_mat": {},
            "passed": False,
            "name": rows[0][0][0],
            "last_seen_chrom": last_seen_chrom,
            "inputfq": fq,
            "read1_seq": None,  # Some sam records may have seq == '*' , need a record of full seq for adding back in
            "read2_seq": None,
            "read1_q": None,
            "read2_q": None,
            "read1_reverse": False,  # Set to true if aligner has reverse complemented the sequence
            "read2_reverse": False,
            "fq_read1_seq": None,
            "fq_read2_seq": None,
            "fq_read1_q": None,
            "fq_read2_q": None
            }


def get_start_end(cigar):
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    end = 0
    start = 0
    for i in range(0, len(c)-1, 2):
        if i == 0 and (c[i+1] == "S" or c[i+1] == "H"):
            start += int(c[i])
            end += int(c[i])
        elif c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += int(c[i])
    return start, end


def rev_comp(s):
    r = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "T", "c": "G", "g": "C", "t": "A", "n": "N"}
    return "".join(r[i] for i in s)[::-1]


def sam_to_array(template):

    data, overlaps = zip(*template["inputdata"])

    template["inputdata"] = [[i[1], i[2], i[3]] + i[4].strip().split("\t") for i in data]
    # [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, num_mis-matches]
    a = []
    chrom_ids = {}
    cc = 0
    if template["inputdata"][0][1] == "*":
        template["inputdata"][0][1] = template["last_seen_chrom"]

    for idx, l in enumerate(template['inputdata']):
        r = [0] * 9

        if l[1] != "*":
            chromname = l[1]
            template["last_seen_chrom"] = chromname
        else:
            l[1] = template["last_seen_chrom"]
            chromname = l[1]
        if chromname not in chrom_ids:
            chrom_ids[chromname] = cc
            cc += 1
        r[0] = chrom_ids[chromname]  # l.rname  # chrom name

        pos = int(l[2])  # Add hard clips
        r[1] = pos  # l.pos
        r[5] = idx
        flag = int(l[0])
        r[6] = -1 if flag & 16 else 1  # Flag

        # Sometimes both first and second are not flagged. Assume first
        if not flag & 64 and not flag & 128:
            r[7] = 1
        else:
            r[7] = 1 if flag & 64 else 2

        tags = [i.split(":") for i in l[11:]]
        seq_len = len(l[8])

        for k, t, v in tags:
            if k == "NM":
                r[8] = int(v)
            elif k == "AS":
                if overlaps[idx]:
                    r[4] = float(v) * template["bias"]
                else:
                    r[4] = int(v)

        cigar = l[4]
        if not cigar:
            query_start = 0  # Unmapped read, no cigar
            query_end = 0
        else:
            query_start, query_end = get_start_end(cigar)

        r[2] = query_start
        r[3] = query_end  # query_start + query_end
        a.append(r)

        if template["paired_end"]:
            if template["read1_seq"] is None and (flag & 64) and (len(l[8]) > 1):
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]

            if template["read2_seq"] is None and (flag & 128) and (len(l[8]) > 1):
                template["read2_seq"] = l[8]
                template["read2_q"] = l[9]

            if template['read1_length'] is None and flag & 64:
                template['read1_length'] = seq_len

            if template['read2_length'] is None and flag & 128:
                template['read2_length'] = seq_len

            if flag & 64 and flag & 16 and not (flag & 256 or flag & 2048):
                template["read1_reverse"] = True

            if flag & 128 and flag & 16 and not (flag & 256 or flag & 2048):
                template["read2_reverse"] = True

        else:
            if template["read1_seq"] is None and not (flag & 256) and (len(l[8]) > 1):
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]
                template["read1_length"] = len(l[8])

    # Save any input fastq information
    fq1, fq2 = template["inputfq"]
    if fq1 is not None:
        template["fq_read1_seq"] = fq1[1]
        template["fq_read1_q"] = fq1[2]
        template["read1_length"] = len(fq1[1])
    if fq2 is not None:
        template["fq_read2_seq"] = fq2[1]
        template["fq_read2_q"] = fq2[2]
        template["read2_length"] = len(fq2[1])

    if template["paired_end"]:  # Increment the contig position of read 2
        for i in range(len(a)):
            if a[i][7] == 2:  # Second in pair
                a[i][2] += template['read1_length']
                a[i][3] += template['read1_length']
            if a[i][3] > template["read1_length"] + template["read2_length"]:
                click.echo((a[i][3], template["read1_length"] + template["read2_length"]), err=True)
                click.echo((len(fq1[1]), len(fq2[1])), err=True)
                raise ValueError

    template['data'] = np.array(sorted(a, key=lambda x: (x[2], -x[4]))).astype(float)
    template['chrom_ids'] = chrom_ids

    del template["inputfq"]


def _choose_supplementary(template):

    actual_rows = [template['ri'][i] for i in template['rows']]
    d = template['data'][actual_rows, :]
    arr = [d[d[:, 7] == 1], d[d[:, 7] == 2]]  # read1, read2

    template['splitters'] = [len(v) > 1 for v in arr]  # Set to True for read1/2 to appear in output
    template['score_mat']["splitter"] = [len(v) - 1 for v in arr]

    ids_to_name = {v: k for k, v in template["chrom_ids"].items()}
    locs = []
    # Find index of highest alignment score
    for read_idx, i in enumerate(arr):
        if len(i) == 0:  # Second read is missing for unpaired
            continue
        argmax = np.argmax(i[:, 4])
        for j in range(len(i)):  # Go through each row of read
            item = i[j]
            # chrom, pos, strand, read
            loc = "{}-{}-{}-{}".format(
                ids_to_name[int(item[0])], int(item[1]), int(item[6]), int(item[7]))
            locs.append(loc)

            if loc not in template['score_mat']:
                template['score_mat'][loc] = []
            # Values are popped when setting supplementary; prevents bug where read contains two identical aligns
            if j == argmax:  # Primary, next best s
                template['score_mat'][loc] += [True, 0]
            else:
                template['score_mat'][loc] += [False, 0]

    template['locs'] = locs


def _score_alignments(template, ri):
    # Scans all alignments for each query, slow for long reads but ok for short read data
    def get_overlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))

    all_xs = []
    for item in template['rows']:
        actual_row = ri[item]
        qstart, qend = template['data'][actual_row, [2, 3]]
        size = qend - qstart
        xs = 0
        for i in range(len(template['data'])):
            if i == actual_row:
                continue
            istart, iend, iscore = template['data'][i, [2, 3, 4]]
            isize = (iend - istart) + 1e-6
            # Check for overlap
            ol = get_overlap((qstart, qend), (istart, iend))

            if ol and (ol / size > .85) and (ol / isize > .85):  # check for 85% reciprocal overlap with current
                if iscore > xs:
                    xs = iscore
        all_xs.append(xs)

    for k, v in zip(template['locs'], all_xs):
        template['score_mat'][k][1] = v


def apply_filter(template, rows, path_score, second_best, dis_to_normal, norm_pairings):

    # The rows correspond to the indexes of the input array, not the original ordering of the data
    template['rows'] = list(map(int, rows))

    # To get the original rows use
    template['ri'] = dict(zip(template['data'][:, 5], range(len(template['data']))))  # Map of row_index and array index

    template['score_mat']["dis_to_next_path"] = path_score - second_best
    template['score_mat']["dis_to_normal"] = dis_to_normal
    template['score_mat']["path_score"] = path_score
    template['score_mat']['normal_pairings'] = norm_pairings
    template['passed'] = True

    _choose_supplementary(template)
    _score_alignments(template, template['ri'])


def to_output(template):

    # Todo make sure all read-pairs have a mapping, otherwise write an unmapped
    #
    # paired = False if template["read2_length"] is None else True
    sam = samclips.fixsam(template)

    if len(sam) == 0:  # Todo fix unmapped reads
        # print("No alignments")
        return None

    if len(sam) == 1:  # Todo deal with these
        # print("Single alignment only")
        return None

    if any(i[0] == "*" or i[4] == "*" for i in sam):  # Unformatted cigar or unformatted cigarstring
        return None

    return "".join(template["name"] + "\t" + "\t".join(i) + "\n" for i in sam)


# def overlap_regions_quicksect(bed):
#     if not bed:
#         return None
#     regions = [i.split("\t")[:3] for i in open(bed, "r") if i[0] != "#"]
#     chrom_intervals = {}
#     for c, s, e in regions:
#         if c not in chrom_intervals:
#             chrom_intervals[c] = quicksect.IntervalTree()
#         chrom_intervals[c].add(int(s), int(e))
#     return chrom_intervals


def overlap_regions(bed):
    if not bed:
        return None
    regions = [i.split("\t")[:3] for i in open(bed, "r") if i[0] != "#"]
    chrom_interval_start = defaultdict(list)
    chrom_interval_end = defaultdict(list)
    for c, s, e in regions:
        chrom_interval_start[c].append(int(s))
        chrom_interval_end[c].append(int(e))

    regions = {k: ncls.NCLS(np.array(chrom_interval_start[k]),
                            np.array(chrom_interval_end[k]),
                            np.array(chrom_interval_start[k])) for k in chrom_interval_start}

    return regions


# def intersecter_quicksect(tree, chrom, start, end):
#     if tree is None:
#         return False
#     elif chrom in tree:
#         if len(tree[chrom].find(quicksect.Interval(start, end))) > 0:
#             return True


def intersecter(tree, chrom, start, end):
    if tree is None:
        return None  # Helpful for debugging
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return True
        else:
            return False
    else:
        return False


def sam_itr(args):

    itr = args["sam"]
    tree = overlap_regions(args["include"])

    # First get header
    header_string = ""
    last_seen_chrom = None
    first_line = None
    for t in itr:

        t = t.decode("utf-8")

        if t[0] == "@":
            header_string += t
            continue

        first_line = t.split("\t", 4)
        last_seen_chrom = first_line[2]

        yield header_string
        break

    pos = int(first_line[3])
    ol = intersecter(tree, first_line[2], pos, pos + 250)

    yield (first_line, last_seen_chrom, ol)

    for t in itr:
        t = t.decode("utf-8")
        line = t.split("\t", 4)

        if line[3] != last_seen_chrom:
            last_seen_chrom = line[2]

        pos = int(line[3])
        ol = intersecter(tree, line[2], pos, pos + 250)

        yield (line, last_seen_chrom, ol)


def fq_reader(args):
    # Iterate the fq files, send back a generator to use, generate only the read name, seq and qual lines
    def readfq(f):
        for l1 in f:
            l1 = l1.strip()[1:-2]  # Strip trailing /1 or /2 and leading @
            l2 = next(f).strip()
            next(f)
            l4 = next(f).strip()
            yield l1, l2, l4

    if args["fq1"] is None:
        yield None, None

    if args["fq1"] and args["fq2"] is None:  # Single end
        with open(args["fq1"]) as fq1:
            for item in readfq(fq1):
                yield (item, None)
    else:
        with open(args["fq1"]) as fq1, open(args["fq2"]) as fq2:
            for item1, item2 in zip(readfq(fq1), readfq(fq2)):
                assert item1[0] == item2[0]
                yield item1, item2


def fq_getter(reader, name, args, fbuffer):

    if args["fq1"] is None:
        return None, None

    if name in fbuffer:
        fqlines = fbuffer[name]
        del fbuffer[fqlines[name]]
        return fqlines

    while True:
        fqlines = next(reader)
        q = fqlines[0][0]
        if q == name:
            return fqlines
        else:
            fbuffer[q] = fqlines
    else:
        click.echo("WARNING: fastq not found", err=True)
        return None, None


def iterate_mappings(args):

    inputstream = sam_itr(args)

    total = 0
    name = ""
    rows = []
    header_string = next(inputstream)
    header_string += "@RG\tID:0\tSM:0\tPU:lane1\tPL:ILLUMINA\tLB:0\n"  # Todo deal with read group information

    yield header_string

    fq_buffer = defaultdict(list)
    fq_iter = fq_reader(args)

    max_d = args["insert_median"] + 4*args["insert_stdev"]  # Separation distance threshold to call a pair discordant
    last_seen_chrom = None
    for m, last_seen_chrom, ol in inputstream:  # Alignment

        nm = m[0]
        if name != nm:
            if len(rows) > 0:
                total += 1
                if total % 10000 == 0:
                    click.echo(str(total), err=True)
                fq = fq_getter(fq_iter, name, args, fq_buffer)
                yield (rows, args, max_d, last_seen_chrom, fq)

            rows = []
            name = nm

        rows.append((m, ol))  # String, ol states if alignment overlaps ROI

    # Deal with last record
    if len(rows) > 0:
        total += 1
        fq = fq_getter(fq_iter, name, args, fq_buffer)
        yield (rows, args, max_d, last_seen_chrom, fq)

    click.echo("Total processed " + str(total) + "\n", err=True)
