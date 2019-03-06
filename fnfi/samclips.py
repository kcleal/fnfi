"""
Utils to generate proper sam output and flag information
"""

import re
import click
import c_io_funcs
# from c_io_funcs import reverse_complement as rev_comp
# from c_io_funcs import get_align_end_offset
from c_samflags import set_bit


def echo(*arg):
    click.echo(arg, err=True)


def set_tlen(out):

    pri_1 = out[0][1]
    pri_2 = out[1][1]

    p1_pos = int(pri_1[2])
    p2_pos = int(pri_2[2])

    flg1 = pri_1[0]
    flg2 = pri_2[0]

    assert flg1 & 64
    assert flg2 & 128

    # Use the end position of the alignment if read is on the reverse strand, or start pos if on the forward
    if flg1 & 16:
        aln_end1 = c_io_funcs.get_align_end_offset(pri_1[4])
        t1 = p1_pos + aln_end1
    else:
        t1 = p1_pos

    if flg2 & 16:
        aln_end2 = c_io_funcs.get_align_end_offset(pri_2[4])
        t2 = p2_pos + aln_end2
    else:
        t2 = p2_pos

    if t1 <= t2:
        tlen1 = t2 - t1  # Positive
        tlen2 = t1 - t2  # Negative
    else:
        tlen1 = t2 - t1  # Negative
        tlen2 = t1 - t2  # Positive

    pri_1[7] = str(tlen1)
    pri_2[7] = str(tlen2)

    out2 = [(out[0][0], pri_1, out[0][2]), (out[1][0], pri_2, out[1][2])]

    # Set tlen's of supplementary
    for sup_tuple in out[2:]:
        sup_tuple = list(sup_tuple)
        sup_flg = sup_tuple[1][0]
        sup_chrom = sup_tuple[1][1]
        sup_pos = int(sup_tuple[1][2])

        sup_end = c_io_funcs.get_align_end_offset(sup_tuple[1][4])
        if sup_flg & 16:  # If reverse strand, count to end
            sup_pos += sup_end

        if sup_flg & 64:  # First in pair, mate is second
            other_end = t2
            other_chrom = pri_2[1]
            other_flag = pri_2[0]
        else:
            other_end = t1
            other_chrom = pri_1[1]
            other_flag = pri_1[0]
        # This is a bit of a hack to make the TLEN identical to bwa
        # Make sure they are both on same chromosome
        if sup_chrom == other_chrom:
            if sup_pos < other_end:
                if bool(sup_flg & 16) != bool(other_flag & 16):  # Different strands
                    tlen = other_end - sup_pos
                else:
                    tlen = sup_pos - other_end
            else:
                if bool(sup_flg & 16) != bool(other_flag & 16):  # Different strands
                    tlen = other_end - sup_pos
                else:
                    tlen = sup_pos - other_end
                # tlen = other_end - sup_pos
            sup_tuple[1][7] = str(tlen)
        out2.append(tuple(sup_tuple))

    return out2


def set_mate_flag(a, b, max_d, read1_rev, read2_rev):

    if not a or not b:  # No alignment, mate unmapped?
        return False, False

    # Make sure chromosome of mate is properly set not "*"
    chrom_a, mate_a = a[2], a[5]
    chrom_b, mate_b = b[2], b[5]
    if chrom_a != mate_b:
        b[5] = chrom_a
    if chrom_b != mate_a:
        a[5] = chrom_b

    aflag = a[0]
    bflag = b[0]

    reverse_A = False
    reverse_B = False

    # If set as not primary, and has been aligned to reverse strand, and primary is mapped on forward
    # the sequence needs to be rev complement
    if aflag & 256:
        if (aflag & 16) and (not read1_rev):
            reverse_A = True
        elif (not aflag & 16) and read1_rev:
            reverse_A = True

    if bflag & 256:
        if (bflag & 16) and (not read2_rev):
            reverse_B = True
        elif (not bflag & 16) and read2_rev:
            reverse_B = True
    # print bflag, read2_rev, reverse_B
    # Turn off proper pair flag, might be erroneously set
    aflag = set_bit(aflag, 1, 0)
    bflag = set_bit(bflag, 1, 0)

    # Set paired
    aflag = set_bit(aflag, 0, 1)
    bflag = set_bit(bflag, 0, 1)

    # Set first and second in pair, in case not set
    aflag = set_bit(aflag, 6, 1)
    bflag = set_bit(bflag, 7, 1)

    # Turn off any mate reverse flags, these should be reset
    aflag = set_bit(aflag, 5, 0)
    bflag = set_bit(bflag, 5, 0)

    # If either read is unmapped
    if aflag & 4:
        bflag = set_bit(bflag, 3, 1)  # Position 3, change to 1
    if bflag & 4:
        aflag = set_bit(aflag, 3, 1)

    # If either read on reverse strand
    if aflag & 16:
        bflag = set_bit(bflag, 5, 1)
    if bflag & 16:
        aflag = set_bit(aflag, 5, 1)

    # Set unmapped
    arname = a[1]
    apos = a[2]
    if apos == "0":  # -1 means unmapped
        aflag = set_bit(aflag, 2, 1)
        bflag = set_bit(bflag, 8, 1)

    brname = b[1]
    bpos = b[2]
    if b[2] == "0":
        bflag = set_bit(bflag, 2, 1)
        aflag = set_bit(aflag, 8, 1)

    # Set RNEXT and PNEXT
    a[5] = brname
    a[6] = bpos

    b[5] = arname
    b[6] = apos

    if not (apos == "-1" or bpos == "-1"):

        if arname == brname:
            # Set TLEN
            p1, p2 = int(apos), int(bpos)

            # Set proper-pair flag
            if (aflag & 16 and not bflag & 16) or (not aflag & 16 and bflag & 16):  # Not on same strand

                if abs(p1 - p2) < max_d:
                    # Check for FR or RF orientation
                    if (p1 < p2 and (not aflag & 16) and (bflag & 16)) or (p2 <= p1 and (not bflag & 16) and (aflag & 16)):
                        aflag = set_bit(aflag, 1, 1)
                        bflag = set_bit(bflag, 1, 1)

                        # If proper pair, sometimes the mate-reverse-strand flag is set
                        # this subsequently means the sequence should be reverse complemented!
                        if aflag & 16 and not bflag & 32:
                            # Mate-reverse strand not set
                            bflag = set_bit(bflag, 5, 1)
                            # reverse_B = True

                        if not aflag & 16 and bflag & 32:
                            # Mate-reverse should'nt be set
                            bflag = set_bit(bflag, 5, 0)
                            reverse_A = True

                        if bflag & 16 and not aflag & 32:
                            # Mate-reverse strand not set
                            aflag = set_bit(aflag, 5, 1)
                            # reverse_A = True

                        if not bflag & 16 and aflag & 32:
                            # Mate-revsere should'nt be set
                            aflag = set_bit(aflag, 5, 0)
                            reverse_B = True

    a[0] = aflag
    b[0] = bflag
    return reverse_A, reverse_B


def set_supp_flags(sup, pri, ori_primary_reversed, primary_will_be_reversed):

    # Set paired
    supflag = sup[0]
    priflag = pri[0]

    # Set paired and supplementary flag
    if not supflag & 1:
        supflag = set_bit(supflag, 0, 1)
    if not supflag & 2048:
        supflag = set_bit(supflag, 11, 1)

    # If primary is on reverse strand, set the mate reverse strand tag
    if priflag & 16 and not supflag & 32:
        supflag = set_bit(supflag, 5, 1)
    # If primary is on forward srand, turn off mate rev strand
    if not priflag & 16 and supflag & 32:
        supflag = set_bit(supflag, 5, 0)

    # Turn off not-primary-alignment
    if supflag & 256:
        supflag = set_bit(supflag, 8, 0)

    rev_sup = False
    if ori_primary_reversed:
        # print "Hi1"
        if not supflag & 16:  # Read on forward strand
            # if primary_will_be_reversed:
            rev_sup = True
            # else:
            #     print "Hi"

    elif supflag & 16:  # Read on reverse strand
        if not ori_primary_reversed:
            rev_sup = True

        # if not primary_will_be_reversed:
        #     rev_sup = True
        # elif primary_will_be_reversed and not ori_primary_reversed:  # Original will be moved from forward to reverse
        #     rev_sup = True
        #     print "Hi"

    elif not supflag & 16:  # Read on forward strand
        if primary_will_be_reversed and not priflag & 16:  # Primary will end up on forward
            rev_sup = True  # Old primary on reverse, so needs rev comp

    # if not ori_primary_reversed and supflag & 16:
    #     rev_sup = True
    # elif ori_primary_reversed and not supflag & 16 and not primary_will_be_reversed:
    #     rev_sup = True
    #
    # elif not priflag & 16 and primary_will_be_reversed and not supflag & 16:
    #     rev_sup = True
    # else:
    #     rev_sup = False
    # print ori_primary_reversed, priflag, supflag, primary_will_be_reversed, rev_sup
    sup[0] = supflag

    sup[5] = pri[1]
    sup[6] = pri[2]

    return rev_sup


def add_sequence_back(item, reverse_me, template):
    # item is the alignment
    flag = item[0]
    c = re.split(r'(\d+)', item[4])[1:]  # Drop leading empty string
    start = 0

    cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "DH"])
    if flag & 64:  # Read1
        seq = template["read1_seq"]
    else:
        seq = template["read2_seq"]

    if len(seq) != cigar_length:
        return item  # Cigar length is not set properly by mapper

    # Occasionally the H is missing, means its impossible to add sequence back in

    total_cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1]])
    if (flag & 64 and len(template["read1_seq"]) > total_cigar_length) or \
            (flag & 128 and len(template["read2_seq"]) > total_cigar_length):
        return item

    if flag & 64 and template["read1_seq"]:
        name = "read1"
        if template["fq_read1_seq"] != 0:
            end = len(template["fq_read1_seq"])
        else:
            end = len(template["read1_seq"])

    elif flag & 128 and template["read2_seq"]:
        name = "read2"
        if template["fq_read2_seq"] != 0:
            end = len(template["fq_read2_seq"])
        else:
            end = len(template["read2_seq"])
    else:
        return item  # read sequence is None or bad flag

    # Try and replace H with S
    if c[1] == "H" or c[-1] == "H":
        # Replace hard with soft-clips
        non_hardclipped_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "D"])
        if non_hardclipped_length == end and template["replace_hard"]:
            item[4] = item[4].replace("H", "S")
            cigar_length = non_hardclipped_length
        else:
            # Remove seq
            if c[1] == "H":
                start += int(c[0])
            if c[-1] == "H":
                end -= int(c[-2])

    # Might need to collect from the reverse direction; swap end and start
    if flag & 256 or flag & 2048:
        if flag & 64 and template["read1_reverse"] != bool(flag & 16):
            # Different strand to primary, count from end
            new_end = template["read1_length"] - start
            new_start = template["read1_length"] - end
            start = new_start
            end = new_end

        elif flag & 128 and (template["read2_reverse"] != bool(flag & 16)):
            new_end = template["read2_length"] - start
            new_start = template["read2_length"] - end
            start = new_start
            end = new_end

    # Try and use the primary sequence to replace hard-clips
    if item[9] == "*" or len(item[9]) < abs(end - start) or len(item[9]) == 0:
        if template["replace_hard"] and template["fq_%s_q" % name]:
            key = "fq_"
        else:
            key = ""
        s = template["%s%s_seq" % (key, name)][start:end]
        q = template["%s%s_q" % (key, name)][start:end]

        # print "so to add", s, len(item[9]), abs(end - start), reverse_me
        if len(s) == cigar_length:
            if reverse_me:
                item[8] = c_io_funcs.reverse_complement(s, len(s))
                item[9] = q[::-1]
            else:
                item[8] = s
                item[9] = q

    # Try and use the supplied fq file to replace the sequence
    elif template["fq_%s_q" % name] != 0 and len(template["fq_%s_q" % name]) > len(item[9]):
        if item[9] in template["fq_%s_q" % name]:
            item[8] = template["fq_%s_seq" % name][start:end]
            item[9] = template["fq_%s_q" % name][start:end]

        elif item[9] in template["fq_%s_q" % name][::-1]:
            sqn = "fq_%s_seq" % name
            s = c_io_funcs.reverse_complement(template[sqn], len(template[sqn]))[start:end]
            q = template["fq_%s_q" % name][::-1][start:end]
            if len(s) == cigar_length:
                item[8] = s
                item[9] = q

        else:
            echo("---")
            echo(item[9], flag)
            echo(name)
            echo(template["read1_q"])
            echo(template["read2_q"])
            echo(item)
            raise ValueError

    if len(item[8]) != cigar_length:
        echo(len(item[8]), cigar_length, len(item[9]), start, end)
        echo(template)

    assert len(item[8]) == cigar_length
    return item


def replace_sa_tags(alns):

    if any([i[0] == "sup" for i in alns]):
        sa_tags = {}  # Read1: tag, might be multiple split alignments
        alns2 = []
        for i, j, k in alns:
            # Remove any SA tags in alignment, might be wrong
            j = [item for idx, item in enumerate(j) if idx <= 9 or (idx > 9 and item[:2] != "SA")]
            flag = j[0]
            mapq = j[3]
            nm = 0
            chrom = j[1]
            pos = j[2]
            for tg in j[10:]:
                if tg[:2] == "NM":
                    nm = tg[5:]
                    break

            strand = "-" if flag & 16 else "+"
            cigar = j[4]
            sa = "%s,%s,%s,%s,%s,%s,%s" % (chrom, pos, strand, cigar, j[0], mapq, nm)
            key = (flag & 64, 1 if flag & 2048 else 0)
            if key in sa_tags:
                sa_tags[key] += ";" + sa
            else:
                sa_tags[key] = sa
            alns2.append([i, j, k])

        # Now add back in
        out = []
        for i, j, k in alns2:
            flag = j[0]
            key = (flag & 64, 0 if flag & 2048 else 1)
            if key in sa_tags:
                j.insert(14, "SA:Z:" + sa_tags[key])
            out.append((i, j, k))
        return out
    else:
        # Might need to remove SA tags
        return [(i, [item for idx, item in enumerate(j) if idx <= 9 or (idx > 9 and item[:2] != "SA")], ii) for i, j, ii in alns]


def fixsam(template):
    # if template["name"] == "HISEQ2500-10:541:CATW5ANXX:6:1309:7860:23396":
    #     click.echo(template, err=True)

    sam = [template['inputdata'][i] for i in template['rows']]
    max_d = template['max_d']

    # Todo make sure all read-pairs have a mapping, otherwise write an unmapped

    paired = False if template["read2_length"] is 0 else True
    score_mat = template['score_mat']

    out = []
    primary1 = 0
    primary2 = 0
    rev_A = False
    rev_B = False
    for l in sam:

        l[0] = int(l[0])  # Convert flag to int
        t = score_mat
        strand = "-1" if l[0] & 16 else "1"
        rid = str(2 if l[0] & 128 else 1)
        key = "{}-{}-{}-{}".format(l[1], l[2], strand, rid)

        if len(t[key]) > 2:
            # Prevent bug where two identical alignments possible
            aln_info_0, aln_info_1 = t[key].pop(0), t[key].pop(0)  # Remove first two items from list
        else:
            aln_info_0, aln_info_1 = t[key]
        #print key, aln_info_0, aln_info_1  # Todo no need for "splitter" field
        # if rid == "1":
        #     if t["splitter"][0]:
        #         split = "1"
        #     else:
        #         split = "0"
        # elif rid == "2":
        #     if t["splitter"][1]:
        #         split = "1"
        #     else:
        #         split = "0"

        xs = int(aln_info_1)

        l += [#"SP:Z:" + split,
              "DA:i:" + str(xs),
              "DP:Z:" + str(round(t["dis_to_next_path"], 0)),
              "DN:Z:" + str(round(t["dis_to_normal"], 2)),
              "PS:Z:" + str(round(t["path_score"], 2)),
              "NP:Z:" + str(round(t["normal_pairings"], 1)),
              ]

        if aln_info_0:
            if rid == "1":
                primary1 = l
            else:
                primary2 = l
        else:
            out.append(['sup', l, False])  # Supplementary, False to decide if rev comp
    # print
    # print primary1
    # print primary2
    if primary1 is 0 or primary2 is 0 and template["paired_end"]:
        return []  # Todo deal with unmapped read or unpaired

    if paired and template["paired_end"]:

        rev_A, rev_B = set_mate_flag(primary1, primary2, max_d, template["read1_reverse"], template["read2_reverse"])
        # print template["read1_reverse"], template["read2_reverse"]
        # print template["read1_reverse"], template["read1_seq"]
        # print rev_A
        # Check if supplementary needs reverse complementing
        for i in range(len(out)):
            if out[i][1][0] & 64:  # First in pair
                revsup = set_supp_flags(out[i][1], primary1, template["read1_reverse"], rev_A)
                # print revsup
            else:
                revsup = set_supp_flags(out[i][1], primary2, template["read2_reverse"], rev_B)
                # print revsup
                # print
            if revsup:
                out[i][2] = True

    if template["paired_end"]:
        out = [('pri', primary1, rev_A), ('pri', primary2, rev_B)] + out
        out = set_tlen(out)

    else:
        out = [('pri', primary1, rev_A)] + out

    # Add read seq info back in if necessary, before reverse complementing. Check for hard clips and clip as necessary

    for a_type, aln, reverse_me in out:

        # print a_type, reverse_me

        if aln:  # None here means no alignment for primary2
            # print ">", aln
            # Do for all supplementary
            if len(aln[8]) <= 1 or "H" in aln[4] or aln[0] & 2048:  # Sequence is set as "*", needs adding back in

                aln = add_sequence_back(aln, reverse_me, template)

            elif reverse_me:
                # print aln[8]
                aln[8] = c_io_funcs.reverse_complement(str(aln[8]), len(aln[8]))
                aln[9] = aln[9][::-1]
                # print reverse_me, rev_B, template["read2_reverse"], template["read2_seq"]
                # print aln[8]
            # Turn off not primary here
            aln[0] = set_bit(aln[0], 8, 0)

    out = replace_sa_tags(out)

    # Convert flag back to string
    for j in range(len(out)):
        out[j][1][0] = str(out[j][1][0])

    # if template["name"] == "HISEQ2500-10:539:CAV68ANXX:7:2211:10642:81376":
    #     click.echo("----", err=True)
    #     for i in out:
    #         click.echo(i, err=True)
    return [i[1] for i in out if i[1] != 0]


if __name__ == "__main__":
    import numpy as np

    array = np.array

    template = {'splitters': [False, True], 'paired_end': 1, 'locs': ['chr21-46696719-1-2', 'chr9-138234514--1-2', 'chr9-138234477-1-1'], 'bias': 1.15, 'fq_read2_seq': 0, 'isize': (250.0, 50.0), 'read2_q': 'CCCBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGG', 'max_d': 450.0, 'read2_seq': 'GTCTTTTTGTTCTCTTAACAAGGTCTCTTCCAGAGTATAAACTGTCCAGGTGCCCCCAGCCATCAATCATTTCATTAGCGTACGAAAGACACATTACTTCGTACATTCCCAAGGCTTAGCGCTCT', 'chrom_ids': {'chr7': 3, 'chr6': 14, 'chr4': 13, 'chr3': 4, 'chr2': 11, 'chr1': 16, 'chr7_KI270899v1_alt': 1, 'chr9': 0, 'chr8': 17, 'chr13': 10, 'chr3_KI270784v1_alt': 5, 'chr11': 6, 'chr10': 9, 'chr17_GL383563v3_alt': 7, 'chr16': 2, 'chr20': 19, 'chr21': 15, 'chr17': 8, 'chr19': 12, 'chr22': 18}, 'read2_length': 125, 'passed': True, 'replace_hard': 0, 'read2_reverse': 0, 'inputdata': [['65', 'chr9', '138234477', '0', '121M4S', 'chr7', '19034', '0', 'GAGTTAATGTGTGGTCTCTGCTGCAGTGTCCTGAAACAGAGCGCTAAGCCTTGGGAATGTACGAAGTAATGTGTCTTTCGTACGCTAATGAAATGATTGATGGCTGGGGGCACCTGGACAGTTTA', '=BBCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGG', 'NM:i:0', 'MD:Z:121', 'AS:i:121', 'XS:i:121'], ['337', 'chr7_KI270899v1_alt', '9205', '0', '4S121M', 'chr7', '19034', '0', '*', '*', 'NM:i:0', 'MD:Z:121', 'AS:i:121'], ['321', 'chr16', '90119483', '0', '121M4S', 'chr7', '19034', '0', '*', '*', 'NM:i:0', 'MD:Z:121', 'AS:i:121'], ['337', 'chr7', '19034', '0', '4S121M', '=', '19034', '-121', '*', '*', 'NM:i:0', 'MD:Z:121', 'AS:i:121'], ['321', 'chr3', '198170532', '0', '121M4S', 'chr7', '19034', '0', '*', '*', 'NM:i:1', 'MD:Z:82T38', 'AS:i:116'], ['337', 'chr3_KI270784v1_alt', '64908', '0', '4S121M', 'chr7', '19034', '0', '*', '*', 'NM:i:1', 'MD:Z:38A82', 'AS:i:116'], ['337', 'chr11', '178593', '0', '4S121M', 'chr7', '19034', '0', '*', '*', 'NM:i:2', 'MD:Z:38A7C74', 'AS:i:111'], ['337', 'chr17_GL383563v3_alt', '55635', '0', '5S120M', 'chr7', '19034', '0', '*', '*', 'NM:i:2', 'MD:Z:37A7C74', 'AS:i:110'], ['337', 'chr17', '115635', '0', '5S120M', 'chr7', '19034', '0', '*', '*', 'NM:i:2', 'MD:Z:37A7C74', 'AS:i:110'], ['129', 'chr7', '19034', '0', '41S84M', 'chr9', '138234477', '0', 'GTCTTTTTGTTCTCTTAACAAGGTCTCTTCCAGAGTATAAACTGTCCAGGTGCCCCCAGCCATCAATCATTTCATTAGCGTACGAAAGACACATTACTTCGTACATTCCCAAGGCTTAGCGCTCT', 'CCCBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGG', 'NM:i:0', 'MD:Z:84', 'AS:i:84', 'XS:i:84', 'SA:Z:chr10,133784666,+,45M80S,0,0;'], ['385', 'chr7_KI270899v1_alt', '9205', '0', '41S84M', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:84', 'AS:i:84'], ['401', 'chr16', '90119520', '0', '84M41S', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:84', 'AS:i:84'], ['401', 'chr9', '138234514', '0', '84M41S', '=', '138234477', '-121', '*', '*', 'NM:i:0', 'MD:Z:84', 'AS:i:84'], ['401', 'chr3', '198170569', '0', '84M41S', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:45T38', 'AS:i:79'], ['385', 'chr3_KI270784v1_alt', '64908', '0', '41S84M', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:38A45', 'AS:i:79'], ['385', 'chr11', '178593', '0', '41S84M', 'chr9', '138234477', '0', '*', '*', 'NM:i:2', 'MD:Z:38A7C37', 'AS:i:74'], ['385', 'chr17_GL383563v3_alt', '55635', '0', '42S83M', 'chr9', '138234477', '0', '*', '*', 'NM:i:2', 'MD:Z:37A7C37', 'AS:i:73'], ['385', 'chr17', '115635', '0', '42S83M', 'chr9', '138234477', '0', '*', '*', 'NM:i:2', 'MD:Z:37A7C37', 'AS:i:73'], ['2177', 'chr10', '133784666', '0', '45M80S', 'chr9', '138234477', '0', 'GTCTTTTTGTTCTCTTAACAAGGTCTCTTCCAGAGTATAAACTGTCCAGGTGCCCCCAGCCATCAATCATTTCATTAGCGTACGAAAGACACATTACTTCGTACATTCCCAAGGCTTAGCGCTCT', 'CCCBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGG', 'NM:i:0', 'MD:Z:45', 'AS:i:45', 'XS:i:45', 'SA:Z:chr7,19034,+,41S84M,0,0;'], ['385', 'chr13', '114351506', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:45', 'AS:i:45'], ['385', 'chr2', '242180757', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:45', 'AS:i:45'], ['385', 'chr19', '58605115', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:45', 'AS:i:45'], ['401', 'chr2', '113605877', '0', '80S45M', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:45', 'AS:i:45'], ['385', 'chr4', '190201964', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:45', 'AS:i:45'], ['401', 'chr17', '135875', '0', '82S43M', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:43', 'AS:i:43'], ['401', 'chr11', '191410', '0', '82S43M', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:43', 'AS:i:43'], ['401', 'chr17_GL383563v3_alt', '75875', '0', '82S43M', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:43', 'AS:i:43'], ['385', 'chr6', '170607610', '0', '43M82S', 'chr9', '138234477', '0', '*', '*', 'NM:i:0', 'MD:Z:43', 'AS:i:43'], ['385', 'chr21', '46696719', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:16G28', 'AS:i:40'], ['385', 'chr1', '248942822', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:16G28', 'AS:i:40'], ['401', 'chr8', '208357', '0', '80S45M', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:37C7', 'AS:i:40'], ['401', 'chr19', '248564', '0', '80S45M', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:37C7', 'AS:i:40'], ['385', 'chr22', '50805076', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:16G28', 'AS:i:40'], ['385', 'chr20', '64284427', '0', '43M82S', 'chr9', '138234477', '0', '*', '*', 'NM:i:1', 'MD:Z:18G24', 'AS:i:38'], ['385', 'chr17', '83246894', '0', '45M80S', 'chr9', '138234477', '0', '*', '*', 'NM:i:2', 'MD:Z:2G13G28', 'AS:i:37']], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'rows': array([28, 12,  0]), 'read1_q': '=BBCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGG', 'read1_length': 125, 'data': array([[ 8.00000000e+00,  1.15635000e+05,  0.00000000e+00,
         1.20000000e+02,  1.26499997e+02,  8.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.10000000e+02],
       [ 0.00000000e+00,  1.38234477e+08,  0.00000000e+00,
         1.21000000e+02,  1.21000000e+02,  0.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.21000000e+02],
       [ 1.00000000e+00,  9.20500000e+03,  0.00000000e+00,
         1.21000000e+02,  1.21000000e+02,  1.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.21000000e+02],
       [ 2.00000000e+00,  9.01194830e+07,  0.00000000e+00,
         1.21000000e+02,  1.21000000e+02,  2.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.21000000e+02],
       [ 3.00000000e+00,  1.90340000e+04,  0.00000000e+00,
         1.21000000e+02,  1.21000000e+02,  3.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.21000000e+02],
       [ 4.00000000e+00,  1.98170532e+08,  0.00000000e+00,
         1.21000000e+02,  1.16000000e+02,  4.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.16000000e+02],
       [ 5.00000000e+00,  6.49080000e+04,  0.00000000e+00,
         1.21000000e+02,  1.16000000e+02,  5.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.16000000e+02],
       [ 6.00000000e+00,  1.78593000e+05,  0.00000000e+00,
         1.21000000e+02,  1.11000000e+02,  6.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.11000000e+02],
       [ 7.00000000e+00,  5.56350000e+04,  0.00000000e+00,
         1.20000000e+02,  1.10000000e+02,  7.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.10000000e+02],
       [ 1.50000000e+01,  4.66967190e+07,  1.25000000e+02,
         1.70000000e+02,  4.59999990e+01,  2.80000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.00000000e+01],
       [ 9.00000000e+00,  1.33784666e+08,  1.25000000e+02,
         1.70000000e+02,  4.50000000e+01,  1.80000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.50000000e+01],
       [ 1.00000000e+01,  1.14351506e+08,  1.25000000e+02,
         1.70000000e+02,  4.50000000e+01,  1.90000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.50000000e+01],
       [ 1.10000000e+01,  2.42180757e+08,  1.25000000e+02,
         1.70000000e+02,  4.50000000e+01,  2.00000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.50000000e+01],
       [ 1.20000000e+01,  5.86051150e+07,  1.25000000e+02,
         1.70000000e+02,  4.50000000e+01,  2.10000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.50000000e+01],
       [ 1.10000000e+01,  1.13605877e+08,  1.25000000e+02,
         1.70000000e+02,  4.50000000e+01,  2.20000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.50000000e+01],
       [ 1.30000000e+01,  1.90201964e+08,  1.25000000e+02,
         1.70000000e+02,  4.50000000e+01,  2.30000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.50000000e+01],
       [ 8.00000000e+00,  1.35875000e+05,  1.25000000e+02,
         1.68000000e+02,  4.30000000e+01,  2.40000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.30000000e+01],
       [ 6.00000000e+00,  1.91410000e+05,  1.25000000e+02,
         1.68000000e+02,  4.30000000e+01,  2.50000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.30000000e+01],
       [ 7.00000000e+00,  7.58750000e+04,  1.25000000e+02,
         1.68000000e+02,  4.30000000e+01,  2.60000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.30000000e+01],
       [ 1.40000000e+01,  1.70607610e+08,  1.25000000e+02,
         1.68000000e+02,  4.30000000e+01,  2.70000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.30000000e+01],
       [ 1.60000000e+01,  2.48942822e+08,  1.25000000e+02,
         1.70000000e+02,  4.00000000e+01,  2.90000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.00000000e+01],
       [ 1.70000000e+01,  2.08357000e+05,  1.25000000e+02,
         1.70000000e+02,  4.00000000e+01,  3.00000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.00000000e+01],
       [ 1.20000000e+01,  2.48564000e+05,  1.25000000e+02,
         1.70000000e+02,  4.00000000e+01,  3.10000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.00000000e+01],
       [ 1.80000000e+01,  5.08050760e+07,  1.25000000e+02,
         1.70000000e+02,  4.00000000e+01,  3.20000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         4.00000000e+01],
       [ 1.90000000e+01,  6.42844270e+07,  1.25000000e+02,
         1.68000000e+02,  3.80000000e+01,  3.30000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         3.80000000e+01],
       [ 8.00000000e+00,  8.32468940e+07,  1.25000000e+02,
         1.70000000e+02,  3.70000000e+01,  3.40000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         3.70000000e+01],
       [ 3.00000000e+00,  1.90340000e+04,  1.66000000e+02,
         2.50000000e+02,  8.40000000e+01,  9.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         8.40000000e+01],
       [ 1.00000000e+00,  9.20500000e+03,  1.66000000e+02,
         2.50000000e+02,  8.40000000e+01,  1.00000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         8.40000000e+01],
       [ 2.00000000e+00,  9.01195200e+07,  1.66000000e+02,
         2.50000000e+02,  8.40000000e+01,  1.10000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         8.40000000e+01],
       [ 0.00000000e+00,  1.38234514e+08,  1.66000000e+02,
         2.50000000e+02,  8.40000000e+01,  1.20000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         8.40000000e+01],
       [ 4.00000000e+00,  1.98170569e+08,  1.66000000e+02,
         2.50000000e+02,  7.90000000e+01,  1.30000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.90000000e+01],
       [ 5.00000000e+00,  6.49080000e+04,  1.66000000e+02,
         2.50000000e+02,  7.90000000e+01,  1.40000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.90000000e+01],
       [ 6.00000000e+00,  1.78593000e+05,  1.66000000e+02,
         2.50000000e+02,  7.40000000e+01,  1.50000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.40000000e+01],
       [ 8.00000000e+00,  1.15635000e+05,  1.67000000e+02,
         2.50000000e+02,  8.39499983e+01,  1.70000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.30000000e+01],
       [ 7.00000000e+00,  5.56350000e+04,  1.67000000e+02,
         2.50000000e+02,  7.30000000e+01,  1.60000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.30000000e+01]]), 'name': 'HISEQ2500-10:541:CATW5ANXX:6:1309:7860:23396', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 'GAGTTAATGTGTGGTCTCTGCTGCAGTGTCCTGAAACAGAGCGCTAAGCCTTGGGAATGTACGAAGTAATGTGTCTTTCGTACGCTAATGAAATGATTGATGGCTGGGGGCACCTGGACAGTTTA', 'last_seen_chrom': 'chr17', 'ri': {0.0: 1, 1.0: 2, 2.0: 3, 3.0: 4, 4.0: 5, 5.0: 6, 6.0: 7, 7.0: 8, 8.0: 0, 9.0: 26, 10.0: 27, 11.0: 28, 12.0: 29, 13.0: 30, 14.0: 31, 15.0: 32, 16.0: 34, 17.0: 33, 18.0: 10, 19.0: 11, 20.0: 12, 21.0: 13, 22.0: 14, 23.0: 15, 24.0: 16, 25.0: 17, 26.0: 18, 27.0: 19, 28.0: 9, 29.0: 20, 30.0: 21, 31.0: 22, 32.0: 23, 33.0: 24, 34.0: 25}, 'score_mat': {'splitter': [0, 1], 'chr21-46696719-1-2': [False, 45.0], 'dis_to_next_path': 0.0, 'dis_to_normal': 236.0, 'path_score': 236.0, 'chr9-138234514--1-2': [True, 84.0], 'chr9-138234477-1-1': [True, 126.5], 'normal_pairings': 1}, 'pairing_params': (150.0, 17.0, 150.0, 0.1, 1.0, 2.0, 9.0)}

    # for item in template["inputdata"]:
    #     f = int(item[0])
    #     # Print Primary 2
    #     if not f & 2304 and f & 128:
    #         print "primary2", item

    sam = fixsam(template)

    for l in sam:
        print l

