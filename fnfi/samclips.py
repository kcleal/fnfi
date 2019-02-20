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

    if not ori_primary_reversed and supflag & 16:
        rev_sup = True
    elif ori_primary_reversed and not supflag & 16 and not primary_will_be_reversed:
        rev_sup = True
    elif not priflag & 16 and primary_will_be_reversed and not supflag & 16:
        rev_sup = True
    else:
        rev_sup = False

    sup[0] = supflag

    sup[5] = pri[1]
    sup[6] = pri[2]

    return rev_sup


def add_sequence_back(item, reverse_me, template):
    # item is the alignment
    c = re.split(r'(\d+)', item[4])[1:]  # Drop leading empty string
    start = 0

    cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "DH"])
    if len(item[8]) != cigar_length:
        return  # Cigar length is not set properly by mapper

    # Occasionally the H is missing, means its impossible to add sequence back in
    flag = item[0]
    total_cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1]])
    if (flag & 64 and len(template["read1_seq"]) > total_cigar_length) or \
            (flag & 128 and len(template["read2_seq"]) > total_cigar_length):
        return

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
        return  # read sequence is None or bad flag

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

    sam = [template['inputdata'][i] for i in template['rows']]
    r1l, r2l = template["read1_length"], template["read2_length"]
    max_d = template['max_d']

    # if template["name"] == "chr22-1878591":
    #     click.echo(template, err=True)
        # quit()

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

        if rid == "1":
            if t["splitter"][0]:
                split = "1"
            else:
                split = "0"
        elif rid == "2":
            if t["splitter"][1]:
                split = "1"
            else:
                split = "0"

        xs = int(aln_info_1)

        l += ["SP:Z:" + split,
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

    if primary1 is 0 or primary2 is 0 and template["paired_end"]:
        return []  # Todo deal with unmapped read or unpaired

    if paired and template["paired_end"]:

        rev_A, rev_B = set_mate_flag(primary1, primary2, max_d, template["read1_reverse"], template["read2_reverse"])

        # Check if supplementary needs reverse complementing
        for i in range(len(out)):
            if out[i][1][0] & 64:  # Flag
                revsup = set_supp_flags(out[i][1], primary2, template["read1_reverse"], rev_A)
            else:
                revsup = set_supp_flags(out[i][1], primary1, template["read2_reverse"], rev_B)

            if revsup:
                out[i][2] = True

    if template["paired_end"]:
        out = [('pri', primary1, rev_A), ('pri', primary2, rev_B)] + out
        out = set_tlen(out)

    else:
        out = [('pri', primary1, rev_A)] + out

    # Add read seq info back in if necessary, before reverse complementing. Check for hard clips and clip as necessary
    for a_type, aln, reverse_me in out:

        if aln:  # None here means no alignment for primary2

            if len(aln[8]) <= 1 or "H" in aln[4]:  # Sequence is set as "*", needs adding back in
                add_sequence_back(aln, reverse_me, template)
            elif reverse_me:
                aln[8] = c_io_funcs.reverse_complement(str(aln[8]), len(aln[8]))
                aln[9] = aln[9][::-1]

            # Turn off not primary here
            aln[0] = set_bit(aln[0], 8, 0)

    out = replace_sa_tags(out)

    # Convert flag back to string
    for j in range(len(out)):
        out[j][1][0] = str(out[j][1][0])

    # if template["name"] == "chr22-59941":
    #
    #     for item in out:
    #         echo(item)
    #
    #     for item in template["data"]:
    #         echo(list(item.astype(int)))
    #
    #     for item in template["inputdata"]:
    #         echo(item)

    return [i[1] for i in out if i[1] != 0]


if __name__ == "__main__":
    import numpy as np

    array = np.array

    template = {'splitters': [True, False], 'paired_end': 1, 'locs': ['chr4-49139696--1-2', 'chrUn_KI270442v1-85382--1-2', 'chrUn_KI270442v1-92912--1-1', 'chrUn_GL000216v2-23597--1-1'], 'bias': 1.0, 'fq_read2_seq': 0, 'isize': (545.0, 160.0), 'read2_q': '@@>@BBBBBBBBBBBBAAAAAAAAAAAAAAAAA@AAAA@@@A@@@@?@9@?@@A@@@@@?@@@???????@@?>>>?>8>=<===<=<==<====<=<=6<=<==<====<=<==<=<==<=<=>====><=<==????A@@@??<>?', 'max_d': 1185.0, 'read2_seq': 'ATTCCATTCCATTCCATTCCATTCTATTCCATTTCATTCCATTCCACTCGGTTTGATTCCATTCCATTCCATTCCACACGGGTTTATTCCATTTTATTCCGTTCCATCCCATTCCATTCCATTCCTTTTAGTTCCATTCCATTCCATT', 'chrom_ids': {'chr2_KI270772v1_alt': 15, 'chr14_GL000225v1_random': 21, 'chrY': 11, 'chr22_KI270736v1_random': 3, 'chr22_KI270737v1_random': 1, 'chr10': 4, 'chr17': 5, 'chr2_KI270894v1_alt': 14, 'chrUn_KI270442v1': 6, 'chr17_KI270730v1_random': 16, 'chr22': 0, 'chr20': 9, 'chr21': 17, 'chr17_KI270729v1_random': 10, 'chr5': 12, 'chr4': 8, 'chr2': 7, 'chr1': 19, 'chr9': 18, 'chrUn_KI270438v1': 22, 'chrUn_GL000216v2': 2, 'chrUn_KI270757v1': 23, 'chr10_KI270824v1_alt': 20, 'chrUn_KI270756v1': 13}, 'read2_length': 148, 'passed': True, 'replace_hard': 0, 'read2_reverse': 1, 'inputdata': [['97', 'chr22', '16352731', '0', '62M10I21M', 'chr2', '89830945', '0', 'AATGGAATGGAATGGAACTAAAAGGAATGGAATGGAATGGGATGGAACGGAATAAAATGGAATAAACCCGTGTGGAATGGAATGGAATGGAAT', '?><??@@@A????==<=<>====>=<=<==<=<==<=<====<==<=<6=<=<====<==<=<===<=>8>?>>>?@@???????@@@?@@@@', 'NM:i:16', 'MD:Z:4A11C1C21A12G0G28', 'AS:i:37', 'XS:i:37', 'RG:Z:0'], ['369', 'chr22_KI270737v1_random', '32986', '0', '20M10I63M', 'chr2', '89830945', '0', '*', '*', 'NM:i:16', 'MD:Z:28C0C12T21G1G11T4', 'AS:i:37', 'RG:Z:0'], [369, 'chrUn_GL000216v2', '23597', '0', '33M10D37M23H', 'chr4', '49139696', '49116105', u'ATGGAATGGAATGGAATGGAAAGAAATCAACTCGACTGGAATAGAATTGAAGGGAATGGAGTGAAATGGAATG', '9:::86?<<89;;;;<>=<<=>=>==>>?><>3==?=>==<>><==>><<>>=<<<===>;<;>;;=<;;;<=', 'NM:i:14', 'MD:Z:29G3^CCATTCCATG12A6T11T5', 'AS:i:34', 'RG:Z:0', 'SP:Z:1', 'DA:i:0', 'DP:Z:0.0', 'DN:Z:111.5', 'PS:Z:127.05', 'NP:Z:0.0'], ['353', 'chr22_KI270736v1_random', '147858', '0', '53M40H', 'chr2', '89830945', '0', '*', '*', 'NM:i:4', 'MD:Z:17T0G21A6T5', 'AS:i:33', 'RG:Z:0'], ['369', 'chr10', '42315909', '0', '46H47M', 'chr2', '89830945', '0', '*', '*', 'NM:i:3', 'MD:Z:6T10T9C19', 'AS:i:32', 'RG:Z:0'], ['353', 'chr17', '26785195', '0', '40M53H', 'chr2', '89830945', '0', '*', '*', 'NM:i:2', 'MD:Z:17T0G21', 'AS:i:30', 'RG:Z:0'], ['353', 'chr22', '16344196', '0', '19M5I38M10I21M', 'chr2', '89830945', '0', '*', '*', 'NM:i:20', 'MD:Z:35A7T4G0G19C8', 'AS:i:30', 'RG:Z:0'], ['369', 'chrUn_GL000216v2', '35070', '0', '53H40M', 'chr2', '89830945', '0', '*', '*', 'NM:i:2', 'MD:Z:21G0A17', 'AS:i:30', 'RG:Z:0'], ['353', 'chrUn_KI270442v1', '13799', '0', '40M53H', 'chr2', '89830945', '0', '*', '*', 'NM:i:2', 'MD:Z:6T11G21', 'AS:i:30', 'RG:Z:0'], ['369', 'chr22_KI270737v1_random', '41560', '0', '20M10I39M5I19M', 'chr2', '89830945', '0', '*', '*', 'NM:i:20', 'MD:Z:8G19C0C4A7T35', 'AS:i:30', 'RG:Z:0'], ['145', 'chr2', '89830945', '0', '76M72S', 'chr22', '16352731', '0', 'ATTCCATTCCATTCCATTCCATTCTATTCCATTTCATTCCATTCCACTCGGTTTGATTCCATTCCATTCCATTCCACACGGGTTTATTCCATTTTATTCCGTTCCATCCCATTCCATTCCATTCCTTTTAGTTCCATTCCATTCCATT', '@@>@BBBBBBBBBBBBAAAAAAAAAAAAAAAAA@AAAA@@@A@@@@?@9@?@@A@@@@@?@@@???????@@?>>>?>8>=<===<=<==<====<=<=6<=<==<====<=<==<=<==<=<=>====><=<==????A@@@??<>?', 'NM:i:2', 'MD:Z:24C8C42', 'AS:i:66', 'XS:i:66', 'RG:Z:0', 'SA:Z:chr22_KI270736v1_random,147858,+,53M95S,0,4;'], [433, 'chr4', '49139696', '0', '76M72H', 'chrUn_GL000216v2', '23597', '-49116105', '*', '*', 'NM:i:2', 'MD:Z:33C17G24', 'AS:i:66', 'RG:Z:0', 'SP:Z:0', 'DA:i:66', 'DP:Z:0.0', 'DN:Z:111.5', 'PS:Z:127.05', 'NP:Z:0.0'], ['385', 'chr20', '31072998', '0', '72H76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:3', 'MD:Z:41T0G8G24', 'AS:i:61', 'RG:Z:0'], ['385', 'chr17_KI270729v1_random', '27484', '0', '72H76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:3', 'MD:Z:42G8G7T16', 'AS:i:61', 'RG:Z:0'], ['401', 'chr4', '49116911', '0', '39M5D37M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:24C14^CATTG37', 'AS:i:60', 'RG:Z:0'], ['385', 'chr17', '21971340', '0', '55H93M', 'chr22', '16352731', '0', '*', '*', 'NM:i:7', 'MD:Z:8C6A43G6G1G14T0C8', 'AS:i:58', 'RG:Z:0'], ['385', 'chr17', '21974497', '0', '55H93M', 'chr22', '16352731', '0', '*', '*', 'NM:i:7', 'MD:Z:8C6A43G6G1G14T0C8', 'AS:i:58', 'RG:Z:0'], ['385', 'chr4', '49652414', '0', '72H76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:4', 'MD:Z:42G8G0C15T7', 'AS:i:56', 'RG:Z:0'], ['401', 'chrY', '11300932', '0', '76M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:4', 'MD:Z:33C14T2G12A11', 'AS:i:56', 'RG:Z:0'], ['401', 'chr20', '62049', '0', '77M71H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:1C18C3C8C16A26', 'AS:i:55', 'RG:Z:0'], ['401', 'chr4', '49099636', '0', '32M5D43M10I5M5I23M1I29M', 'chr22', '16352731', '0', '*', '*', 'NM:i:28', 'MD:Z:32^TCCAG1C17G33A6T19C0C0A17', 'AS:i:55', 'RG:Z:0'], ['385', 'chr17', '21976005', '0', '27H40M5I33M10I33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:20', 'MD:Z:10A2A12C9C27C41', 'AS:i:54', 'RG:Z:0'], ['385', 'chr10', '41883920', '0', '72H76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:21G0G4A1C29G16', 'AS:i:51', 'RG:Z:0'], ['385', 'chr17', '21972339', '0', '71H31M5D46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:8', 'MD:Z:31^AGAAA12G8T10C13', 'AS:i:51', 'RG:Z:0'], ['401', 'chr17_KI270729v1_random', '1228', '0', '76M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:24C1C1A4C9A32', 'AS:i:51', 'RG:Z:0'], ['401', 'chr4', '49098629', '0', '32M5D43M10I5M5I23M1I29M', 'chr22', '16352731', '0', '*', '*', 'NM:i:29', 'MD:Z:32^TCCAG1C17G10A22A6T19C0C0A17', 'AS:i:50', 'RG:Z:0'], ['401', 'chr4', '49099367', '0', '32M5D43M10I5M5I53M', 'chr22', '16352731', '0', '*', '*', 'NM:i:29', 'MD:Z:32^TCCAG1C17G10A22T6T17A2C0C0A17', 'AS:i:50', 'RG:Z:0'], ['401', 'chr5', '49661739', '0', '76M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:1A3T18C4G3C9G32', 'AS:i:50', 'RG:Z:0'], ['401', 'chr4', '49137646', '0', '73M75H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:12C11C26G3T9G7', 'AS:i:48', 'RG:Z:0'], ['401', 'chr4', '49142584', '0', '73M75H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:12C11C26G3T9G7', 'AS:i:48', 'RG:Z:0'], ['401', 'chr4', '49140145', '0', '73M75H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:12C11C26G3T9G7', 'AS:i:48', 'RG:Z:0'], ['385', 'chr17', '21969267', '0', '15M5I33M2D2M3D25M15D68M', 'chr22', '16352731', '0', '*', '*', 'NM:i:35', 'MD:Z:1G15G13T3A12^GG2^TGG8C6A9^AATGGAATGGAATTC33A0G8G5A18', 'AS:i:48', 'RG:Z:0'], ['401', 'chr4', '49108558', '0', '73M75H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:5T14T3C8C10T28', 'AS:i:48', 'RG:Z:0'], ['401', 'chr2', '89829255', '0', '29M5I21M20D70M23H', 'chr22', '16352731', '0', '*', '*', 'NM:i:32', 'MD:Z:46G3^ATTCCATTCCACTGAACTCA22T6G8C0C12T0T16', 'AS:i:48', 'RG:Z:0'], ['385', 'chr10', '38869401', '0', '71H77M', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:12T5C3G0G2A16G33', 'AS:i:47', 'RG:Z:0'], ['401', 'chrUn_KI270756v1', '46932', '0', '77M71H', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:33C16T2C0C3G5A12', 'AS:i:47', 'RG:Z:0'], ['385', 'chr10', '41864285', '0', '72H76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:5C4C10G2C26G5C18', 'AS:i:46', 'RG:Z:0'], ['385', 'chr2', '89838842', '0', '72H76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:31T10G7C0G8G3G11', 'AS:i:46', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '3713', '0', '76M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:33C16T0G3C3T6G9', 'AS:i:46', 'RG:Z:0'], ['385', 'chr2', '89838278', '0', '78H60M10H', 'chr22', '16352731', '0', '*', '*', 'NM:i:3', 'MD:Z:14A3C26G14', 'AS:i:45', 'RG:Z:0'], ['401', 'chr4', '49146628', '0', '36M5D40M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:9', 'MD:Z:33C2^TTCCT4T5G22A6', 'AS:i:45', 'RG:Z:0'], ['385', 'chr10', '38919075', '0', '99H49M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:15G33', 'AS:i:44', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '32941', '0', '45M10I20M10I63M', 'chr22', '16352731', '0', '*', '*', 'NM:i:31', 'MD:Z:24C0G7C0T0G37C0C12T21G1G11T4', 'AS:i:44', 'RG:Z:0'], ['401', 'chr2', '89832103', '0', '29M5I42M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:9', 'MD:Z:1A3T15A17A31', 'AS:i:44', 'RG:Z:0'], ['385', 'chr22', '16352731', '0', '62M10I20M10I46M', '=', '16352731', '0', '*', '*', 'NM:i:31', 'MD:Z:4A11C1C21A12G0G37C0A0G7C0G24', 'AS:i:44', 'RG:Z:0'], ['385', 'chr17', '21982348', '0', '69H79M', 'chr22', '16352731', '0', '*', '*', 'NM:i:8', 'MD:Z:27C1A2A12G5T2G0T3A19', 'AS:i:43', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '23557', '0', '5H40M10I33M10D37M23H', 'chr22', '16352731', '0', '*', '*', 'NM:i:27', 'MD:Z:13G4T16C33G3^CCATTCCATG12A6T11T5', 'AS:i:43', 'RG:Z:0'], ['385', 'chr4', '49642629', '0', '74H18M10I46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:11', 'MD:Z:30G33', 'AS:i:43', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '6387', '0', '77M71H', 'chr22', '16352731', '0', '*', '*', 'NM:i:7', 'MD:Z:24C1A6A0G11A2C0A26', 'AS:i:42', 'RG:Z:0'], ['385', 'chr17', '21979224', '0', '80H30M5I28M5H', 'chr22', '16352731', '0', '*', '*', 'NM:i:6', 'MD:Z:47G10', 'AS:i:42', 'RG:Z:0'], ['385', 'chr2_KI270894v1_alt', '191226', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:12G33', 'AS:i:41', 'RG:Z:0'], ['385', 'chr22_KI270736v1_random', '78879', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:12G33', 'AS:i:41', 'RG:Z:0'], ['385', 'chrY', '11739202', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:12G33', 'AS:i:41', 'RG:Z:0'], ['401', 'chr4', '49110868', '0', '25M10I41M7I2M2D40M23H', 'chr22', '16352731', '0', '*', '*', 'NM:i:25', 'MD:Z:41G10A6T8^CC8C0A12T17', 'AS:i:41', 'RG:Z:0'], ['385', 'chr17', '21896095', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:12C33', 'AS:i:41', 'RG:Z:0'], ['385', 'chrUn_KI270442v1', '90613', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:12G33', 'AS:i:41', 'RG:Z:0'], ['401', 'chrY', '11297983', '0', '46M102H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C12', 'AS:i:41', 'RG:Z:0'], ['385', 'chrY', '11723214', '0', '19H44M2D2M7I76M', 'chr22', '16352731', '0', '*', '*', 'NM:i:21', 'MD:Z:21A6T5G0G8^GG23G0G2T0A2C12G7A0G24', 'AS:i:41', 'RG:Z:0'], ['385', 'chr2', '90379579', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:12G33', 'AS:i:41', 'RG:Z:0'], ['401', 'chr2_KI270772v1_alt', '22888', '0', '46M102H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C12', 'AS:i:41', 'RG:Z:0'], ['385', 'chr20', '31059729', '0', '72H19M3D6M3I48M', 'chr22', '16352731', '0', '*', '*', 'NM:i:9', 'MD:Z:19^CTG1A16G1G33', 'AS:i:40', 'RG:Z:0'], ['385', 'chr17', '21884040', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:37A6C1', 'AS:i:39', 'RG:Z:0'], ['385', 'chr17_KI270729v1_random', '227339', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:32C11G1', 'AS:i:39', 'RG:Z:0'], ['385', 'chr17_KI270730v1_random', '27731', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:32C11G1', 'AS:i:39', 'RG:Z:0'], ['385', 'chr17', '21878885', '0', '103H38M7H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:38', 'AS:i:38', 'RG:Z:0'], ['401', 'chr21', '10713392', '0', '40M108H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:1C38', 'AS:i:38', 'RG:Z:0'], ['401', 'chr2', '89818565', '0', '43M105H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C9', 'AS:i:38', 'RG:Z:0'], ['385', 'chr9', '116362015', '0', '102H37M9H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:37', 'AS:i:37', 'RG:Z:0'], ['385', 'chr20', '31063693', '0', '72H22M4D4M4I46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:11', 'MD:Z:8T13^AACA14G1G33', 'AS:i:37', 'RG:Z:0'], ['385', 'chr10', '41916213', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:7C4G33', 'AS:i:36', 'RG:Z:0'], ['401', 'chr4', '49114496', '0', '41M107H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C7', 'AS:i:36', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '23392', '0', '41M107H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C7', 'AS:i:36', 'RG:Z:0'], ['385', 'chr17', '21907705', '0', '102H36M10H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:36', 'AS:i:36', 'RG:Z:0'], ['401', 'chr4', '49128053', '0', '20M9D50M15I39M24H', 'chr22', '16352731', '0', '*', '*', 'NM:i:32', 'MD:Z:1C18^TCTCCGGTT4C8C17G3T36T5A4A5', 'AS:i:36', 'RG:Z:0'], ['385', 'chr22', '16352413', '0', '107H41M', '=', '16352731', '319', '*', '*', 'NM:i:1', 'MD:Z:7G33', 'AS:i:36', 'RG:Z:0'], ['401', 'chrY', '10801633', '0', '45M10I20M10I5M5I53M', 'chr22', '16352731', '0', '*', '*', 'NM:i:34', 'MD:Z:14G9G8C41A6T21C0A5C0A10', 'AS:i:36', 'RG:Z:0'], ['385', 'chr20', '31166659', '0', '107H41M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:7G33', 'AS:i:36', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '65441', '0', '46M102H', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:13T0G31', 'AS:i:36', 'RG:Z:0'], ['401', 'chrY', '10778942', '0', '46M102H', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:33C5T6', 'AS:i:36', 'RG:Z:0'], ['401', 'chr4', '49135880', '0', '20M9D50M15I39M24H', 'chr22', '16352731', '0', '*', '*', 'NM:i:32', 'MD:Z:1C18^TCTCCGGTT4C8C17G3T36T5A4A5', 'AS:i:36', 'RG:Z:0'], ['401', 'chr21', '10658726', '0', '76M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:8', 'MD:Z:33C0A1G9T2C0A2C0C21', 'AS:i:36', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '33361', '0', '41M107H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C7', 'AS:i:36', 'RG:Z:0'], ['385', 'chr4', '49642607', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:9G2G33', 'AS:i:36', 'RG:Z:0'], ['385', 'chr17', '21978213', '0', '74H26M5D23M4D4M6D21M', 'chr22', '16352731', '0', '*', '*', 'NM:i:16', 'MD:Z:26^AGTGT14T8^GCAA4^ATGCAG21', 'AS:i:36', 'RG:Z:0'], ['385', 'chr10', '41893701', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:34C1A9', 'AS:i:36', 'RG:Z:0'], ['385', 'chr4', '49655365', '0', '107H41M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:7G33', 'AS:i:36', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '26764', '0', '32M5D15M1D2M6I21M72H', 'chr22', '16352731', '0', '*', '*', 'NM:i:13', 'MD:Z:32^TCAAG1C13^A23', 'AS:i:35', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '22617', '0', '40M108H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:33C6', 'AS:i:35', 'RG:Z:0'], ['385', 'chrY', '11533469', '0', '98H40M10H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:31A8', 'AS:i:35', 'RG:Z:0'], ['401', 'chrY', '10992553', '0', '49M99H', 'chr22', '16352731', '0', '*', '*', 'NM:i:3', 'MD:Z:33C2C0A11', 'AS:i:34', 'RG:Z:0'], ['385', 'chr20', '31073345', '0', '74H74M', 'chr22', '16352731', '0', '*', '*', 'NM:i:8', 'MD:Z:8A13C0A0A0A4T7G1G33', 'AS:i:34', 'RG:Z:0'], ['385', 'chr20', '31188246', '0', '104H44M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:9T0T33', 'AS:i:34', 'RG:Z:0'], ['385', 'chrY', '11523204', '0', '80H59M9H', 'chr22', '16352731', '0', '*', '*', 'NM:i:5', 'MD:Z:9T7T0G2A26A10', 'AS:i:34', 'RG:Z:0'], ['385', 'chr17', '21926947', '0', '109H39M', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:5G33', 'AS:i:34', 'RG:Z:0'], ['385', 'chrY', '11739552', '0', '102H46M', 'chr22', '16352731', '0', '*', '*', 'NM:i:3', 'MD:Z:8C2C32C1', 'AS:i:34', 'RG:Z:0'], ['401', 'chr2', '89830860', '0', '75M10I40M5I18M', 'chr22', '16352731', '0', '*', '*', 'NM:i:30', 'MD:Z:9G3T4G5C26A1G9A5A1C11C0C5A6T17A15C1', 'AS:i:34', 'RG:Z:0'], ['401', 'chr22', '10728258', '0', '39M109H', '=', '16352731', '5624436', '*', '*', 'NM:i:1', 'MD:Z:33C5', 'AS:i:34', 'RG:Z:0'], ['401', 'chr10', '42308459', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr17', '21971944', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr17', '21893850', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr2_KI270894v1_alt', '191209', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr22_KI270736v1_random', '79937', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr17', '21971835', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '38320', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10', '42295963', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr22', '16359013', '0', '115H33M', '=', '16352731', '-6283', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chrY', '11730262', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr17_KI270729v1_random', '32528', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrY', '10777535', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chrUn_KI270442v1', '95994', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr22', '10724332', '0', '33M115H', '=', '16352731', '5628368', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr22', '16346675', '0', '115H33M', '=', '16352731', '6057', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr1', '790128', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '58732', '0', '34M114H', 'chr22', '16352731', '0', '*', '*', 'NM:i:1', 'MD:Z:0T33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr21', '10655249', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10_KI270824v1_alt', '178520', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr2', '90379562', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr21', '10662430', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '18071', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr20', '31074296', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10_KI270824v1_alt', '179130', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '27844', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr17', '21974997', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr4', '49124185', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chrY', '11664057', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr14_GL000225v1_random', '199780', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr20', '31234194', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr4', '49656128', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr10', '41904931', '0', '115H33M', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrY', '10785969', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['385', 'chr22', '16347482', '0', '115H33M', '=', '16352731', '5250', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr22', '10718999', '0', '33M115H', '=', '16352731', '5633701', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr5', '49658814', '0', '51M1I23M10I44M5D19M', 'chr22', '16352731', '0', '*', '*', 'NM:i:30', 'MD:Z:16C7C8C2C9T1A2G30C0C5A6T13C2T0G3^GATTG19', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '75424', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '39127', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_KI270756v1', '1691', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10', '42318644', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10', '42296573', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '23647', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr4', '49131972', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr2_KI270772v1_alt', '22918', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrY', '10777957', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['2177', 'chr22_KI270736v1_random', '147858', '0', '53M95H', 'chr22', '16352731', '0', 'AATGGAATGGAATGGAACTAAAAGGAATGGAATGGAATGGGATGGAACGGAAT', '?><??@@@A????==<=<>====>=<=<==<=<==<=<====<==<=<6=<=<', 'NM:i:4', 'MD:Z:17T0G21A6T5', 'AS:i:33', 'XS:i:32', 'RG:Z:0', 'SA:Z:chr2,89830945,-,76M72S,0,2;'], ['401', 'chrY', '11300379', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrY', '10865463', '0', '33M115H', 'chr22', '16352731', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10', '42315909', '0', '101H47M', 'chr22', '16352731', '0', '*', '*', 'NM:i:3', 'MD:Z:6T10T9C19', 'AS:i:32', 'RG:Z:0'], ['385', 'chr17', '26785195', '0', '40M108H', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:17T0G21', 'AS:i:30', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '41573', '0', '108H40M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:17A2C19', 'AS:i:30', 'RG:Z:0'], ['385', 'chr22', '16344221', '0', '40M108H', '=', '16352731', '8511', '*', '*', 'NM:i:2', 'MD:Z:19G2T17', 'AS:i:30', 'RG:Z:0'], ['385', 'chrUn_KI270442v1', '13799', '0', '40M108H', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:6T11G21', 'AS:i:30', 'RG:Z:0'], ['401', 'chrUn_GL000216v2', '35070', '0', '108H40M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:21G0A17', 'AS:i:30', 'RG:Z:0'], ['385', 'chr22', '16344196', '0', '40M108H', '=', '16352731', '8536', '*', '*', 'NM:i:2', 'MD:Z:19G2T17', 'AS:i:30', 'RG:Z:0'], ['401', 'chr22_KI270737v1_random', '41598', '0', '108H40M', 'chr22', '16352731', '0', '*', '*', 'NM:i:2', 'MD:Z:17A2C19', 'AS:i:30', 'RG:Z:0'], ['97', 'chr4', '49153231', '0', '96M', 'chr10', '38883902', '0', 'TTGCATTCCATTCCATTCCATTGCATTCCATTTCACTCCATTCCCTTCAATTCTATTCCAGTCGAGTTGATTTCTTTCCATTCCATTCCATTCCAT', '>><>???>?@?>><<;;;<<;;>=<;;;<=;;>;<;>===<<<=>><<>>==<>><==>=?==3><>?>>==>=>=<<=><;;;;98<<?68:::9', 'NM:i:11', 'MD:Z:32C2T8A3C4C2A2C0T3G7C5T17', 'AS:i:41', 'XS:i:39', 'RG:Z:0'], ['369', 'chrY', '11690375', '0', '96M', 'chr10', '38883902', '0', '*', '*', 'NM:i:12', 'MD:Z:21T1G3G0G1A1G2A6G0C3G3T8A35', 'AS:i:39', 'RG:Z:0'], ['369', 'chr17_KI270730v1_random', '82379', '0', '56H40M', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:37G2', 'AS:i:37', 'RG:Z:0'], ['369', 'chrUn_KI270442v1', '85382', '0', '54H42M', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:6A35', 'AS:i:37', 'RG:Z:0'], ['353', 'chr21', '7921622', '0', '38M5D22M36H', 'chr10', '38883902', '0', '*', '*', 'NM:i:8', 'MD:Z:2C35^GATGT6A3G11', 'AS:i:36', 'RG:Z:0'], ['353', 'chrY', '10657312', '0', '40M56H', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:32C7', 'AS:i:35', 'RG:Z:0'], ['353', 'chr17_KI270729v1_random', '11607', '0', '4H40M52H', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:18C21', 'AS:i:35', 'RG:Z:0'], ['369', 'chr20', '31241837', '0', '96M', 'chr10', '38883902', '0', '*', '*', 'NM:i:13', 'MD:Z:18C2T1G3G0G1A1G2A6G4G3T8A2G32', 'AS:i:35', 'RG:Z:0'], ['369', 'chrUn_KI270442v1', '27148', '0', '56H35M5H', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:35', 'AS:i:35', 'RG:Z:0'], ['369', 'chr20', '31241712', '0', '96M', 'chr10', '38883902', '0', '*', '*', 'NM:i:13', 'MD:Z:3C17T1C3G0G1A1G2A6G4G3T8A2G32', 'AS:i:35', 'RG:Z:0'], ['353', 'chrY', '10776002', '0', '60M36H', 'chr10', '38883902', '0', '*', '*', 'NM:i:5', 'MD:Z:32C2T8A3C4C6', 'AS:i:35', 'RG:Z:0'], ['353', 'chr21', '10692904', '0', '44M52H', 'chr10', '38883902', '0', '*', '*', 'NM:i:2', 'MD:Z:32C2T8', 'AS:i:34', 'RG:Z:0'], ['369', 'chr20', '31054488', '0', '52H44M', 'chr10', '38883902', '0', '*', '*', 'NM:i:2', 'MD:Z:8A2G32', 'AS:i:34', 'RG:Z:0'], ['369', 'chr20', '31241642', '0', '96M', 'chr10', '38883902', '0', '*', '*', 'NM:i:13', 'MD:Z:21T1G3G0G1A1G2A2C3G4G3T8A2G32', 'AS:i:34', 'RG:Z:0'], ['369', 'chr17_KI270729v1_random', '57000', '0', '52H44M', 'chr10', '38883902', '0', '*', '*', 'NM:i:2', 'MD:Z:8A2G32', 'AS:i:34', 'RG:Z:0'], ['353', 'chr21', '10692979', '0', '96M', 'chr10', '38883902', '0', '*', '*', 'NM:i:13', 'MD:Z:32C2T8A3C3G0C6T2C1T1C0C3C1A21', 'AS:i:34', 'RG:Z:0'], ['353', 'chrUn_GL000216v2', '80442', '0', '44M52H', 'chr10', '38883902', '0', '*', '*', 'NM:i:2', 'MD:Z:32A2T8', 'AS:i:34', 'RG:Z:0'], ['353', 'chrUn_GL000216v2', '80412', '0', '44M52H', 'chr10', '38883902', '0', '*', '*', 'NM:i:2', 'MD:Z:32G2T8', 'AS:i:34', 'RG:Z:0'], ['369', 'chr10', '38517416', '0', '52H44M', '=', '38883902', '366500', '*', '*', 'NM:i:3', 'MD:Z:21G16G2G2', 'AS:i:33', 'RG:Z:0'], ['369', 'chrY', '11682277', '0', '53H43M', 'chr10', '38883902', '0', '*', '*', 'NM:i:2', 'MD:Z:7A2G32', 'AS:i:33', 'RG:Z:0'], ['369', 'chrUn_KI270442v1', '92621', '0', '58H38M', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:5G32', 'AS:i:33', 'RG:Z:0'], ['369', 'chr20', '31186086', '0', '61H35M', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:32G2', 'AS:i:32', 'RG:Z:0'], [2161, 'chrUn_KI270442v1', '92912', '0', '64H32M', 'chr4', '49139696', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0', 'SP:Z:1', 'DA:i:32', 'DP:Z:0.0', 'DN:Z:111.5', 'PS:Z:127.05', 'NP:Z:0.0'], ['369', 'chr17_KI270729v1_random', '57032', '0', '64H32M', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['353', 'chr5', '49658420', '0', '31M5I60M', 'chr10', '38883902', '0', '*', '*', 'NM:i:15', 'MD:Z:39T3C4C6A3G7C1G3T4T0T11', 'AS:i:32', 'RG:Z:0'], ['353', 'chrY', '10979355', '0', '32M64H', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['369', 'chr20', '31059847', '0', '64H32M', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['369', 'chr1', '224014161', '0', '64H32M', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['353', 'chrY', '11294605', '0', '32M64H', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['369', 'chrY', '11664199', '0', '61H35M', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:32G2', 'AS:i:32', 'RG:Z:0'], ['369', 'chr10', '38579142', '0', '36H57M3H', '=', '38883902', '304761', '*', '*', 'NM:i:5', 'MD:Z:24A2G6C2G1C17', 'AS:i:32', 'RG:Z:0'], ['369', 'chr17_KI270729v1_random', '57052', '0', '64H32M', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['369', 'chr10', '41909243', '0', '64H32M', '=', '38883902', '-3025317', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['353', 'chr21', '10671877', '0', '3H32M61H', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['369', 'chr20', '31241542', '0', '26M10I60M', 'chr10', '38883902', '0', '*', '*', 'NM:i:18', 'MD:Z:21T1C8G4G3T6C1A2C32', 'AS:i:32', 'RG:Z:0'], ['369', 'chrUn_KI270438v1', '70324', '0', '63H33M', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:32T0', 'AS:i:32', 'RG:Z:0'], ['353', 'chr4', '49132517', '0', '32M64H', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['369', 'chr20', '31073466', '0', '64H32M', 'chr10', '38883902', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['353', 'chr22', '10722703', '0', '44M52H', 'chr10', '38883902', '0', '*', '*', 'NM:i:3', 'MD:Z:2C17A1C21', 'AS:i:31', 'RG:Z:0'], ['369', 'chrUn_KI270757v1', '16310', '0', '52H44M', 'chr10', '38883902', '0', '*', '*', 'NM:i:3', 'MD:Z:21G1G17G2', 'AS:i:31', 'RG:Z:0'], ['369', 'chr22_KI270736v1_random', '66757', '0', '54H42M', 'chr10', '38883902', '0', '*', '*', 'NM:i:3', 'MD:Z:24A3A10G2', 'AS:i:30', 'RG:Z:0'], ['369', 'chr4', '49651308', '0', '50M2D46M', 'chr10', '38883902', '0', '*', '*', 'NM:i:14', 'MD:Z:6A24C3A0G2T2G4G2^AT1T8A2G9G19G2', 'AS:i:30', 'RG:Z:0'], ['353', 'chrY', '10937109', '0', '9H35M52H', 'chr10', '38883902', '0', '*', '*', 'NM:i:1', 'MD:Z:13A21', 'AS:i:30', 'RG:Z:0'], ['145', 'chr10', '38883902', '0', '16S57M75S', 'chr4', '49153231', '0', 'GAGTGGAATGCAATGGAATGGAATCAAACGGATTGGAATTGAAAGGAAGGGAATGGAATGGAATGGAATGGAAAGAAATCAACTCGACTGGAATAGAATTGAAGGGAATGGAGTGAAATGGAATGCAATGGAATGGAATGGAATGCAA', '@?*@BBCBC>7CCCBBCBB@@@@?/6@>5A@>8>A?@A>8?A@A<A?B=@@89:::86?<<89;;;;<>=<<=>=>==>>?><>3==?=>==<>><==>><<>>=<<<===>;<;>;;=<;;;<=>;;<<;;;<<>>?@?>???><>>', 'NM:i:1', 'MD:Z:36G20', 'AS:i:52', 'XS:i:51', 'RG:Z:0', 'SA:Z:chr4,49153231,+,96M52S,0,11;'], ['385', 'chrUn_KI270756v1', '32334', '0', '76H56M16H', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:19C36', 'AS:i:51', 'RG:Z:0'], ['385', 'chr10', '38784647', '0', '77H55M16H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:18C21A14', 'AS:i:45', 'RG:Z:0'], ['385', 'chr5', '49658420', '0', '31M5I112M', 'chr4', '49153231', '0', '*', '*', 'NM:i:23', 'MD:Z:39T3C4C6A3G7C1G3T4T0T14A4C3C6C2G1G13T2C10', 'AS:i:42', 'RG:Z:0'], ['401', 'chr2', '89836223', '0', '34H114M', 'chr4', '49153231', '0', '*', '*', 'NM:i:15', 'MD:Z:5G33T1G3G3C1A1G6G4G0C2T8A2G9G19G2', 'AS:i:41', 'RG:Z:0'], ['2177', 'chr4', '49153231', '0', '96M52H', '=', '49153231', '0', 'TTGCATTCCATTCCATTCCATTGCATTCCATTTCACTCCATTCCCTTCAATTCTATTCCAGTCGAGTTGATTTCTTTCCATTCCATTCCATTCCAT', '>><>???>?@?>><<;;;<<;;>=<;;;<=;;>;<;>===<<<=>><<>>==<>><==>=?==3><>?>>==>=>=<<=><;;;;98<<?68:::9', 'NM:i:11', 'MD:Z:32C2T8A3C4C2A2C0T3G7C5T17', 'AS:i:41', 'XS:i:37', 'RG:Z:0', 'SA:Z:chr10,38883902,-,16S57M75S,0,1;'], ['401', 'chrY', '11690328', '0', '23M10I63M5D52M', 'chr4', '49153231', '0', '*', '*', 'NM:i:30', 'MD:Z:2A6A16T11T24T1G3G0G1A1G2A6G1^AATGC3G3T8A35', 'AS:i:40', 'RG:Z:0'], ['401', 'chr17', '21867006', '0', '16H59M73H', 'chr4', '49153231', '0', '*', '*', 'NM:i:4', 'MD:Z:9G11G1C3T31', 'AS:i:39', 'RG:Z:0'], ['401', 'chr20', '31241670', '0', '23M10I115M', 'chr4', '49153231', '0', '*', '*', 'NM:i:27', 'MD:Z:2A26G3T4T6C17T1C3G0G1A1G2A6G4G3T8A2G32', 'AS:i:39', 'RG:Z:0'], ['401', 'chr20', '31241795', '0', '73M10I65M', 'chr4', '49153231', '0', '*', '*', 'NM:i:27', 'MD:Z:2A21G0G2T3A2C3G3T4T11C13G2A6G4G3T8A2G32', 'AS:i:39', 'RG:Z:0'], ['401', 'chr20', '31241615', '0', '23M10I40M15I60M', 'chr4', '49153231', '0', '*', '*', 'NM:i:35', 'MD:Z:2A26G3T4T26C3G4G3T8A2G32', 'AS:i:39', 'RG:Z:0'], ['401', 'chr10', '38848228', '0', '16H52M80H', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:9G22C5T13', 'AS:i:37', 'RG:Z:0'], ['401', 'chr17_KI270730v1_random', '82379', '0', '108H40M', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:37G2', 'AS:i:37', 'RG:Z:0'], ['385', 'chrUn_KI270756v1', '68100', '0', '80H52M16H', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:13A5G22C9', 'AS:i:37', 'RG:Z:0'], [2225, 'chrUn_KI270442v1', '85382', '0', '106H42M', 'chrUn_GL000216v2', '23597', '0', '*', '*', 'NM:i:1', 'MD:Z:6A35', 'AS:i:37', 'RG:Z:0', 'SP:Z:0', 'DA:i:37', 'DP:Z:0.0', 'DN:Z:111.5', 'PS:Z:127.05', 'NP:Z:0.0'], ['385', 'chrY', '11302715', '0', '73H75M', 'chr4', '49153231', '0', '*', '*', 'NM:i:8', 'MD:Z:31A3C6C2G1A0A15C7T2', 'AS:i:37', 'RG:Z:0'], ['385', 'chrY', '10845623', '0', '60M15I40M33H', 'chr4', '49153231', '0', '*', '*', 'NM:i:24', 'MD:Z:2C16T2C9C2T12C4C35A3C6', 'AS:i:36', 'RG:Z:0'], ['385', 'chr21', '7921622', '0', '38M5D22M88H', 'chr4', '49153231', '0', '*', '*', 'NM:i:8', 'MD:Z:2C35^GATGT6A3G11', 'AS:i:36', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '27148', '0', '108H35M5H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:35', 'AS:i:35', 'RG:Z:0'], ['385', 'chrY', '10657312', '0', '40M108H', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:32C7', 'AS:i:35', 'RG:Z:0'], ['385', 'chr17_KI270729v1_random', '11607', '0', '4H40M10I94M', 'chr4', '49153231', '0', '*', '*', 'NM:i:27', 'MD:Z:18C24A2C3G4T2C0A5C3T15G4A3C6T3A2C0C13C10', 'AS:i:35', 'RG:Z:0'], ['385', 'chrY', '10776002', '0', '60M88H', 'chr4', '49153231', '0', '*', '*', 'NM:i:5', 'MD:Z:32C2T8A3C4C6', 'AS:i:35', 'RG:Z:0'], ['385', 'chrY', '10845641', '0', '73H46M5I20M4H', 'chr4', '49153231', '0', '*', '*', 'NM:i:9', 'MD:Z:31A3C6T16C6', 'AS:i:35', 'RG:Z:0'], ['401', 'chrY', '11568227', '0', '28M5I40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:10', 'MD:Z:2A7A13G24C0A17', 'AS:i:34', 'RG:Z:0'], ['401', 'chrY', '11579436', '0', '28M5I40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:10', 'MD:Z:2A7A13G24C0A17', 'AS:i:34', 'RG:Z:0'], ['401', 'chr17', '21970278', '0', '33H79M36H', 'chr4', '49153231', '0', '*', '*', 'NM:i:9', 'MD:Z:6G3T29T1G7C3G6G4G3T8', 'AS:i:34', 'RG:Z:0'], ['385', 'chr21', '10692904', '0', '59M10I79M', 'chr4', '49153231', '0', '*', '*', 'NM:i:28', 'MD:Z:32C2T8A3C3G0C8C1A12G11A4A2G0C6T3A2C0C21T2', 'AS:i:34', 'RG:Z:0'], ['401', 'chr20', '31054488', '0', '104H44M', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:8A2G32', 'AS:i:34', 'RG:Z:0'], ['385', 'chrUn_GL000216v2', '80412', '0', '44M104H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:32G2T8', 'AS:i:34', 'RG:Z:0'], ['385', 'chrUn_GL000216v2', '80442', '0', '44M104H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:32A2T8', 'AS:i:34', 'RG:Z:0'], ['401', 'chr17_KI270729v1_random', '57000', '0', '104H44M', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:8A2G32', 'AS:i:34', 'RG:Z:0'], ['401', 'chrY', '11573829', '0', '28M5I40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:10', 'MD:Z:2A7A13G24C0A17', 'AS:i:34', 'RG:Z:0'], ['401', 'chrY', '11590644', '0', '28M5I40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:10', 'MD:Z:2A7A13G24C0A17', 'AS:i:34', 'RG:Z:0'], ['401', 'chrY', '11551377', '0', '28M5I40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:10', 'MD:Z:2A7A13G24C0A17', 'AS:i:34', 'RG:Z:0'], ['401', 'chr22_KI270736v1_random', '80369', '0', '48H25M15I49M11H', 'chr4', '49153231', '0', '*', '*', 'NM:i:19', 'MD:Z:31G4G3T21A11', 'AS:i:33', 'RG:Z:0'], ['401', 'chr4', '49640102', '0', '24M5I40M79H', '=', '49153231', '-486935', '*', '*', 'NM:i:9', 'MD:Z:27A6G3T4T20', 'AS:i:33', 'RG:Z:0'], ['401', 'chr10', '38517416', '0', '104H44M', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:21G16G2G2', 'AS:i:33', 'RG:Z:0'], ['385', 'chrY', '10856519', '0', '75H33M40H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:33', 'AS:i:33', 'RG:Z:0'], ['401', 'chrY', '11682277', '0', '105H43M', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:7A2G32', 'AS:i:33', 'RG:Z:0'], ['385', 'chr4', '49150436', '0', '75H73M', '=', '49153231', '2796', '*', '*', 'NM:i:8', 'MD:Z:14G9A4A3C6T3A2C0C24', 'AS:i:33', 'RG:Z:0'], ['401', 'chr20', '31072184', '0', '44H27M5I7M10I18M5I32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:23', 'MD:Z:40G3A36G2', 'AS:i:33', 'RG:Z:0'], ['385', 'chr4', '49092464', '0', '3M1I56M15I33M40H', '=', '49153231', '60768', '*', '*', 'NM:i:23', 'MD:Z:7T11C1C21A3C4C30G8', 'AS:i:33', 'RG:Z:0'], ['385', 'chr10', '38785493', '0', '94H43M11H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:5A8G28', 'AS:i:33', 'RG:Z:0'], ['401', 'chr20', '31051549', '0', '39M10I41M5I53M', 'chr4', '49153231', '0', '*', '*', 'NM:i:30', 'MD:Z:2A25G0C2G17A12T1G7A3G6G12A3G4G3G14A7', 'AS:i:33', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '92621', '0', '110H38M', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:5G32', 'AS:i:33', 'RG:Z:0'], ['401', 'chr17_KI270729v1_random', '57052', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['385', 'chrY', '11294605', '0', '32M116H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr17_KI270729v1_random', '57032', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr10', '41909243', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr1', '224014161', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['385', 'chrY', '10687556', '0', '69H50M5I24M', 'chr4', '49153231', '0', '*', '*', 'NM:i:12', 'MD:Z:30A4A3C6T11T4C9G0', 'AS:i:32', 'RG:Z:0'], ['401', 'chrY', '11664199', '0', '113H35M', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:32G2', 'AS:i:32', 'RG:Z:0'], ['401', 'chr10', '41911695', '0', '33H42M73H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:6G3T31', 'AS:i:32', 'RG:Z:0'], ['401', 'chr20', '31073466', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '92912', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr20', '31241596', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['385', 'chr21', '10671877', '0', '3H32M113H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr10', '38579142', '0', '88H57M3H', 'chr4', '49153231', '0', '*', '*', 'NM:i:5', 'MD:Z:24A2G6C2G1C17', 'AS:i:32', 'RG:Z:0'], ['401', 'chrUn_KI270438v1', '70324', '0', '115H33M', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:32T0', 'AS:i:32', 'RG:Z:0'], ['401', 'chr17', '21916263', '0', '16H32M100H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr20', '31186086', '0', '113H35M', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:32G2', 'AS:i:32', 'RG:Z:0'], ['401', 'chr20', '31241042', '0', '40H37M71H', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:28C8', 'AS:i:32', 'RG:Z:0'], ['401', 'chr20', '31059847', '0', '116H32M', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['385', 'chr4', '49132517', '0', '32M116H', '=', '49153231', '20715', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['385', 'chrY', '10979355', '0', '32M116H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:32', 'AS:i:32', 'RG:Z:0'], ['401', 'chr10', '38578982', '0', '33H42M73H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:6G3T31', 'AS:i:32', 'RG:Z:0'], ['401', 'chr22_KI270736v1_random', '59938', '0', '11H13M5I44M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:8', 'MD:Z:21G1G3T29', 'AS:i:31', 'RG:Z:0'], ['401', 'chrUn_KI270757v1', '16310', '0', '104H44M', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:21G1G17G2', 'AS:i:31', 'RG:Z:0'], ['401', 'chr5', '49602668', '0', '24M5I12M5D32M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:13', 'MD:Z:34G1^AATGC2T14G14', 'AS:i:31', 'RG:Z:0'], ['385', 'chrUn_KI270756v1', '39204', '0', '75H33M5D24M16H', 'chr4', '49153231', '0', '*', '*', 'NM:i:8', 'MD:Z:29C3^CGTAT0G13A9', 'AS:i:31', 'RG:Z:0'], ['385', 'chr22', '10722703', '0', '44M104H', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:2C17A1C21', 'AS:i:31', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '84572', '0', '44M5D29M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:12', 'MD:Z:0C9G12G0G0G6A6G4^GGAAT29', 'AS:i:31', 'RG:Z:0'], ['401', 'chr20', '31055928', '0', '39M10I41M5I53M', 'chr4', '49153231', '0', '*', '*', 'NM:i:31', 'MD:Z:2A25G0C2G17A12T1G7A3G6G12A2G0G4G3G14A7', 'AS:i:31', 'RG:Z:0'], ['401', 'chr20', '31066852', '0', '69M79H', 'chr4', '49153231', '0', '*', '*', 'NM:i:8', 'MD:Z:2A25G0C2G6G3T4T1A18', 'AS:i:31', 'RG:Z:0'], ['401', 'chr10', '38877036', '0', '16H23M5D34M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:8', 'MD:Z:9T13^CATAC0G3G29', 'AS:i:31', 'RG:Z:0'], ['401', 'chr10', '41908203', '0', '33H46M69H', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:6G8A0T29', 'AS:i:31', 'RG:Z:0'], ['401', 'chr2', '89839907', '0', '18M10D6M5I44M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:18', 'MD:Z:18^TCAACCTGAG9A6G3T29', 'AS:i:30', 'RG:Z:0'], ['401', 'chr4', '49636550', '0', '33H40M75H', '=', '49153231', '-483359', '*', '*', 'NM:i:2', 'MD:Z:6G3T29', 'AS:i:30', 'RG:Z:0'], ['401', 'chr22_KI270736v1_random', '66757', '0', '106H42M', 'chr4', '49153231', '0', '*', '*', 'NM:i:3', 'MD:Z:24A3A10G2', 'AS:i:30', 'RG:Z:0'], ['385', 'chrY', '10869929', '0', '70H49M5I24M', 'chr4', '49153231', '0', '*', '*', 'NM:i:12', 'MD:Z:29A4A1C1C6T16C7T2', 'AS:i:30', 'RG:Z:0'], ['385', 'chrUn_KI270756v1', '79274', '0', '69H30M49H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:30', 'AS:i:30', 'RG:Z:0'], ['385', 'chrUn_GL000216v2', '1825', '0', '53M1I5M10I46M33H', 'chr4', '49153231', '0', '*', '*', 'NM:i:22', 'MD:Z:2C19C8G3T8A3G1C10C1A29A3C6', 'AS:i:30', 'RG:Z:0'], ['385', 'chrUn_GL000216v2', '12132', '0', '75H40M33H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:29A3C6', 'AS:i:30', 'RG:Z:0'], ['385', 'chr4', '49118158', '0', '73H30M45H', '=', '49153231', '35074', '*', '*', 'NM:i:0', 'MD:Z:30', 'AS:i:30', 'RG:Z:0'], ['401', 'chr10', '41889869', '0', '49H30M69H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:30', 'AS:i:30', 'RG:Z:0'], ['385', 'chr22', '10718524', '0', '69H30M49H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:30', 'AS:i:30', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '85925', '0', '33H40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:6G3T29', 'AS:i:30', 'RG:Z:0'], ['385', 'chrY', '10820966', '0', '69H35M44H', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:5A29', 'AS:i:30', 'RG:Z:0'], ['401', 'chr10', '38837076', '0', '49H30M69H', 'chr4', '49153231', '0', '*', '*', 'NM:i:0', 'MD:Z:30', 'AS:i:30', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '88137', '0', '33H40M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:2', 'MD:Z:6G3T29', 'AS:i:30', 'RG:Z:0'], ['401', 'chr17', '21969351', '0', '73M75H', 'chr4', '49153231', '0', '*', '*', 'NM:i:9', 'MD:Z:2A6T19C2G6G3T0A3T11A12', 'AS:i:30', 'RG:Z:0'], ['401', 'chr4', '49651308', '0', '52H50M2D46M', '=', '49153231', '-498175', '*', '*', 'NM:i:14', 'MD:Z:6A24C3A0G2T2G4G2^AT1T8A2G9G19G2', 'AS:i:30', 'RG:Z:0'], ['401', 'chrUn_KI270442v1', '84661', '0', '44H35M69H', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:29T5', 'AS:i:30', 'RG:Z:0'], ['385', 'chrY', '10937109', '0', '9H35M104H', 'chr4', '49153231', '0', '*', '*', 'NM:i:1', 'MD:Z:13A21', 'AS:i:30', 'RG:Z:0']], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'rows': array([ 11, 210, 175,   2]), 'read1_q': '>><>???>?@?>><<;;;<<;;>=<;;;<=;;>;<;>===<<<=>><<>>==<>><==>=?==3><>?>>==>=>=<<=><;;;;98<<?68:::9', 'read1_length': 96, 'data': array([]), 'name': 'HISEQ1:11:H8GV6ADXX:2:1104:4363:8803', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 'TTGCATTCCATTCCATTCCATTGCATTCCATTTCACTCCATTCCCTTCAATTCTATTCCAGTCGAGTTGATTTCTTTCCATTCCATTCCATTCCAT', 'last_seen_chrom': 'chrY', 'ri': {0.0: 2, 1.0: 3, 2.0: 52, 3.0: 18, 4.0: 22, 5.0: 40, 6.0: 41, 7.0: 42, 8.0: 43, 9.0: 44, 10.0: 53, 11.0: 54, 12.0: 55, 13.0: 56, 14.0: 57, 15.0: 58, 16.0: 59, 17.0: 60, 18.0: 61, 19.0: 62, 20.0: 63, 21.0: 64, 22.0: 65, 23.0: 66, 24.0: 67, 25.0: 68, 26.0: 69, 27.0: 70, 28.0: 71, 29.0: 72, 30.0: 73, 31.0: 74, 32.0: 75, 33.0: 76, 34.0: 77, 35.0: 78, 36.0: 79, 37.0: 80, 38.0: 81, 39.0: 209, 40.0: 82, 41.0: 83, 42.0: 84, 43.0: 85, 44.0: 86, 45.0: 87, 46.0: 204, 47.0: 88, 48.0: 89, 49.0: 205, 50.0: 91, 51.0: 92, 52.0: 93, 53.0: 94, 54.0: 95, 55.0: 96, 56.0: 97, 57.0: 98, 58.0: 99, 59.0: 100, 60.0: 101, 61.0: 103, 62.0: 104, 63.0: 105, 64.0: 206, 65.0: 109, 66.0: 110, 67.0: 207, 68.0: 111, 69.0: 113, 70.0: 114, 71.0: 115, 72.0: 210, 73.0: 116, 74.0: 117, 75.0: 118, 76.0: 119, 77.0: 120, 78.0: 121, 79.0: 122, 80.0: 123, 81.0: 124, 82.0: 125, 83.0: 126, 84.0: 127, 85.0: 128, 86.0: 129, 87.0: 130, 88.0: 211, 89.0: 132, 90.0: 133, 91.0: 134, 92.0: 208, 93.0: 135, 94.0: 136, 95.0: 137, 96.0: 138, 97.0: 145, 98.0: 146, 99.0: 147, 100.0: 148, 101.0: 149, 102.0: 150, 103.0: 151, 104.0: 152, 105.0: 153, 106.0: 154, 107.0: 155, 108.0: 156, 109.0: 157, 110.0: 158, 111.0: 159, 112.0: 160, 113.0: 161, 114.0: 162, 115.0: 163, 116.0: 164, 117.0: 165, 118.0: 166, 119.0: 167, 120.0: 168, 121.0: 169, 122.0: 170, 123.0: 171, 124.0: 172, 125.0: 173, 126.0: 174, 127.0: 175, 128.0: 176, 129.0: 177, 130.0: 178, 131.0: 179, 132.0: 180, 133.0: 181, 134.0: 182, 135.0: 183, 136.0: 184, 137.0: 185, 138.0: 186, 139.0: 187, 140.0: 188, 141.0: 189, 142.0: 251, 143.0: 190, 144.0: 191, 145.0: 252, 146.0: 267, 147.0: 268, 148.0: 269, 149.0: 270, 150.0: 271, 151.0: 272, 152.0: 273, 153.0: 0, 154.0: 1, 155.0: 4, 156.0: 5, 157.0: 6, 158.0: 7, 159.0: 49, 160.0: 8, 161.0: 50, 162.0: 9, 163.0: 10, 164.0: 11, 165.0: 12, 166.0: 13, 167.0: 14, 168.0: 15, 169.0: 16, 170.0: 17, 171.0: 19, 172.0: 20, 173.0: 21, 174.0: 23, 175.0: 24, 176.0: 25, 177.0: 26, 178.0: 27, 179.0: 28, 180.0: 29, 181.0: 30, 182.0: 31, 183.0: 47, 184.0: 32, 185.0: 33, 186.0: 48, 187.0: 34, 188.0: 35, 189.0: 36, 190.0: 37, 191.0: 38, 192.0: 39, 193.0: 45, 194.0: 46, 195.0: 51, 196.0: 214, 197.0: 215, 198.0: 216, 199.0: 90, 200.0: 233, 201.0: 246, 202.0: 102, 203.0: 217, 204.0: 106, 205.0: 107, 206.0: 108, 207.0: 218, 208.0: 264, 209.0: 219, 210.0: 262, 211.0: 112, 212.0: 223, 213.0: 248, 214.0: 265, 215.0: 266, 216.0: 131, 217.0: 249, 218.0: 203, 219.0: 139, 220.0: 140, 221.0: 224, 222.0: 141, 223.0: 253, 224.0: 254, 225.0: 255, 226.0: 256, 227.0: 142, 228.0: 143, 229.0: 144, 230.0: 241, 231.0: 192, 232.0: 257, 233.0: 234, 234.0: 261, 235.0: 193, 236.0: 237, 237.0: 235, 238.0: 212, 239.0: 194, 240.0: 274, 241.0: 279, 242.0: 280, 243.0: 281, 244.0: 282, 245.0: 283, 246.0: 195, 247.0: 275, 248.0: 225, 249.0: 284, 250.0: 285, 251.0: 286, 252.0: 276, 253.0: 250, 254.0: 278, 255.0: 220, 256.0: 277, 257.0: 236, 258.0: 287, 259.0: 288, 260.0: 289, 261.0: 226, 262.0: 213, 263.0: 258, 264.0: 196, 265.0: 221, 266.0: 259, 267.0: 197, 268.0: 198, 269.0: 199, 270.0: 222, 271.0: 227, 272.0: 200, 273.0: 228, 274.0: 263, 275.0: 201, 276.0: 242, 277.0: 229, 278.0: 230, 279.0: 240, 280.0: 243, 281.0: 244, 282.0: 231, 283.0: 238, 284.0: 245, 285.0: 232, 286.0: 202, 287.0: 247, 288.0: 239, 289.0: 260}, 'score_mat': {'splitter': [1, 0], 'chr4-49139696--1-2': [True, 66.0], 'dis_to_next_path': 0.0, 'dis_to_normal': 111.5, 'path_score': 127.05000305175781, 'chrUn_KI270442v1-92912--1-1': [False, 32.0], 'normal_pairings': 0, 'chrUn_GL000216v2-23597--1-1': [True, 0.0], 'chrUn_KI270442v1-85382--1-2': [False, 37.0]}, 'pairing_params': (150.0, 17.0, 150.0, 0.05, 1.05, 2.0, 9.0)}
    sam = fixsam(template)
    print type(sam)
    print len(sam)
    for l in sam:
        print l

