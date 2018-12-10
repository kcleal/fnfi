"""
Utils to generate proper sam output and flag information
"""

import re
import click
from c_io_funcs import reverse_complement as rev_comp
from c_samflags import set_bit


def echo(*arg):
    click.echo(arg, err=True)


# def rev_comp(s):
#     r = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "T", "c": "G", "g": "C", "t": "A", "n": "N"}
#     return "".join(r[i] for i in s)[::-1]

#
# def set_bit(v, index, x):
#     """Set the index:th bit of v to 1 if x is truthy, else to 0, and return the new value."""
#     mask = 1 << index  # Compute mask, an integer with just bit 'index' set.
#     v &= ~mask  # Clear the bit indicated by the mask (if x is False)
#     if x:
#         v |= mask  # If x was True, set the bit indicated by the mask.
#
#     assert v == set_bit2(v, index, x)
#     return v


def add_softclips(s, fq, readid):
    fq_seq = fq[readid][0].upper()
    fq_qual = fq[readid][1]
    if s.flag & 16:
        fq_seq = rev_comp(fq_seq, len(fq_seq))
        fq_qual = fq_qual[::-1]
    s.seq = fq_seq
    s.qual = fq_qual


def set_mate_flag(a, b, r1_len, r2_len, max_d, read1_rev, read2_rev, tempplate):

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

    # Turn off secondary flag as these are primary pair
    # aflag = set_bit(aflag, 8, 0)
    # bflag = set_bit(bflag, 8, 0)

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
            if p1 < p2:
                tl = p2 + r2_len - p1
                a[7] = str(tl)
                b[7] = str(-1 * tl)
            else:
                tl = p2 - p1 - r1_len
                a[7] = str(tl)
                b[7] = str(-1 * tl)

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

    # Set paird and supplementary flag
    if not supflag & 1:
        supflag = set_bit(supflag, 0, 1)
    if not supflag & 2048:
        supflag = set_bit(supflag, 11, 1)

    # If primary is on reverse strand, set the mate reverse strand tag
    if priflag & 16 and not supflag & 32:
        supflag = set_bit(supflag, 5, 1)

    # Turn off not-primary-alignment
    if supflag & 256:
        supflag = set_bit(supflag, 8, 0)

    # Sometimes supplementary needs to be reverse complemented too
    # Occurs when primary has been rev-comp and supp is forward strand or vice versa
    # if primary_reversed and not supflag & 16:
    #     rev_sup = True

    if not ori_primary_reversed and supflag & 16:
        rev_sup = True
    elif ori_primary_reversed and not supflag & 16 and not primary_will_be_reversed:
        rev_sup = True
    else:
        rev_sup = False

    sup[0] = supflag

    sup[5] = pri[1]
    sup[6] = pri[2]

    return rev_sup


def add_sequence_back(item, reverse_me, template):
    # item is the alignment
    # if template["name"] == u'HISEQ1:9:H8962ADXX:1:1101:4291:82991':
    #     echo(template)
    #     echo()

    c = re.split(r'(\d+)', item[4])[1:]  # Drop leading empty string
    start = 0

    cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "DH"])

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

        if reverse_me:
            item[8] = rev_comp(s, len(s))
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
            item[8] = rev_comp(template[sqn], len(template[sqn]))[start:end]
            item[9] = template["fq_%s_q" % name][::-1][start:end]
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
        # echo(template)

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
            out.append(['sup', l, False])  # Supplementary

    if primary1 is 0 or primary2 is 0 and template["paired_end"]:
        return []  # Todo deal with unmapped read or unpaired

    if paired and template["paired_end"]:

        rev_A, rev_B = set_mate_flag(primary1, primary2, r1l, r2l, max_d,
                                     template["read1_reverse"], template["read2_reverse"], template)

        # Check if supplementary needs reverse complementing
        for i in range(len(out)):
            if out[i][1][0] & 64:  # Flag
                revsup = set_supp_flags(out[i][1], primary1, template["read1_reverse"], rev_A)
            else:
                revsup = set_supp_flags(out[i][1], primary2, template["read2_reverse"], rev_B)

            if revsup:
                out[i][2] = True

    out = [('pri', primary1, rev_A), ('pri', primary2, rev_B)] + out

    # Add read seq info back in if necessary, before reverse complementing. Check for hard clips and clip as necessary
    for a_type, aln, reverse_me in out:

        if aln:

            if len(aln[8]) <= 1 or "H" in aln[4]:  # Sequence is set as "*", needs adding back in, or hard-clipped
                add_sequence_back(aln, reverse_me, template)
            elif reverse_me:
                aln[8] = rev_comp(aln[8], len(aln[8]))
                aln[9] = aln[9][::-1]

            # Turn off not primary here
            aln[0] = set_bit(aln[0], 8, 0)

            # Convert flag back to string
            #aln[0] = str(aln[0])
        else:
            # None here means no alignment for primary2
            pass

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

    template = {'splitters': [False, False], 'paired_end': 1, 'locs': ['chr17-116696-1-1', 'chr17-116846--1-2'], 'bias': 1.15, 'fq_read2_seq': 0, 'isize': (210.0, 175.0), 'read2_q': 'FCGBFGDFFDAA6FFGGEG<FFGGGGGGGGD$7GG<FFC:G?GG;DAEEB->AFGG@CDFGEFGFGG@FBGEGFGFGGFFFEFGGGGGGGGFGCGF&.(BDFGDEGGGGGEEGFGFGGCEGF>FG', 'max_d': 910.0, 'read2_seq': 'ATTAGTCTTAGGGCATGGTGTGGATATTATTACATTAGTATTGGAAGCGATGGTGTGGACTAGATCAGTGATAGGGCATGGTGTGGATATTATTACCTTAGTATTGGAAGCGATGGTGTGGACTA', 'chrom_ids': {'chr17_GL383563v3_alt': 1, 'chr11': 2, 'chr17': 0}, 'read2_length': 125, 'passed': True, 'replace_hard': 0, 'read2_reverse': 1, 'inputdata': [['97', 'chr17', '116696', '0', '125M', '=', '116846', '275', 'ACAGTCCACGTAGTAGATATTATTATATTTATTAGTATTAGAGGCCATTGTGTCAATATATTAGTGTTAGGGCATGGTATGGGTGTCTAGATTAGTGTTAGGGAATGGTGTGGATATTGTATTGG', 'F#GGGGEFGBG:GGG=GEGAGDGD@EFGGDFGFF;FBBGGF4FAFGDG3.&FGCGG@DG::GGCGGGFGFGGGGD=EGGGFGF<GEG;G?C3EGFG<GEFG9GGD>-D?7=7DGGF5<@EAFFGC', 'NM:i:1', 'MD:Z:1G123', 'AS:i:123', 'XS:i:123', 'RG:Z:0'], ['353', 'chr17_GL383563v3_alt', '56696', '0', '125M', 'chr17', '116846', '0', '*', '*', 'NM:i:1', 'MD:Z:1G123', 'AS:i:123', 'RG:Z:0'], ['353', 'chr11', '179656', '0', '125M', 'chr17', '116846', '0', '*', '*', 'NM:i:3', 'MD:Z:1G80A3T38', 'AS:i:113', 'RG:Z:0'], ['145', 'chr17', '116846', '0', '125M', '=', '116696', '-275', 'ATTAGTCTTAGGGCATGGTGTGGATATTATTACATTAGTATTGGAAGCGATGGTGTGGACTAGATCAGTGATAGGGCATGGTGTGGATATTATTACCTTAGTATTGGAAGCGATGGTGTGGACTA', 'FCGBFGDFFDAA6FFGGEG<FFGGGGGGGGD$7GG<FFC:G?GG;DAEEB->AFGG@CDFGEFGFGG@FBGEGFGFGGFFFEFGGGGGGGGFGCGF&.(BDFGDEGGGGGEEGFGFGGCEGF>FG', 'NM:i:1', 'MD:Z:96A28', 'AS:i:120', 'XS:i:120', 'RG:Z:0'], ['401', 'chr17_GL383563v3_alt', '56846', '0', '125M', 'chr17', '116696', '0', '*', '*', 'NM:i:1', 'MD:Z:96A28', 'AS:i:120', 'RG:Z:0'], ['401', 'chr17', '117161', '0', '125M', '=', '116696', '-590', '*', '*', 'NM:i:3', 'MD:Z:2C3G89A28', 'AS:i:113', 'RG:Z:0'], ['401', 'chr17_GL383563v3_alt', '57161', '0', '125M', 'chr17', '116696', '0', '*', '*', 'NM:i:3', 'MD:Z:2C3G89A28', 'AS:i:113', 'RG:Z:0'], ['401', 'chr17', '117035', '0', '125M', '=', '116696', '-464', '*', '*', 'NM:i:5', 'MD:Z:2C3G89A9A15T2', 'AS:i:105', 'RG:Z:0'], ['401', 'chr17_GL383563v3_alt', '57035', '0', '125M', 'chr17', '116696', '0', '*', '*', 'NM:i:5', 'MD:Z:2C3G89A9A15T2', 'AS:i:105', 'RG:Z:0'], ['401', 'chr17_GL383563v3_alt', '57106', '0', '8H117M', 'chr17', '116696', '0', '*', '*', 'NM:i:4', 'MD:Z:35A15T10T25A28', 'AS:i:97', 'RG:Z:0'], ['401', 'chr17', '117106', '0', '8H117M', '=', '116696', '-527', '*', '*', 'NM:i:4', 'MD:Z:35A15T10T25A28', 'AS:i:97', 'RG:Z:0']], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'rows': array([0, 3]), 'read1_q': 'F#GGGGEFGBG:GGG=GEGAGDGD@EFGGDFGFF;FBBGGF4FAFGDG3.&FGCGG@DG::GGCGGGFGFGGGGD=EGGGFGF<GEG;G?C3EGFG<GEFG9GGD>-D?7=7DGGF5<@EAFFGC', 'read1_length': 125, 'data': array([[ 0.00000000e+00,  1.16696000e+05,  0.00000000e+00,
         1.25000000e+02,  1.41449997e+02,  0.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  5.66960000e+04,  0.00000000e+00,
         1.25000000e+02,  1.23000000e+02,  1.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 2.00000000e+00,  1.79656000e+05,  0.00000000e+00,
         1.25000000e+02,  1.13000000e+02,  2.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.16846000e+05,  1.25000000e+02,
         2.50000000e+02,  1.37999997e+02,  3.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.17161000e+05,  1.25000000e+02,
         2.50000000e+02,  1.29949997e+02,  5.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.17035000e+05,  1.25000000e+02,
         2.50000000e+02,  1.20749997e+02,  7.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  5.68460000e+04,  1.25000000e+02,
         2.50000000e+02,  1.20000000e+02,  4.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  5.71610000e+04,  1.25000000e+02,
         2.50000000e+02,  1.13000000e+02,  6.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  5.70350000e+04,  1.25000000e+02,
         2.50000000e+02,  1.05000000e+02,  8.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.17106000e+05,  1.33000000e+02,
         2.50000000e+02,  1.11549998e+02,  1.00000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  5.71060000e+04,  1.33000000e+02,
         2.50000000e+02,  9.70000000e+01,  9.00000000e+00,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00]]), 'name': 'simulated_genome_reads.2-chr17-3', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 'ACAGTCCACGTAGTAGATATTATTATATTTATTAGTATTAGAGGCCATTGTGTCAATATATTAGTGTTAGGGCATGGTATGGGTGTCTAGATTAGTGTTAGGGAATGGTGTGGATATTGTATTGG', 'last_seen_chrom': 'chr17', 'ri': {0.0: 0, 1.0: 1, 2.0: 2, 3.0: 3, 4.0: 6, 5.0: 4, 6.0: 7, 7.0: 5, 8.0: 8, 9.0: 10, 10.0: 9}, 'score_mat': {'splitter': [0, 0], 'dis_to_next_path': 11.086029052734375, 'chr17-116846--1-2': [True, 129.9499969482422], 'dis_to_normal': 279.4217224121094, 'path_score': 279.4217224121094, 'normal_pairings': 1, 'chr17-116696-1-1': [False, 123.0]}, 'pairing_params': (100.0, 17.0, 100.0, 1.0, 3.0, 2.0, 9.0)}

    sam = fixsam(template)
    print type(sam)
    print len(sam)
    for l in sam:
        print l

