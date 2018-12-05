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

    template = {'splitters': [True, True], 'paired_end': 1, 'locs': ['chr10-128464470--1-1', 'chr2-32916476-1-1', 'chr10-128464365-1-2', 'chr10-128464505--1-2'], 'bias': 1.0, 'fq_read2_seq': 0, 'isize': (545.0, 160.0), 'read2_q': u"??<@?@@?A???@<===><>=>><><<><=>=>>>=?===<=<<>1<>>><@=9<<=;=+>/=9/;90?;>?>;?'/4;'2;=?-9&3((*'8(+'337:&16/6=@=+&666==A=A:&9(86''&&.&.:0&55005>?B+,6;6>", 'max_d': 1185.0, 'read2_seq': u'TCCTCCCCATCCTCCCCTCTAGACTTCTCCTCCTCCTCCACCATCTTCCCCTCTAGACTTCTCCTCCTCCTCCCCCTCCCCCCCCCTAAACCTTCCCCCCTCCTCCCCCCCCCCCCCCCTAAACCCCCCCCCCCCCCCCCCCCCCCCC', 'chrom_ids': {u'chr10': 0, u'chr2': 1}, 'read2_length': 148, 'passed': True, 'replace_hard': 0, 'read2_reverse': 0, 'inputdata': [[u'69', u'chr10', u'128464365', u'0', u'*', u'=', u'128464365', u'0', u'CCCCCCTCCTCCCCCCCCCCCCCCCTAAACCCCCCCCCCCCCCCCCCCCCCCCC', u"-*3358&+7.0<?<.&554<<;;>8&8(45'&&&)&)80&54114:;?114;5;", u'AS:i:0', u'XS:i:0', u'RG:Z:0'], ['137', u'chr10', u'128464365', u'0', u'6M1D69M73S', u'chr2', u'32916476', u'0', u'TCCTCCCCATCCTCCCCTCTAGACTTCTCCTCCTCCTCCACCATCTTCCCCTCTAGACTTCTCCTCCTCCTCCCCCTCCCCCCCCCTAAACCTTCCCCCCTCCTCCCCCCCCCCCCCCCTAAACCCCCCCCCCCCCCCCCCCCCCCCC', u"??<@?@@?A???@<===><>=>><><<><=>=>>>=?===<=<<>1<>>><@=9<<=;=+>/=9/;90?;>?>;?'/4;'2;=?-9&3((*'8(+'337:&16/6=@=+&666==A=A:&9(86''&&.&.:0&55005>?B+,6;6>", u'NM:i:2', u'MD:Z:6^A39C29', u'AS:i:64', u'XS:i:64', u'RG:Z:0'], [u'393', u'chr10', u'128464510', u'0', u'6M1D69M73H', u'=', u'128464510', u'0', u'*', u'*', u'NM:i:2', u'MD:Z:6^A39C29', u'AS:i:64', u'RG:Z:0'], [u'393', u'chr10', u'128464476', u'0', u'6M1D67M75H', u'=', u'128464476', u'0', u'*', u'*', u'NM:i:2', u'MD:Z:6^A39C27', u'AS:i:62', u'RG:Z:0'], [u'393', u'chr10', u'128464544', u'0', u'6M1D33M1I35M73H', u'=', u'128464544', u'0', u'*', u'*', u'NM:i:3', u'MD:Z:6^A38C29', u'AS:i:56', u'RG:Z:0'], [u'393', u'chr10', u'128464399', u'0', u'6M1D33M1I34M74H', u'=', u'128464399', u'0', u'*', u'*', u'NM:i:4', u'MD:Z:6^A38C1T26', u'AS:i:50', u'RG:Z:0'], [u'393', u'chr10', u'128464362', u'0', u'30H43M75H', u'=', u'128464362', u'0', u'*', u'*', u'NM:i:1', u'MD:Z:15C27', u'AS:i:38', u'RG:Z:0'], ['65', u'chr2', u'32916476', u'0', u'59S11M3D78M', u'chr10', u'128464365', u'0', u'GATGGTGGAGGAGGAGGAGAAGTCTAGAGGGGAGGATGGGGAGGAGGAGGGGAAGGCTAGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGTCGGGGGGGGGGTGGGGGGGGGGGGGGGGAGGGGCGGGGGGGGGGGGAGGGGGGGGGG', u">><>??@??@???===<==;<===>=>>><<=>>>>=>=>=';?8*:09:%05<<-&,(1%.9:?&3:-&53<<<:>*&4&&-48)&&).7&,&*8>;5692&&-/5&/5&&-&&-/8/.0&/&-///9<?99<<'/&&/4&/5&..&", u'NM:i:8', u'MD:Z:11^GGG20G0G10G16G4G23', u'AS:i:55', u'XS:i:54', u'RG:Z:0', u'SA:Z:chr10,128464573,-,98S44M1I5M,0,1;'], [u'353', u'chr2', u'32916463', u'0', u'59H89M', u'chr10', u'128464573', u'0', u'*', u'*', u'NM:i:7', u'MD:Z:23A7G0G10G16G4G12G10', u'AS:i:54', u'RG:Z:0'], [u'353', u'chr2', u'32916391', u'0', u'59H47M4D42M', u'chr10', u'128464573', u'0', u'*', u'*', u'NM:i:9', u'MD:Z:31G0G10G3^GGGA13G4G23', u'AS:i:54', u'RG:Z:0'], [u'353', u'chr2', u'32916431', u'0', u'59H28M5I56M', u'chr10', u'128464573', u'0', u'*', u'*', u'NM:i:9', u'MD:Z:38G3A17G12G10', u'AS:i:53', u'RG:Z:0'], [u'2161', u'chr10', u'128464573', u'0', u'98H44M1I5M', u'=', u'128464573', u'0', u'CCTCCTCCTCCCCATCCTCCCCTCTAGACTTCTCCTCCTCCTCCACCATC', u":90:*8?;'=>=>=>>>>=<<>>>=>===<;==<===???@??@??><>>", u'NM:i:1', u'MD:Z:49', u'AS:i:44', u'XS:i:44', u'RG:Z:0', u'SA:Z:chr2,32916476,+,59S11M3D78M,0,8;'], [u'353', u'chr2', u'32916547', u'0', u'59H24M9I32M1I23M', u'chr10', u'128464573', u'0', u'*', u'*', u'NM:i:13', u'MD:Z:34C5C27G10', u'AS:i:44', u'RG:Z:0'], ['2161', u'chr10', u'128464470', u'0', u'97H12M1D39M', u'chr2', u'32916476', u'101', u'*', u'*', u'NM:i:1', u'MD:Z:12^A39', u'AS:i:44', u'RG:Z:0'], [u'369', u'chr10', u'128464505', u'0', u'98H11M1D39M', u'=', u'128464573', u'67', u'*', u'*', u'NM:i:1', u'MD:Z:11^A39', u'AS:i:43', u'RG:Z:0'], [u'353', u'chr2', u'32916324', u'0', u'59H14M1D17M2I56M', u'chr10', u'128464573', u'0', u'*', u'*', u'NM:i:9', u'MD:Z:14^A8T1A1A14G16G4G23', u'AS:i:42', u'RG:Z:0'], [u'369', u'chr10', u'128464362', u'0', u'100H9M1D39M', u'=', u'128464573', u'212', u'*', u'*', u'NM:i:1', u'MD:Z:9^A39', u'AS:i:41', u'RG:Z:0'], [u'369', u'chr10', u'128464539', u'0', u'98H11M1D33M1I5M', u'=', u'128464573', u'34', u'*', u'*', u'NM:i:2', u'MD:Z:11^A38', u'AS:i:37', u'RG:Z:0'], [u'369', u'chr10', u'128464394', u'0', u'98H11M1D33M1I5M', u'=', u'128464573', u'179', u'*', u'*', u'NM:i:2', u'MD:Z:11^A38', u'AS:i:37', u'RG:Z:0'], [u'145', u'chr10', u'128464573', u'0', u'44M1I5M', u'chr2', u'32916476', u'0', u'CCTCCTCCTCCCCATCCTCCCCTCTAGACTTCTCCTCCTCCTCCACCATC', u":90:*8?;'=>=>=>>>>=<<>>>=>===<;==<===???@??@??><>>", u'NM:i:1', u'MD:Z:49', u'AS:i:44', u'XS:i:43', u'RG:Z:0'], ['2193', u'chr10', u'128464505', u'0', u'11M1D39M', u'chr10', u'128464365', u'0', u'*', u'*', u'NM:i:1', u'MD:Z:11^A39', u'AS:i:43', u'RG:Z:0'], [u'401', u'chr10', u'128464471', u'0', u'11M1D39M', u'chr2', u'32916476', u'0', u'*', u'*', u'NM:i:1', u'MD:Z:11^A39', u'AS:i:43', u'RG:Z:0'], [u'401', u'chr10', u'128464360', u'0', u'11M1D39M', u'chr2', u'32916476', u'0', u'*', u'*', u'NM:i:2', u'MD:Z:1T9^A39', u'AS:i:41', u'RG:Z:0'], [u'401', u'chr10', u'128464539', u'0', u'11M1D33M1I5M', u'chr2', u'32916476', u'0', u'*', u'*', u'NM:i:2', u'MD:Z:11^A38', u'AS:i:37', u'RG:Z:0'], [u'401', u'chr10', u'128464394', u'0', u'11M1D33M1I5M', u'chr2', u'32916476', u'0', u'*', u'*', u'NM:i:2', u'MD:Z:11^A38', u'AS:i:37', u'RG:Z:0']], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'rows': [1, 20, 13, 7], 'read1_q': u">><>??@??@???===<==;<===>=>>><<=>>>>=>=>=';?8*:09:%05<<-&,(1%.9:?&3:-&53<<<:>*&4&&-48)&&).7&,&*8>;5692&&-/5&/5&&-&&-/8/.0&/&-///9<?99<<'/&&/4&/5&..&", 'read1_length': 148, 'data': array([[ 0.00000000e+00,  1.28464573e+08,  0.00000000e+00,
         5.00000000e+01,  4.40000000e+01,  1.10000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464470e+08,  0.00000000e+00,
         5.10000000e+01,  4.40000000e+01,  1.30000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464505e+08,  0.00000000e+00,
         5.00000000e+01,  4.30000000e+01,  1.40000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464362e+08,  0.00000000e+00,
         4.80000000e+01,  4.10000000e+01,  1.60000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464539e+08,  0.00000000e+00,
         5.00000000e+01,  3.70000000e+01,  1.70000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464394e+08,  0.00000000e+00,
         5.00000000e+01,  3.70000000e+01,  1.80000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464365e+08,  0.00000000e+00,
         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  3.29164760e+07,  5.90000000e+01,
         1.48000000e+02,  5.50000000e+01,  7.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  3.29164630e+07,  5.90000000e+01,
         1.48000000e+02,  5.40000000e+01,  8.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  3.29163910e+07,  5.90000000e+01,
         1.48000000e+02,  5.40000000e+01,  9.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  3.29164310e+07,  5.90000000e+01,
         1.48000000e+02,  5.30000000e+01,  1.00000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  3.29165470e+07,  5.90000000e+01,
         1.48000000e+02,  4.40000000e+01,  1.20000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 1.00000000e+00,  3.29163240e+07,  5.90000000e+01,
         1.48000000e+02,  4.20000000e+01,  1.50000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464365e+08,  1.48000000e+02,
         2.23000000e+02,  6.40000000e+01,  1.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464510e+08,  1.48000000e+02,
         2.23000000e+02,  6.40000000e+01,  2.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464476e+08,  1.48000000e+02,
         2.21000000e+02,  6.20000000e+01,  3.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464544e+08,  1.48000000e+02,
         2.23000000e+02,  5.60000000e+01,  4.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464399e+08,  1.48000000e+02,
         2.22000000e+02,  5.00000000e+01,  5.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464573e+08,  1.48000000e+02,
         1.98000000e+02,  4.40000000e+01,  1.90000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464362e+08,  1.78000000e+02,
         2.21000000e+02,  3.80000000e+01,  6.00000000e+00,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464505e+08,  2.46000000e+02,
         2.96000000e+02,  4.30000000e+01,  2.00000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464471e+08,  2.46000000e+02,
         2.96000000e+02,  4.30000000e+01,  2.10000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464360e+08,  2.46000000e+02,
         2.96000000e+02,  4.10000000e+01,  2.20000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464539e+08,  2.46000000e+02,
         2.96000000e+02,  3.70000000e+01,  2.30000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00],
       [ 0.00000000e+00,  1.28464394e+08,  2.46000000e+02,
         2.96000000e+02,  3.70000000e+01,  2.40000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00]]), 'name': u'HISEQ1:9:H8962ADXX:1:1101:4291:82991', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': u'GATGGTGGAGGAGGAGGAGAAGTCTAGAGGGGAGGATGGGGAGGAGGAGGGGAAGGCTAGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGTCGGGGGGGGGGTGGGGGGGGGGGGGGGGAGGGGCGGGGGGGGGGGGAGGGGGGGGGG', 'last_seen_chrom': u'chr10', 'ri': {0.0: 6, 1.0: 13, 2.0: 14, 3.0: 15, 4.0: 16, 5.0: 17, 6.0: 19, 7.0: 7, 8.0: 8, 9.0: 9, 10.0: 10, 11.0: 0, 12.0: 11, 13.0: 1, 14.0: 2, 15.0: 12, 16.0: 3, 17.0: 4, 18.0: 5, 19.0: 18, 20.0: 20, 21.0: 21, 22.0: 22, 23.0: 23, 24.0: 24}, 'score_mat': {'splitter': [1, 1], 'dis_to_next_path': 0.0, 'chr2-32916476-1-1': [True, 43.0], 'chr10-128464470--1-1': [False, 64.0], 'dis_to_normal': 12.0, 'path_score': 146.0, 'chr10-128464365-1-2': [True, 44.0], 'normal_pairings': 0, 'chr10-128464505--1-2': [False, 54.0]}, 'pairing_params': (100.0, 17.0, 100.0, 1.0, 3.0, 2.0, 9.0)}

    sam = fixsam(template)

    for l in sam:
        print l

