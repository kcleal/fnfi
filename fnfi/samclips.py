"""
Utils to generate proper sam output and flag information
"""

import re
import click


def echo(*arg):
    click.echo(arg, err=True)


def rev_comp(s):
    r = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "T", "c": "G", "g": "C", "t": "A", "n": "N"}
    return "".join(r[i] for i in s)[::-1]


def set_bit(v, index, x):
    """Set the index:th bit of v to 1 if x is truthy, else to 0, and return the new value."""
    mask = 1 << index  # Compute mask, an integer with just bit 'index' set.
    v &= ~mask  # Clear the bit indicated by the mask (if x is False)
    if x:
        v |= mask  # If x was True, set the bit indicated by the mask.
    return v


def add_softclips(s, fq, readid):
    fq_seq = fq[readid][0].upper()
    fq_qual = fq[readid][1]
    if s.flag & 16:
        fq_seq = rev_comp(fq_seq)
        fq_qual = fq_qual[::-1]
    s.seq = fq_seq
    s.qual = fq_qual


def set_mate_flag(a, b, r1_len, r2_len, max_d, read1_rev, read2_rev, template):

    if not a or not b:  # No alignment, mate unmapped?
        return False, False

    # Make sure chromosome of mate is properly set not "*"
    chrom_a, mate_a = a[2], a[5]
    chrom_b, mate_b = b[2], b[5]
    if chrom_a != mate_b:
        b[5] = chrom_a
    if chrom_b != mate_a:
        a[5] = chrom_b

    aflag = int(a[0])
    bflag = int(b[0])

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
    aflag = set_bit(aflag, 8, 0)
    bflag = set_bit(bflag, 8, 0)

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

    a[0] = str(aflag)
    b[0] = str(bflag)
    return reverse_A, reverse_B


def set_supp_flags(sup, pri, primary_reversed, template):
    # Set paired
    supflag = int(sup[0])
    priflag = int(pri[0])

    # Set supplementary flag
    if not supflag & 1:
        supflag = set_bit(supflag, 0, 1)

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
    if not primary_reversed and supflag & 16:
        rev_sup = True
    else:
        rev_sup = False

    sup[0] = str(supflag)

    sup[5] = pri[1]
    sup[6] = pri[2]
    return rev_sup


def add_sequence_back(atype, item, reverse_me, template):

    c = re.split(r'(\d+)', item[4])[1:]  # Drop leading empty string
    start = 0

    cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "DH"])
    flag = int(item[0])

    if flag & 64 and template["read1_seq"]:
        name = "read1"
        if template["fq_read1_seq"] is not None:
            end = len(template["fq_read1_seq"])
        else:
            end = len(template["read1_seq"])

    elif flag & 128 and template["read2_seq"]:
        name = "read2"
        if template["fq_read2_seq"] is not None:
            end = len(template["fq_read2_seq"])
        else:
            end = len(template["read2_seq"])
    else:
        return  # read sequence is None or bad flag

    # Try and replace H with S
    if c[1] == "H" or c[-1] == "H":
        # Replace hard with soft-clips
        non_hardclipped_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "D"])
        if non_hardclipped_length == end:
            item[4] = item[4].replace("H", "S")
            cigar_length = non_hardclipped_length
        else:
            # Remove seq
            if c[1] == "H":
                start += int(c[0])
            if c[-1] == "H":
                end -= int(c[-2])

    if item[9] == "*":
        if reverse_me:
            item[8] = rev_comp(template["%s_seq" % name])[start:end]
            item[9] = template["%s_q" % name][::-1][start:end]
        else:
            item[8] = template["%s_seq" % name][start:end]
            item[9] = template["%s_q" % name][start:end]

    elif item[9] in template["fq_%s_q" % name]:
        item[8] = template["fq_%s_seq" % name][start:end]
        item[9] = template["fq_%s_q" % name][start:end]

    elif item[9] in template["fq_%s_q" % name][::-1]:
        item[8] = rev_comp(template["fq_%s_seq" % name])[start:end]
        item[9] = template["fq_%s_q" % name][::-1][start:end]
    else:
        echo("---")
        echo(item[9], flag)
        echo(name)
        echo(template["read1_q"])
        echo(template["read2_q"])
        echo(item)
        echo(rev_comp(template["fq_%s_seq" % name]))
        raise ValueError

    assert len(item[8]) == cigar_length


def fixsam(template):

    sam = [template['inputdata'][i] for i in template['rows']]
    r1l, r2l = template["read1_length"], template["read2_length"]
    max_d = template['max_d']

    # Todo make sure all read-pairs have a mapping, otherwise write an unmapped

    paired = False if template["read2_length"] is None else True
    score_mat = template['score_mat']

    out = []
    primary1 = None
    primary2 = None
    rev_A = False
    rev_B = False
    for l in sam:
        t = score_mat
        strand = "-1" if int(l[0]) & 16 else "1"
        rid = 2 if int(l[0]) & 128 else 1
        key = "{}-{}-{}-{}".format(l[1], str(int(l[2])), strand, rid)

        if len(t[key]) > 2:
            # Prevent bug where two identical alignments possible
            aln_info_0, aln_info_1 = t[key].pop(0), t[key].pop(0)  # Remove first two items from list
        else:
            aln_info_0, aln_info_1 = t[key]

        if not aln_info_0:  # Not primary, set supplementary flag
            l[0] = str(set_bit(int(l[0]), 11, 1))

        if rid == 1:
            if t["splitter"][0]:
                split = "1"
            else:
                split = "0"
        elif rid == 2:
            if t["splitter"][1]:
                split = "1"
            else:
                split = "0"
        # "@RG\tID:0\tSM:0\tPU:0\tPL:0\tLB:0\n"
        xs = int(aln_info_1)

        l += ["SP:Z:" + split,
              "DA:i:" + str(xs),
              "DP:Z:" + str(round(t["dis_to_next_path"], 0)),
              "DN:Z:" + str(round(t["dis_to_normal"], 2)),
              "PS:Z:" + str(round(t["path_score"], 2)),
              "NP:Z:" + str(round(t["normal_pairings"], 1)),
              "RG:Z:0", "SM:Z:0", "PU:Z:lane1", "PL:Z:ILLUMINA", "LB:Z:0",
              ]

        if aln_info_0:
            if rid == 1:
                primary1 = l
            else:
                primary2 = l
        else:
            out.append(['sup', l, False])  # Supplementary

    if primary1 is None or primary2 is None and template["paired_end"]:
        return []  # Todo deal with unmapped read or unpaired

    if paired and template["paired_end"]:

        rev_A, rev_B = set_mate_flag(primary1, primary2, r1l, r2l, max_d,
                                     template["read1_reverse"], template["read2_reverse"], template)

        # Check if supplementary needs reverse complementing
        for i in range(len(out)):
            if int(out[i][1][0]) & 64:  # Flag
                revsup = set_supp_flags(out[i][1], primary1, template["read1_reverse"], template)
            else:
                revsup = set_supp_flags(out[i][1], primary2, template["read2_reverse"], template)
            if revsup:
                out[i][2] = True

        # if template["name"] == "simulated_reads.3.10-id293_A_chr21:46699688_B_chr1:38378863-38631":
        #     echo(rev_A, rev_B)
        #     echo(primary1)
        #     echo(primary2)
            # quit()

    out = [('pri', primary1, rev_A), ('pri', primary2, rev_B)] + out

    # Add read seq info back in if necessary, before reverse complementing. Check for hard clips and clip as necessary
    for a_type, item, reverse_me in out:
        if item:
            if len(item[8]) <= 1 or "H" in item[4]:  # Sequence is set as "*", needs adding back in, or hard-clipped
                add_sequence_back(a_type, item, reverse_me, template)
            elif reverse_me:
                item[8] = rev_comp(item[8])
                item[9] = item[9][::-1]
        else:
            # None here means no alignment for primary2
            pass

    return [i[1] for i in out if i[1] is not None]
