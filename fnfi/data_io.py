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
            "paired_end": int(args["paired"] == "True"),
            "inputdata": rows,
            "bias": args["bias"],
            "read1_length": 0,
            "read2_length": 0,
            "score_mat": {},
            "passed": 0,
            "name": rows[0][0][0],
            "last_seen_chrom": last_seen_chrom,
            "inputfq": fq,
            "read1_seq": 0,  # Some sam records may have seq == '*' , need a record of full seq for adding back in
            "read2_seq": 0,
            "read1_q": 0,
            "read2_q": 0,
            "read1_reverse": 0,  # Set to true if aligner has reverse complemented the sequence
            "read2_reverse": 0,
            "replace_hard": int(args["replace_hardclips"] == "True"),
            "fq_read1_seq": 0,
            "fq_read2_seq": 0,
            "fq_read1_q": 0,
            "fq_read2_q": 0
            }


def to_output(template):

    # Todo make sure all read-pairs have a mapping, otherwise write an unmapped
    #
    # paired = False if template["read2_length"] is None else True
    sam = samclips.fixsam(template)

    if len(sam) == 0:  # Todo fix unmapped reads
        # print("No alignments")
        # click.echo(template, err=True)
        # click.echo("Unmapped read error", err=True)
        # quit()
        return None

    # if len(sam) == 1:  # Todo deal with these
    #     click.echo("Single alignment only", err=True)
    #    return None

    if any(i[0] == "*" or i[4] == "*" for i in sam):  # Unformatted cigar or unformatted cigarstring
        # click.echo("Unformatted cigar error", err=True)
        # click.echo(sam, err=True)
        # quit()
        return None

    return "".join(template["name"] + "\t" + "\t".join(i) + "\n" for i in sam)


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


def intersecter(tree, chrom, start, end):
    if tree is None:
        return 0
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return 1
        else:
            return 0
    else:
        return 0


def sam_itr(args):

    itr = args["sam"]
    tree = overlap_regions(args["include"])

    # First get header
    header_string = ""
    last_seen_chrom = ""
    first_line = ""
    for t in itr:

        t = str(t.decode("ascii"))

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
        t = str(t.decode("ascii"))
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
    # else:
    #     click.echo("WARNING: fastq not found", err=True)
    #     return None, None


def iterate_mappings(args):

    inputstream = sam_itr(args)

    total = 0
    name = ""
    rows = []
    header_string = next(inputstream)
    #header_string += "@RG\tID:0\tSM:0\tPU:lane1\tPL:ILLUMINA\tLB:0\n"  # Todo deal with read group information
    # Todo add somthing to RG header, indicating fnfi alignment

    yield header_string

    fq_buffer = defaultdict(list)
    fq_iter = fq_reader(args)

    max_d = args["insert_median"] + 4*args["insert_stdev"]  # Separation distance threshold to call a pair discordant
    last_seen_chrom = ""
    for m, last_seen_chrom, ol in inputstream:  # Alignment

        nm = m[0]
        if name != nm:
            if len(rows) > 0:
                total += 1

                fq = fq_getter(fq_iter, name, args, fq_buffer)
                yield (rows, args, max_d, last_seen_chrom, fq)
                # if name == "chr22-163734":
                #     break
            rows = []
            name = nm

        rows.append((m, ol))  # String, ol states if alignment overlaps ROI

    # Deal with last record
    if len(rows) > 0:
        total += 1
        fq = fq_getter(fq_iter, name, args, fq_buffer)
        yield (rows, args, max_d, last_seen_chrom, fq)

    click.echo("Total processed " + str(total) + "\n", err=True)


if __name__ == "__main__":
    import time
    import c_io_funcs
    import c_samflags
    import pairing

    sam = [[u't', u'97', u'chr22', u'50629342', u'0', u'100M', u'=', u'50630825', u'1539', u'CCTACTTGCCTGATCTGGGAAGAGTAACAGTCTGACGCCTTTAAAGGCCTGAAGGGAACATTCACCATCTGTTCTCTCTCAGGGCTGCTACCTGTGAGGT', u'AADBDDCDAAACBACBA@?ABA=>C?>@@B?B<BB?@A>B??BAB@?AC?@?CBC@>=@@@>>??=A>=@B@><@????B?B<>>;<<><=@:;=??<:=', u'NM:i:1', u'MD:Z:3C96', u'AS:i:96', u'XS:i:96'],
           [u't', u'353', u'chr22', u'50639946', u'0', u'100M', u'=', u'50630825', u'-9067', u'*', u'*', u'NM:i:1', u'MD:Z:3C96', u'AS:i:96'],
           [u't', u'145', u'chr22', u'50630825', u'60', u'44S56M', u'=', u'50629342', u'-1539', u'CTCTCCTGGGAGCCTCCGTCAGGGGAACCCAGGAGCCGACGTTCGTCTTAAGAAGTCCCTGCCAGGCGAGCTGTCAGAGCCCCGGCACTGGGAGTGGTGG', u'==;:=@=<<><=<>><>=A>>><>=<==@@==<=>A??A?>=@>?@>@A=A?A???>@AA?@AB?@??@A@A@C@@>BC@BCB@@CBA@AAB@BB?A?AD', u'NM:i:0', u'MD:Z:56', u'AS:i:56', u'XS:i:0', u'SA:Z:chr22,50629781,-,40M60S,60,0;'],
           [u't', u'2193', u'chr22', u'50629781', u'60', u'40M60H', u'=', u'50629342', u'-479', u'CTCTCCTGGGAGCCTCCGTCAGGGGAACCCAGGAGCCGAC', u'==;:=@=<<><=<>><>=A>>><>=<==@@==<=>A??A?', u'NM:i:0', u'MD:Z:40', u'AS:i:40', u'XS:i:0', u'SA:Z:chr22,50630825,-,44S56M,60,0;']]

    sam = [

[u't', u'97', u'chrUn_KI270442v1', u'25262', u'60', u'76M10D72M', u'chr4', u'49141816', u'0', u'GGAATGTAGTGGAGTGTATTGGAATGGAAGGGAATGGAATGGAATGGAATTGAGTGAAGTTGAGTGGAGTGGAATGGAATGGAATGGAAAGGAATGCAACAGAATGACGTGGAGTGGAATGGAATGGAACGGAACAGAATGGAAGTGA', u'??>?@@=@A?<??==>>>>0>>=>>><??:><?=>?<=>==?>=>?<>>>????<??=>==?===>?;>=>>>>>???@@??>????>>?=;@??<?@><@?@@@=@=75?>@?;>0=??@@@?=?2?A?90AA@B>AB@@?BA><=@', u'NM:i:18', u'MD:Z:76^GAGTGTAAGA18C0A0G4A1G4A20T14T3', u'AS:i:93', u'XS:i:41', u'RG:Z:0'],
[u't', u'145', u'chr4', u'49141816', u'0', u'30M', u'chrUn_KI270442v1', u'25262', u'0', u'TTCCATTCTGTTCCGTTCCATTCCATTCCA', u'==;=<==;?=>?/7;==0;9=>7?<></>6', u'NM:i:0', u'MD:Z:30', u'AS:i:30', u'XS:i:30', u'RG:Z:0'],
[u't', u'401', u'chr4', u'49139377', u'0', u'30M', u'chrUn_KI270442v1', u'25262', u'0', u'*', u'*', u'NM:i:0', u'MD:Z:30', u'AS:i:30', u'RG:Z:0'],
[u't', u'385', u'chrUn_KI270442v1', u'1073', u'0', u'30M', u'=', u'25262', u'24190', u'*', u'*', u'NM:i:0', u'MD:Z:30', u'AS:i:30', u'RG:Z:0'],
[u't', u'401', u'chr22', u'10724538', u'0', u'30M', u'chrUn_KI270442v1', u'25262', u'0', u'*', u'*', u'NM:i:0', u'MD:Z:30', u'AS:i:30', u'RG:Z:0'],
[u't', u'99', u'chrUn_KI270442v1', u'25262', u'60', u'50M', u'=', u'25623', u'509', u'GGAATGTAGTGGAGTGTATTGGAATGGAAGGGAATGGAATGGAATGGAAT', u'??>?@@=@A?<??==>>>>0>>=>>><??:><?=>?<=>==?>=>?<>>>', u'NM:i:0', u'MD:Z:50', u'AS:i:50', u'XS:i:33', u'RG:Z:0'],
[u't', u'147', u'chrUn_KI270442v1', u'25623', u'60', u'148M', u'=', u'25262', u'-509', u'AAATCGAATGTCATCAAATGGAATAAAATGGAATATGGACAAGTGTAGAAGAGTGAATTGGAATACAGTGGAAGGGAATTGGGTGGATTGGAATTTAATGTACTGGAGTGGAGTGGAACGGAGTGGATGGAATGGAATGGGGAGAAAT', u'7><>6AB?>=)????BB?<BA>@<ABBAAA@A@@>A@?><@@=?>@A@>?=@=????=A??@A?>?@>@>?>?==<=<??>>=?>>>??>>>=??=>>?=>==?=>>=>==?<><===4><?<>==>>=<==?>>>?@?=?>@>><>>', u'NM:i:7', u'MD:Z:24C13T0G33C8A19G15A29', u'AS:i:113', u'XS:i:26', u'RG:Z:0']
]
    sam_ol = [["\t".join(i).split("\t", 4), False] for i in sam]  # Identifier for overlapping region of interest

    args = {"max_insertion": 100, "min_aln": 17, "max_overlap": 100, "ins_cost": 1, "ol_cost": 3, "inter_cost":2,
            "u": 9, "match_score": 1, "insert_median": 500, "insert_stdev": 50, "bias": 1.15, "paired": "True", "replace_hardclips": False, "fq1": None, "fq2": None}
    template = make_template(sam_ol, args, 2000, "chr22", (0, 0))

    # t0 = time.time()
    # sam_to_array(template)
    # v0 = time.time() - t0
    # print template["data"].astype(int)
    template = {'inputfq': (None, None), 'paired_end': 1, 'bias': 1.0, 'fq_read2_seq': 0, 'isize': (545.0, 160.0), 'read2_q': 0, 'max_d': 1185.0, 'read2_seq': 0, 'read2_length': 0, 'passed': 0, 'replace_hard': 0, 'read2_reverse': 0, 'inputdata': [([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'81', u'chr10', u'41899773', u'17\t52S96M\tchr20\t31185060\t0\tGAATGCACAAGTGCAATGGAATGGAATGGAATGAAATGGTATCAACCCGGTTGGAATGGAATGGAATGGAATGGAATGGAATGGATTCAACCCGAGTGCAATGGAATGGATTGCAATGGAATGGAATCAACCCGATTGGAATGGAATG\t@=<?<?@@BB?@@8BBBAAAAA?A@AA>?AA@A@@@?@>=@@@>?A@7=<4A?@>@@??@@?????@=??==>=>>?>>?>?==>>=?><<;4<@==?>>>?=<?<?>>==>>?<<>=<==>=====>==??6>>>@????@@??<@?\tNM:i:7\tMD:Z:46G11A2G7C5A3A3A12\tAS:i:61\tXS:i:52\tRG:Z:0\tSA:Z:chr20,31185784,-,14S71M63S,0,4;\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr4', u'49643176', u'0\t51H58M6D4M1I34M\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:16A4T12A9T2G10^AATAGA38\tAS:i:52\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'2129', u'chr20', u'31185784', u'0\t14H71M63H\t=\t31185060\t-795\tAATGGAATGGAATGGAATGAAATGGTATCAACCCGGTTGGAATGGAATGGAATGGAATGGAATGGAATGGA\tBBBAAAAA?A@AA>?AA@A@@@?@>=@@@>?A@7=<4A?@>@@??@@?????@=??==>=>>?>>?>?==>\tNM:i:4\tMD:Z:19C5A9A0A34\tAS:i:51\tXS:i:51\tRG:Z:0\tSA:Z:chr10,41899773,-,52S96M,17,7;\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr5', u'49659969', u'0\t31M20D66M51H\tchr20\t31185060\t0\t*\t*\tNM:i:25\tMD:Z:22G8^ATTCCACTCCATTCCACACA3C2T11C2A44\tAS:i:49\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'11449', u'0\t21H77M50H\tchr20\t31185060\t0\t*\t*\tNM:i:6\tMD:Z:12A0C2T3T7C12T35\tAS:i:47\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'29157', u'0\t51M10I36M51H\tchr20\t31185060\t0\t*\t*\tNM:i:15\tMD:Z:34C2T3G2T4C37\tAS:i:46\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31070378', u'0\t14H89M5D45M\t=\t31185060\t114545\t*\t*\tNM:i:21\tMD:Z:26G4A3A0C16A17A1A0C3A5G4^GAATA7A2G6G11A0C1G12\tAS:i:46\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31052141', u'0\t50H32M5D44M22H\t=\t31185060\t132840\t*\t*\tNM:i:9\tMD:Z:32^GGAGA3G0G3A7G27\tAS:i:45\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17_KI270729v1_random', u'27570', u'0\t22H63M63H\tchr20\t31185060\t0\t*\t*\tNM:i:4\tMD:Z:11G5A5A3A35\tAS:i:43\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr4', u'49642058', u'0\t50H72M26H\tchr20\t31185060\t0\t*\t*\tNM:i:6\tMD:Z:35A1A10G4C6A2A8\tAS:i:42\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41880880', u'0\t50H46M15D52M\tchr20\t31185060\t0\t*\t*\tNM:i:22\tMD:Z:35A1A8^TGGAATGTAAAGGAA2G4C6A21A2G12\tAS:i:42\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrUn_KI270442v1', u'96019', u'0\t14H71M63H\tchr20\t31185060\t0\t*\t*\tNM:i:6\tMD:Z:9C3A6T4A6T2A35\tAS:i:41\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10997113', u'0\t63H33M10I34M8H\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:35T5C25\tAS:i:41\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'19389', u'0\t51M10I36M51H\tchr20\t31185060\t0\t*\t*\tNM:i:16\tMD:Z:34C2T3G2T4C16G20\tAS:i:41\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'11690416', u'0\t11H31M1D3M6I33M15I28M21H\tchr20\t31185060\t0\t*\t*\tNM:i:25\tMD:Z:28A2^G20C26A16\tAS:i:40\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr21', u'10670478', u'0\t21H27M5I44M51H\tchr20\t31185060\t0\t*\t*\tNM:i:9\tMD:Z:29C0T0T4G34\tAS:i:40\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31066660', u'0\t50H35M6D3M1I20M5I34M\t=\t31185060\t118304\t*\t*\tNM:i:17\tMD:Z:35^ATGGAC5A9C19C8G9G2\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56753276', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56703386', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17', u'21987214', u'0\t53H74M21H\tchr20\t31185060\t0\t*\t*\tNM:i:7\tMD:Z:9T8C0A9T2A5G6G28\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56721172', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49093874', u'0\t33M5I10M15I78M7H\tchr20\t31185060\t0\t*\t*\tNM:i:30\tMD:Z:12T23T12T11T15C0T9T1T0G2G26\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr10', u'42321421', u'0\t58H39M51H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:39\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56760389', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17', u'21988030', u'0\t53H74M21H\tchr20\t31185060\t0\t*\t*\tNM:i:7\tMD:Z:9T8C0A9T2A5G6G28\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56689090', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'17792', u'0\t26H64M58H\tchr20\t31185060\t0\t*\t*\tNM:i:5\tMD:Z:8A2T14T2C0A33\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56739054', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56685532', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56699807', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49095483', u'0\t33M5I10M15I78M7H\tchr20\t31185060\t0\t*\t*\tNM:i:30\tMD:Z:12T23T12T11T15C0T9T1T0G2G26\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56749802', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56710469', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56724760', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56735510', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31055501', u'0\t51H39M58H\t=\t31185060\t129522\t*\t*\tNM:i:0\tMD:Z:39\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr2', u'89820762', u'0\t58H39M51H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:39\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56714053', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56742633', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41881276', u'0\t11H27M10I36M10D43M21H\tchr20\t31185060\t0\t*\t*\tNM:i:27\tMD:Z:22G5A34^AGTGCAATGA1A9T2G11A1A14\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56696228', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56763963', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56728354', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:11\tMD:Z:24G41\tAS:i:39\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49105323', u'0\t33M5I10M15I78M7H\tchr20\t31185060\t0\t*\t*\tNM:i:31\tMD:Z:4A7T23T12T11T4C10C0T9T2G2G26\tAS:i:38\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49145152', u'0\t63H37M48H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:37\tAS:i:37\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10673906', u'0\t63H34M4D3M9I30M9H\tchr20\t31185060\t0\t*\t*\tNM:i:14\tMD:Z:34^TTGC8C24\tAS:i:37\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr17_KI270729v1_random', u'3476', u'0\t21H75M10I33M9H\tchr20\t31185060\t0\t*\t*\tNM:i:21\tMD:Z:13C2T11C2T2C0A5T13T3G17T5G24\tAS:i:37\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49144832', u'0\t63H37M48H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:37\tAS:i:37\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'16124', u'0\t63H35M1D2M6I28M14H\tchr20\t31185060\t0\t*\t*\tNM:i:9\tMD:Z:35^T4T5C19\tAS:i:36\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41915043', u'0\t14H71M63H\tchr20\t31185060\t0\t*\t*\tNM:i:7\tMD:Z:17C1G5A2G0T4A0A35\tAS:i:36\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17', u'21981677', u'0\t51H97M\tchr20\t31185060\t0\t*\t*\tNM:i:13\tMD:Z:11T8C0A9T2A5G0T5G28G0G2T0G2A12\tAS:i:36\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41911833', u'0\t14H71M63H\tchr20\t31185060\t0\t*\t*\tNM:i:7\tMD:Z:19G0C3T0A2A3A2A35\tAS:i:36\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr4', u'49643130', u'0\t59H29M1I20M16D4M1I34M\tchr20\t31185060\t0\t*\t*\tNM:i:21\tMD:Z:26A3T7G10^AATGGAATGGAATAGA38\tAS:i:36\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49108255', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr1', u'224012042', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrUn_KI270442v1', u'90992', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_KI270756v1', u'15521', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr10', u'38805127', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'49523', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'11706730', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17', u'21901259', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49102375', u'0\t27M10D24M10I37M50H\tchr20\t31185060\t0\t*\t*\tNM:i:25\tMD:Z:12C14^GTCCGTTCCT7T2T11C2T35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41898002', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr22_KI270736v1_random', u'147674', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10925512', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10869945', u'0\t63H33M10I33M9H\tchr20\t31185060\t0\t*\t*\tNM:i:13\tMD:Z:20C14T5C24\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'11022129', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10770511', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr18', u'6688464', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17', u'21910071', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'11298185', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10794699', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41880407', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'41891562', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10887351', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10867988', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31075665', u'0\t50H35M63H\t=\t31185060\t109362\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr17_KI270729v1_random', u'1145', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31092116', u'0\t50H35M63H\t=\t31185060\t92911\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10814486', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10997722', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31209581', u'0\t50H35M63H\t=\t31185060\t-24556\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'38864408', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr22_KI270737v1_random', u'32024', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr16', u'34063105', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr17_KI270729v1_random', u'2783', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10809126', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31091422', u'0\t50H35M63H\t=\t31185060\t93605\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49120493', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'10935595', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'11643726', u'0\t50H37M5I35M5I16M\tchr20\t31185060\t0\t*\t*\tNM:i:17\tMD:Z:35A1G2A2G11A2G13A15\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr20', u'31235175', u'0\t50H35M63H\t=\t31185060\t-50150\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr10', u'38502809', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrY', u'11303458', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr21', u'7922716', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr22', u'16353746', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr2', u'89838691', u'0\t50H35M63H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chrUn_GL000216v2', u'6379', u'0\t63H35M50H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:35\tAS:i:35\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56771100', u'0\t9H29M4I3M6I34M63H\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:24G29T11\tAS:i:34\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr4', u'49150794', u'0\t50H44M54H\tchr20\t31185060\t0\t*\t*\tNM:i:2\tMD:Z:32A1G9\tAS:i:34\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56717627', u'0\t9H29M4I3M6I33M64H\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:24G30G9\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56692669', u'0\t9H29M4I3M6I33M64H\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:24G30G9\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr17', u'21960369', u'0\t58H31M5I33M21H\tchr20\t31185060\t0\t*\t*\tNM:i:9\tMD:Z:14C12A4A2G28\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56676395', u'0\t9H29M4I3M6I33M64H\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:24G30G9\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56706955', u'0\t9H29M4I3M6I33M64H\tchr20\t31185060\t0\t*\t*\tNM:i:12\tMD:Z:24G30G9\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr5', u'49601887', u'0\t7H35M6I2M1D45M53H\tchr20\t31185060\t0\t*\t*\tNM:i:13\tMD:Z:26G5A4^A0A2C7C23A9\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr2', u'89839075', u'0\t114H34M\tchr20\t31185060\t0\t*\t*\tNM:i:1\tMD:Z:33A0\tAS:i:33\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chr4', u'49640470', u'0\t62H22M4I5M1I32M22H\tchr20\t31185060\t0\t*\t*\tNM:i:7\tMD:Z:28A2G27\tAS:i:32\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr10', u'42316894', u'0\t68H32M48H\tchr20\t31185060\t0\t*\t*\tNM:i:0\tMD:Z:32\tAS:i:32\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr5', u'49661039', u'0\t66H28M15I30M9H\tchr20\t31185060\t0\t*\t*\tNM:i:16\tMD:Z:33C24\tAS:i:32\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'321', u'chr5', u'49658463', u'0\t22M15D60M66H\tchr20\t31185060\t0\t*\t*\tNM:i:21\tMD:Z:22^TTCCGTTCTATTCTT12C2T3C7C12T8T10\tAS:i:31\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'337', u'chrY', u'56756840', u'0\t9H32M10I35M62H\tchr20\t31185060\t0\t*\t*\tNM:i:14\tMD:Z:24G5A14G9T11\tAS:i:31\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'161', u'chr20', u'31185060', u'44\t148M\tchr10\t41899773\t0\tAGTGGAATGGAATGGAATGGAAAGGAATGGAGTGGAATTCAATAGAATCAACCCTAGTGGAATGGAATGGAAAGGAATTGAATGGAATGGAAAGGAATGGAATGGAATTCAATCCAGTGTAAAGGAATGGAATGGAATGGAATGAAAT\t??>?????@??????>>=><;<>=>;<;>>=>=>><><>==>==?<=<>??<==>>?==<?=>??>=??>>?>>??=>=>?@????@??@@@??A@@??=?@?@@@@@@=?>>=4?>>=3=@?@@A<>?<A??A?A??@>BB?A>=>+\tNM:i:9\tMD:Z:19T30T3G20T6C9T15C23C11G3\tAS:i:104\tXS:i:78\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'417', u'chr20', u'31053250', u'0\t108M40H\tchr10\t41899773\t0\t*\t*\tNM:i:6\tMD:Z:19T4A15T11A1G37T15\tAS:i:78\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'417', u'chr20', u'31057609', u'0\t108M40H\tchr10\t41899773\t0\t*\t*\tNM:i:6\tMD:Z:19T4A15T11A1G37T15\tAS:i:78\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'417', u'chr20', u'31055049', u'0\t108M40H\tchr10\t41899773\t0\t*\t*\tNM:i:8\tMD:Z:1A29A20G1G7G6C8G13T15\tAS:i:71\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'433', u'chr2', u'89828605', u'0\t40H108M\tchr10\t41899773\t0\t*\t*\tNM:i:10\tMD:Z:15A13C5A20T7C3C0C6T8A20T1\tAS:i:61\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'417', u'chr10', u'41883362', u'0\t110M5D38M\t=\t41899773\t16507\t*\t*\tNM:i:21\tMD:Z:22T1C6A6G0G3G28T5G13T4C10G1^AATGG4A1C2G2T21G3\tAS:i:58\tRG:Z:0\n'], None), ([u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', u'417', u'chr20', u'31070975', u'0\t57H91M\tchr10\t41899773\t0\t*\t*\tNM:i:8\tMD:Z:15T2T6G9T10A4C23C11G3\tAS:i:52\tRG:Z:0\n'], None)], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'read1_q': 0, 'read1_length': 0, 'name': u'HISEQ1:9:H8962ADXX:1:1101:1144:90311', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 0, 'last_seen_chrom': u'chr17', 'score_mat': {}, 'pairing_params': (100.0, 17.0, 100.0, 1.0, 3.0, 2.0, 9.0)}


    # template2 = make_template(sam_ol, args, 2000, "chr22", (0, 0))

    t0 = time.time()
    c_io_funcs.sam_to_array(template)
    v = time.time() - t0

    # print v0, v, v0 / v
    # print [v == v2 for v, v2 in zip(template.values(), template2.values())]

    # res = pairing.process(template)
    #
    # t0 = time.time()
    # apply_filter(template, *res)
    # v0 = time.time() - t0
    #
    # res2 = pairing.process(template2)
    # t0 = time.time()
    # c_io_funcs.apply_filter(template2, *res2)
    # v1 = time.time() - t0
    #
    # # print v0 / v1
    #
    # t0 = time.time()
    # sam = samclips.fixsam(template)
    # v0 = time.time() - t0
    #
    #
    # t0 = time.time()
    # sam = c_samflags.fixsam(template2)
    # v1 = time.time() - t0

    # print v0 / v1
