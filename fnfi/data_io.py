from __future__ import absolute_import
from collections import defaultdict
import os
import click
import ncls
import numpy as np
from . import samclips


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory {}".format(d))


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

    if "outstr" in template:
        return template["outstr"]

    # Todo make sure all read-pairs have a mapping, otherwise write an unmapped
    #
    # paired = False if template["read2_length"] is None else True
    sam = samclips.fixsam(template)

    if len(sam) == 0:  # Todo fix unmapped reads
        template["passed"] = 0
        # print("No alignments")
        # click.echo(template, err=True)
        # click.echo("Unmapped read error", err=True)
        # click.echo(template["inputdata"], err=True)
        # quit()
        return ""

    # if len(sam) == 1:  # Todo deal with these
    #     click.echo("Single alignment only", err=True)
    #    return None

    if any(i[0] == "*" or i[4] == "*" for i in sam):
        template["passed"] = 0
        # click.echo("Unformatted cigar error", err=True)
        # click.echo(sam, err=True)
        # quit()
        return ""

    return "".join(template["name"] + "\t" + "\t".join(i) + "\n" for i in sam)


def get_bed_regions(bed):
    return [tuple([int(j) if j.isdigit() else j for j in i.strip().split("\t")[:3]]) for i in open(bed, "r") if i[0] != "#"]


def overlap_regions(bed):
    if not bed:
        return None
    regions = get_bed_regions(bed)
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


def get_include_reads(include, bam):

    if not include:
        for r in bam:
            yield r

    regions = [i.strip().split("\t")[:3] for i in open(include, "r") if i[0] != "#"]
    for c, s, e in regions:
        click.echo("Reading {}:{}-{}".format(c, s, e), err=True)
        for r in bam.fetch(c, int(s), int(e)):
            yield r


def sam_itr(args):

    itr = args["sam"]
    tree = overlap_regions(args["include"])

    # First get header
    header_string = ""
    last_seen_chrom = ""
    first_line = ""
    for t in itr:

        # t = str(t.decode("ascii"))

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
        # t = str(t.decode("ascii"))
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
            l1 = l1.strip()
            last2 = l1[-2:]
            if last2 == "/2" or last2 == "/1":
                l1 = l1[1:-2]  # Strip trailing /1 or /2 and leading @
            l2 = next(f).strip()
            next(f)  # Skip "+"
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


def iterate_mappings(args, version):
    arg_str = ", ".join(["{}={}".format(i, j) for i, j in args.items() if j is not None and "sam" not in i])
    inputstream = sam_itr(args)

    total = 0
    name = ""
    rows = []
    header_string = next(inputstream)
    header_string += "@PG\tID:fnfi\tPN:fnfi align\tVN:{}\tCL:{}\n".format(version, arg_str)

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

                # if name == "HISEQ2500-10:539:CAV68ANXX:7:2211:10642:81376":
                #     click.echo("io", err=True)

                yield (rows, args, max_d, last_seen_chrom, fq)

            rows = []
            name = nm

        rows.append((m, ol))  # String, ol states if alignment overlaps ROI

    # Deal with last record
    if len(rows) > 0:
        total += 1
        fq = fq_getter(fq_iter, name, args, fq_buffer)
        yield (rows, args, max_d, last_seen_chrom, fq)

    click.echo("Total processed " + str(total), err=True)


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
    sam_ol = [["\t".join(map(str, i)).split("\t", 4), False] for i in sam]  # Identifier for overlapping region of interest

    args = {"max_insertion": 100, "min_aln": 17, "max_overlap": 100, "ins_cost": 1, "ol_cost": 3, "inter_cost":2,
            "u": 9, "match_score": 1, "insert_median": 500, "insert_stdev": 50, "bias": 1.15, "paired": "True", "replace_hardclips": False, "fq1": None, "fq2": None}


    template = make_template(sam_ol, args, 2000, "chr22", (0, 0))

    # t0 = time.time()
    # sam_to_array(template)
    # v0 = time.time() - t0
    # print template["data"].astype(int)
    template = {'inputfq': (('chr5.76630308-76630338.5272014-5272109:76630645-76630770.76630645-76630770', 'ACACTTGCAATGAAAAGGGGGGAGCAATTTCCTTTTCTGTAGTTGCTTATCTCTTTTCGATCTGAGGCCGTGAAAATACTAATCACCCAGTAATGGCACAGCACAATTTTAATGACAAGGACTCT', 'EFFGG>FAGFFGFFDGGFGGEEGBGF8GFGGEGGDGEDGAGFECF/DFGEGFGGF$;GG??EDGFGEGFGGGF9EGFDGG<@AFGF?BFBGFGGGGGGGG=1AED>F:EE+=GEGCA6(4>C<GC'), ('chr5.76630308-76630338.5272014-5272109:76630645-76630770.76630645-76630770', 'GTAGACATAGCCCTTTCCTCTTCTCCAGGAAAGCTTTGTAGGATAGAGAAAATAAGAGTACTGGCAATTGAGAAGTTTGTGGTCCAAGAAAAGTAATAGGAACAAGAAAATGACAAAAATTATCA', 'FGGGDEGGGDCGGCAGGF4;GFEGFFGAFGGCGFGGDGF=FGDGEGEBFGDGGE2GGGFG=GGFGE??E@GFGGFCGGGBGGFGGFD;GDE,2>FEA=GGF.,=AG4C?.FGD@F=:;;A4DEF$')), 'paired_end': 1, 'bias': 1.0, 'fq_read2_seq': 0, 'isize': (210.0, 175.0), 'read2_q': 0, 'max_d': 1085.0, 'read2_seq': 0, 'read2_length': 0, 'passed': 0, 'replace_hard': 1, 'read2_reverse': 0, 'inputdata': [(['chr5.76630308-76630338.5272014-5272109:76630645-76630770.76630645-76630770', '99', 'chr5', '5272015', '254\t95=30H\t*\t0\t0\tACACTTGCAATGAAAAGGGGGGAGCAATTTCCTTTTCTGTAGTTGCTTATCTCTTTTCGATCTGAGGCCGTGAAAATACTAATCACCCAGTAATG\tEFFGG>FAGFFGFFDGGFGGEEGBGF8GFGGEGGDGEDGAGFECF/DFGEGFGGF$;GG??EDGFGEGFGGGF9EGFDGG<@AFGF?BFBGFGGG\tNM:i:0\tAS:i:94\tEV:Z:4.6e-46\n'], 0), (['chr5.76630308-76630338.5272014-5272109:76630645-76630770.76630645-76630770', '99', 'chr5', '76630309', '0\t95H30=\t*\t0\t0\tGCACAGCACAATTTTAATGACAAGGACTCT\tGGGGG=1AED>F:EE+=GEGCA6(4>C<GC\tNM:i:0\tAS:i:30\tEV:Z:4.1e-07\n'], 0), (['chr5.76630308-76630338.5272014-5272109:76630645-76630770.76630645-76630770', '147', 'chr5', '76630647', '254\t1H124=\t*\t0\t0\tGATAATTTTTGTCATTTTCTTGTTCCTATTACTTTTCTTGGACCACAAACTTCTCAATTGCCAGTACTCTTATTTTCTCTATCCTACAAAGCTTTCCTGGAGAAGAGGAAAGGGCTATGTCTAC\tFED4A;;:=F@DGF.?C4GA=,.FGG=AEF>2,EDG;DFGGFGGBGGGCFGGFG@E??EGFGG=GFGGG2EGGDGFBEGEGDGF=FGDGGFGCGGFAGFFGEFG;4FGGACGGCDGGGEDGGGF\tNM:i:0\tAS:i:124\tEV:Z:1.2e-65\n'], 0)], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'read1_q': 0, 'read1_length': 0, 'name': 'chr5.76630308-76630338.5272014-5272109:76630645-76630770.76630645-76630770', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 0, 'last_seen_chrom': 'chr5', 'score_mat': {}, 'pairing_params': (100.0, 17.0, 100.0, 1.0, 3.0, 2.0, 9.0)}

    t0 = time.time()
    c_io_funcs.sam_to_array(template)
    v = time.time() - t0
