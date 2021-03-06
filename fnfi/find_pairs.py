"""
Get all read pairs which fall within these regions, or have a read anchored in a region.
Outputs paired reads in .fastq format for re-mapping
Note regions are intersected by only the mapping 'position' of the read in the .bam file, the read-length is not used.
"""

from __future__ import absolute_import
import pysam
import time
import datetime
import numpy as np
import os
import click
from subprocess import call, Popen, PIPE
from . import data_io


def iter_bam(bam, search):

    if not search:
        click.echo("Searching input file", err=True)
        for aln in bam.fetch(until_eof=True):  # Also get unmapped reads
            yield aln
    else:
        click.echo("Limiting search to {}".format(search), err=True)
        with open(search, "r") as bed:
            for line in bed:
                if line[0] == "#":
                    continue
                chrom, start, end = line.split("\t")[:3]
                for aln in bam.fetch(reference=chrom, start=int(start), end=int(end)):  # Also get unmapped reads
                    yield aln


def get_reads(args):
    kind = args["bam"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "rs"}
    click.echo("Input file is {}, (.{} format)".format(args["bam"], kind), err=True)

    bam = pysam.AlignmentFile(args["bam"], opts[kind])

    bam_i = iter_bam(bam, args["search"])

    if args["exclude"]:
        click.echo("Excluding {} from search".format(args["exclude"]), err=True)
        exc_tree = data_io.overlap_regions(args["exclude"])

    clip_length = args["clip_length"]
    if args["dest"]:

        if "/" in args["bam"]:
            out_name = args["dest"] + "/" + args["bam"].rsplit("/", 1)[1].rsplit(".", 1)[0] + "." + args["post_fix"]
        else:
            out_name = args["dest"] + "/" + args["bam"].rsplit(".", 1)[0] + "." + args["post_fix"]
        if not os.path.exists(args["dest"]):
            os.makedirs(args["dest"])
    else:
        if "/" in args["bam"]:
            out_name = os.getcwd() + "/" + args["bam"].rsplit("/", 1)[1].rsplit(".", 1)[0] + "." + args["post_fix"]
        else:
            out_name = os.getcwd() + "/" + args["bam"].rsplit(".", 1)[0] + "." + args["post_fix"]

    if args["mapper"] == "last":
        call("samtools view -H -o {}.dict {}".format(out_name, args["bam"]), shell=True)

    read_names = set([])

    insert_size = []
    read_length = []

    for r in bam_i:

        chrom = bam.get_reference_name(r.rname)
        if args["exclude"]:  # Skip exclude regions
            if data_io.intersecter(exc_tree, chrom, r.pos, r.pos + 1):
                continue

        if r.flag & 1280:
            continue  # Read is duplicate or not primary alignment
        if not r.cigartuples:  # Cigar is not formatted
            continue

        if len(insert_size) < 10000:
            if r.flag & 2:
                tmpl = abs(r.template_length)
                if tmpl < 1000:
                    insert_size.append(tmpl)
        if len(read_length) < 10000:
            read_length.append(r.infer_read_length())

        if r.qname not in read_names:
            if not r.flag & 2 or r.flag & 2048:  # Save if read is discordant or supplementary exists
                read_names.add(r.qname)
            else:
                softclips = [i[1] for i in r.cigartuples if i[0] == 4]
                if len(softclips) != 0 and max(softclips) >= clip_length:
                    read_names.add(r.qname)

    if len(read_names) == 0:
        click.echo("No reads found, aborting", err=True)
        quit()

    click.echo("Retrieving {} query names into {}".format(len(read_names), out_name + ".bam"), err=True)
    outbam = pysam.AlignmentFile(out_name + ".bam", "wb", template=bam)

    for r in bam.fetch():
        if r.qname in read_names and not r.flag & 3328:  # Skip not primary, and supplementary reads, and duplicates
            outbam.write(r)

    outbam.close()

    if len(insert_size) == 0:
        insert_median, insert_stdev = args["insert_median"], args["insert_stdev"]
        click.echo("WARNING: could not infer insert size, no 'normal' pairings found. Using arbitrary values", err=True)
    else:
        insert_median, insert_stdev = np.median(insert_size), np.std(insert_size)
        insert_median, insert_stdev = np.round(insert_median, 2), np.round(insert_stdev, 2)

    click.echo("Median insert size: {} (+/- {})".format(insert_median, insert_stdev), err=True)

    return insert_median, insert_stdev, int(np.mean(read_length)), out_name


def convert_to_fastq(outname):

    click.echo("Converting to fastq {}".format(outname), err=True)
    samtools_version = None
    try:
        samtools_version = Popen("samtools --version", shell=True, stderr=PIPE, stdout=PIPE).communicate()[0]
        samtools_version = float(str(samtools_version).split("\\n")[0].split(" ")[1])
    except:
        # Couldnt parse samtools version number
        pass

    if samtools_version is not None and samtools_version >= 1.9:
        call("samtools collate -r 100000 -f -o {o} {i}".format(o=outname + ".srt.bam",
                                                               i=outname + ".bam"), shell=True)
    else:
        call("samtools sort -n -o {o}.srt.bam {i}.bam".format(o=outname, i=outname), shell=True)

    call("samtools fastq -s /dev/null -1 {fq1} -2 {fq2} {bam}".format(fq1=outname + "1.fq",
                                                                      fq2=outname + "2.fq",
                                                                      bam=outname + ".srt.bam"), shell=True)

    os.remove(outname + ".srt.bam")
    os.remove(outname + ".bam")


def process(args):
    t0 = time.time()

    data_io.mk_dest(args["dest"])

    insert_median, insert_stdev, read_length, out_name = get_reads(args)

    if args["bam_only"] is None:
        convert_to_fastq(out_name)
    else:
        call("samtools index {}.bam".format(out_name), shell=True)

    click.echo("Collected reads in {} h:m:s".format(str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)

    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            "fastq": out_name,
            "out_pfix": out_name}

