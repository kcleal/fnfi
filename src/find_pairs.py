"""
Get all read pairs which fall within these regions, or have a read anchored in a region.
Outputs paired reads in .fastq format for re-mapping
Note regions are intersected by only the mapping 'position' of the read in the .bam file, the read-length is not used.
"""

import pysam
import time
import click
import datetime
import numpy as np
from pybedtools import BedTool
import os
import quicksect
from collections import defaultdict
import click
from subprocess import call


def make_tree(bed_path):
    b = BedTool(bed_path)
    tree = defaultdict(lambda: quicksect.IntervalTree())
    for item in b:
        tree[item.chrom].add(item.start, item.end)
    return tree


def get_reads(args):

    if args["search"]:
        inc_tree = make_tree(args["search"])
    if args["exclude"]:
        exc_tree = make_tree(args["exclude"])

    bam = pysam.AlignmentFile(args["bam"], "rb")

    clip_length = args["clip_length"]
    if args["dest"]:
        out_name = args["dest"] + "/" + args["bam"].rsplit("/", 1)[1].rsplit(".", 1)[0] + "." + args["post_fix"]
        if not os.path.exists(args["dest"]):
            os.makedirs(args["dest"])
    else:
        out_name = args["bam"].rsplit(".", 1)[0] + "." + args["post_fix"]

    if args["mapper"] == "last":
        call("samtools view -H -o {}.dict {}".format(out_name, args["bam"]), shell=True)

    read_names = set([])
    insert_size = []
    read_length = []
    for r in bam.fetch(until_eof=True):  # Also get unmapped reads

        chrom = bam.get_reference_name(r.rname)
        if args["exclude"]:  # Skip exclude regions
            if chrom in exc_tree and len(exc_tree[chrom].search(r.pos, r.pos + 1)) > 0:
                continue

        if args["search"]:  # Limit search to regions
            if chrom not in inc_tree:
                continue
            elif len(inc_tree[chrom].search(r.pos, r.pos + 1)) == 0:
                continue

        if r.flag & 1280:
            continue  # Read is duplicate or not primary alignment
        if not r.cigartuples:  # Cigar is not formatted
            continue

        if len(insert_size) < 10000 and r.flag & 2:
            tmpl = abs(r.template_length)
            if tmpl < 10000:
                insert_size.append(tmpl)
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
        if r.qname in read_names and not r.flag & 2304:  # Skip not primary, and supplementary reads
            outbam.write(r)
    outbam.close()

    # Save the insert size and read length for later
    insert_median, insert_stdev = np.median(insert_size), np.std(insert_size)
    read_length = np.mean(read_length)
    click.echo("Median insert size: {} (+/- {})".format(np.round(insert_median, 2), np.round(insert_stdev, 2)), err=True)
    click.echo("Read length: {}".format(np.round(read_length, 1)), err=True)

    return insert_median, insert_stdev, read_length, out_name


def convert_to_fastq(args, outname):

    pysam.sort(*("-n -@{} -o {} {}".format(args["procs"], outname + ".srt.bam", outname + ".bam").split(" ")))
    BedTool(outname + ".srt.bam").bam_to_fastq(fq=outname + "1.fq", fq2=outname + "2.fq")

    os.remove(outname + ".srt.bam")
    os.remove(outname + ".bam")


def process(args):
    t0 = time.time()
    insert_median, insert_stdev, read_length, out_name = get_reads(args)
    convert_to_fastq(args, out_name)
    click.echo("Collected reads in {} h:m:s\n".format(str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)
    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            "fastq": out_name,
            "out_pfix": out_name}

