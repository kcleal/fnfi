import datetime
import os
import time
from multiprocessing import cpu_count
from subprocess import Popen, PIPE, check_call

import click

from src import find_pairs, cluster
from align import input_stream_alignments

cpu_range = click.IntRange(min=1, max=cpu_count())


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory")

# ----------------------------------------------------------------------------------------------------------------------


@click.group(chain=False, invoke_without_command=False)
def cli():
    """Fusion-finder calls structural variants from input .bam or .fastq files."""
    pass


# ----------------------------------------------------------------------------------------------------------------------
@cli.command("run")
@click.argument('reference', required=True, type=click.Path(exists=False))
@click.argument("bam", required=True, type=click.Path(exists=True))
@click.option('--include', help=".bed file, limit calls to regions", default=None, type=click.Path(exists=True))
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
              default=None, type=click.Path(exists=True))
@click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept", default=21, type=int,
              show_default=True)
@click.option('--mapper', help="External mapper to use for re-alignment", default="bwamem",
              type=click.Choice(['bwamem', 'last']), show_default=True)
@click.option('--map-script', help="""External shell script for mappping fastq files. \
Overrides --mapper argument. Script must take positional arguments as: $1 reference genome; \
$2 read1.fastq file $3 read2.fastq if paired-end \
$4 threads to use""", default=None, type=click.Path(exists=True))
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--dest', help="Destination folder to use/create for saving results. Defaults to directory of input bam",
              default=None, type=click.Path())
@click.pass_context
def run_command(ctx, **kwargs):
    """Run the fusion-finder pipeline."""
    t0 = time.time()
    click.echo("Running fufi pipeline", err=True)
    if not os.path.exists(kwargs["bam"] + ".bai"):
        raise IOError("Input .bai index file not found.")
    mk_dest(kwargs["dest"])

    other_kwargs = ctx.forward(find_reads, kwargs)
    kwargs.update(other_kwargs)
    single = True if kwargs["procs"] == 1 else False

    process, used_procs, remaining_procs = launch_external_mapper(kwargs)
    kwargs["procs"] = remaining_procs
    kwargs["fq1"] = kwargs["out_pfix"] + "1.fq"
    kwargs["fq2"] = kwargs["out_pfix"] + "2.fq"

    ctx.invoke(fufi_aligner, sam=process.stdout, output=kwargs["out_pfix"] + ".sam", **kwargs)
    process.kill()
    # Clean up
    os.remove(kwargs["fq1"])
    os.remove(kwargs["fq2"])
    if kwargs["mapper"] == "last":
        os.remove(kwargs["out_pfix"] + ".dict")

    if not single:
        kwargs["procs"] = remaining_procs + used_procs
    sort_and_index(kwargs)

    ctx.forward(call_events, raw_bam=kwargs["bam"], sv_bam=kwargs["out_pfix"] + ".srt.bam", **kwargs)

    click.echo("fufi run completed in {} h:m:s\n".format(str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)


# ----------------------------------------------------------------------------------------------------------------------
@cli.command("find-reads")
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.option('--post-fix', help="Post fix to tag temp files with. Default is to use 'fufi'", default='fufi', type=str)
@click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept", default=21, type=int,
              show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--search', help=".bed file, limit search to regions.", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
              default=None, type=click.Path(exists=True))
@click.option('--dest', help="Destination folder to use/create for saving results. Defaults to directory of input bam",
              default=None, type=click.Path())
def find_reads(**kwargs):
    """Filters input .bam for read-pairs that are discordant or have a soft-clip of length >= '--clip-length'"""
    # insert_median, insert_stdev, read_length, out_name
    mk_dest(kwargs["dest"])
    return find_pairs.process(kwargs)


# ----------------------------------------------------------------------------------------------------------------------
def launch_external_mapper(kwargs):
    """Run a shell script to map interleaved .fastq files to .sam format. Uses a basic bwa-mem by default.
    A custom shell script can be provided but must take positional arguments as:
    $1 reference genome
    $2 .fastq file; interleaved if paired-end, otherwise single end reads"""

    if kwargs["procs"] > 1:
        p = int(kwargs["procs"] / 2)
        other_p = kwargs["procs"] - p
    else:
        p = kwargs["procs"]
        other_p = 1

    if not kwargs["map_script"] and kwargs["mapper"] == "bwamem":
        command = "bwa mem -t {procs} -P -a {ref} {s}1.fq {s}2.fq".format(procs=p,
                                                                             ref=kwargs["reference"],
                                                                             s=kwargs["fastq"])

    elif not kwargs["map_script"] and kwargs["mapper"] == "last":
        # Todo need to exract the sam header from input file, and use this as the dict argument in maf-convert
        command = "fastq-interleave {s}1.fq {s}2.fq \
        | lastal -k2 -l11 -Q1 -D10000 -K8 -C8 -i10M -r1 -q4 -a6 -b1 -P{procs} {ref} \
        | last-map-probs -m 1 -s 1 | maf-convert -f {d}.dict sam".format(procs=p,
                                              ref=kwargs["reference"],
                                              d=kwargs["out_pfix"],
                                              s=kwargs["fastq"])

    else:
        command = "bash {script} {ref} {s}1.fq {s}2.fq {procs}".format(script=kwargs["map_script"],
                                                                       procs=p,
                                                                       ref=kwargs["reference"],
                                                                       s=kwargs["fastq"])

    click.echo("Mapping command:\n" + command, err=True)
    proc = Popen(command, stdout=PIPE, shell=True)

    return proc, p, other_p


# ----------------------------------------------------------------------------------------------------------------------
@cli.command("align")
@click.argument("sam", type=click.File('r'), required=True)
@click.argument("output", required=False, type=click.Path())
@click.option("--insert-median", help="Template insert size", default=300, type=float)
@click.option("--insert-stdev",  help="Template standard-deviation", default=100, type=float)
@click.option("--read-length",  help="Length of a read in base-pairs", default=100, type=float)
@click.option("--fq1",  help="Fastq reads 1, used to add soft-clips to all hard-clipped read 1 alignments",
              default=None, type=click.Path(exists=True))
@click.option("--fq2",  help="Fastq reads 2, used to add soft-clips to all hard-clipped read 2 alignments",
              default=None, type=click.Path(exists=True))
@click.option("--max_insertion", help="Maximum insertion within read", default=100, type=int)
@click.option("--min-aln", help="Minimum alignment length", default=17, type=int)
@click.option("--max-overlap", help="Maximum overlap between successive alignments", default=100, type=int)
@click.option("--ins-cost", help="Insertion cost", default=1, type=float)
@click.option("--ol-cost", help="Overlapping alignment cost", default=3, type=int)
@click.option("--inter-cost", help="Cost of inter-chromosomal jump", default=2, type=int)
@click.option("--u", help="Pairing heuristic cost", default=9, type=int)
@click.option("--match-score", help="Matched base score used for input sam reads", default=1, type=int)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1)
@click.option('--include', help=".bed file, elevate alignment scores in these regions. Determined by '--bias'",
              default=None, type=click.Path(exists=True))
@click.option("--bias", help="""Multiply match score by bias if alignment falls within regions .bed file.
Unused if .bed not provided.""", default=1.15, type=float)
def fufi_aligner(**kwargs):
    """Choose an optimal set of alignments from from a collection of candidate alignments.
    If reads are paired, alignments must be sorted by read-name while the bit flag
    designates read 1 vs read 2."""
    if not kwargs["include"]:
        kwargs["bias"] = 1.0
    else:
        click.echo("Elevating alignments in --include with --bias {}".format(kwargs["bias"]), err=True)
    input_stream_alignments.process_reads(kwargs)


def sort_and_index(kwargs):
    """Convenience function to sort and index a sam file, then remove the input sam file"""
    c = "samtools view -uh {fix}.sam | samtools sort -@ {p} -o {fix}.srt.bam - ; samtools index -@ {p} {fix}.srt.bam"
    c = c.format(fix=kwargs["out_pfix"], p=kwargs["procs"])
    click.echo(c, err=True)
    check_call(c, shell=True)
    os.remove(kwargs["out_pfix"] + ".sam")


# ----------------------------------------------------------------------------------------------------------------------
@cli.command("call-events")
@click.argument('raw-bam', required=True, type=click.Path(exists=True))
@click.argument('sv-bam', required=True, type=click.Path(exists=True))
@click.argument("output", required=False, type=click.Path())
@click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept.", default=21, type=int,
              show_default=True)
@click.option("--insert-median", help="Template insert size", default=300, type=float)
@click.option("--insert-stdev",  help="Template standard-deviation", default=100, type=float)
@click.option("--read-length",  help="Length of a read in base-pairs", default=100, type=float)
@click.option('--verbose', type=click.Choice(["True", "False"]), default="True", help="If set to 'True' output is directed to a folder with the same name as the input file. Otherwise a .vcf file is generated.")
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--include', help=".bed file, limit calls to regions.", default=None, type=click.Path(exists=True))
@click.option('--dest', help="Folder to use/create for saving results. Defaults to directory of input bam directory",
              default=None, type=click.Path())
def call_events(**kwargs):
    """Clusters reads into SV-events. Takes as input the original .bam file, and a .bam file with only sv-like reads."""
    # Create dest in not done so already
    mk_dest(kwargs["dest"])
    cluster.cluster_reads(kwargs)

