import click
import datetime
import os
import time
from multiprocessing import cpu_count
from subprocess import Popen, PIPE, check_call, call
import find_pairs
import cluster
import input_stream_alignments
import data_io

cpu_range = click.IntRange(min=1, max=cpu_count())

defaults = {
            "clip_length": 21.,
            "mapper": "bwamem",
            "map_script": None,
            "procs": 1,
            "dest": None,
            "post_fix": "fnfi",
            "search": None,
            "exclude": None,
            "include": None,
            "paired": "True",
            # "insert_median": 210.,
            # "insert_stdev": 175.,
            "read_length": 125.,
            "max_insertion": 150.,
            "min_aln": 17.,
            "max_overlap": 150.,
            "ins_cost": 0.1,
            "ol_cost": 1.,
            "inter_cost": 2.,
            "u": 9.,
            "match_score": 1.,
            "bias": 1.15,
            "output": "-",
            "svs_out": "-",
            "replace_hardclips": "False",
            "fq1": None,
            "fq2": None,
            "max_cov": 150,
            "buffer_size": 1000000,
            "I": "210,175"
            }

align_args = {}


def pipeline(kwargs):  # Todo add a verbosity option (keep some or all of temp files)
    t0 = time.time()
    click.echo("Running fnfi pipeline", err=True)
    if kwargs["bam"] is None:
        raise IOError("Error: Input bam is None")

    if not os.path.exists(kwargs["bam"] + ".bai"):
        raise IOError("Input .bai index file not found.")

    data_io.mk_dest(kwargs["dest"])

    other_kwargs = find_pairs.process(kwargs)
    kwargs.update(other_kwargs)

    single = True if kwargs["procs"] == 1 else False

    process, used_procs, remaining_procs = launch_external_mapper(kwargs)
    kwargs["procs"] = remaining_procs

    if kwargs["mapper"] == "bwamem" and kwargs["map_script"] is None:
        kwargs["fq1"] = None
        kwargs["fq2"] = None
    else:
        kwargs["fq1"] = kwargs["out_pfix"] + "1.fq"
        kwargs["fq2"] = kwargs["out_pfix"] + "2.fq"

    # os.fdopen(process.stdout.fileno(), "r", 1)  # Alyways get this error close failed in file object destructor
    kwargs["sam"] = iter(process.stdout.readline, "")
    kwargs["output"] = kwargs["out_pfix"] + ".sam"

    input_stream_alignments.process_reads(kwargs)

    process.kill()

    if not single:
        kwargs["procs"] = remaining_procs + used_procs - 1

    sort_and_index(kwargs)

    # Clean up
    if kwargs["fq1"] is not None:
        os.remove(kwargs["fq1"])
        os.remove(kwargs["fq2"])

    if kwargs["mapper"] == "last":
        os.remove(kwargs["out_pfix"] + ".dict")

    kwargs["sv_aligns"] = kwargs["out_pfix"] + ".srt.bam"
    kwargs["raw_aligns"] = kwargs["bam"]

    cluster.cluster_reads(kwargs)
    # Todo cleanup of other temp files
    click.echo("fnfi run {} completed in {} h:m:s\n".format(kwargs["bam"],
                                                            str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)


def launch_external_mapper(kwargs):
    """Run a shell script to map interleaved .fastq files to .sam format. Uses a basic bwa-mem by default.
    A custom shell script can be provided but must take positional arguments as:
    $1 reference genome
    $2 .fastq file; interleaved if paired-end, otherwise single end reads"""

    # if kwargs["procs"] > 1:
    #     p = int(kwargs["procs"]) - 1
    #     other_p = 1  # kwargs["procs"] - p
    # else:
    #     p = kwargs["procs"]
    #     other_p = 1

    p = int(kwargs["procs"])
    if not kwargs["map_script"] and kwargs["mapper"] == "bwamem":
        command = "bwa mem -Y -t {procs} -a {ref} {s}1.fq {s}2.fq".format(procs=p - 1 if p > 1 else 1,
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

    proc = Popen(command, stdout=PIPE, shell=True, bufsize=0, universal_newlines=True)

    return proc, p, 1  # Note fnfi align always has 1 core assigned


def sort_and_index(kwargs):
    """Convenience function to sort and index a sam file, then remove the input sam file"""
    c = "samtools view -uh {fix}.sam | samtools sort -@ {p} -o {fix}.srt.bam - ; samtools index -@ {p} {fix}.srt.bam"
    c = c.format(fix=kwargs["out_pfix"], p=kwargs["procs"])
    click.echo(c, err=True)
    check_call(c, shell=True)
    os.remove(kwargs["output"])


# Interface:
# ----------------------------------------------------------------------------------------------------------------------
def apply_ctx(ctx, kwargs):
    ctx.ensure_object(dict)
    if len(ctx.obj) == 0:  # When run is invoked from cmd line, else run was invoked from test function
        for k, v in defaults.items() + kwargs.items():
            ctx.obj[k] = v
    i, j = map(float, ctx.obj["I"].split(","))
    if "insert_median" not in ctx.obj:
        ctx.obj["insert_median"] = i
        ctx.obj["insert_stdev"] = j

    return ctx


@click.group(chain=False, invoke_without_command=False)
@click.version_option()
def cli():
    """Fusion-finder calls structural variants from input .bam or .fastq files."""
    pass


@cli.command("run")
@click.argument('reference', required=True, type=click.Path(exists=False))
@click.argument("bam", required=True, type=click.Path(exists=True))
@click.option('--include', help=".bed file, limit calls to regions", default=None, type=click.Path(exists=True))
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
              default=None, type=click.Path(exists=True))
@click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept", default=defaults["clip_length"], type=int,
              show_default=True)
@click.option('--mapper', help="External mapper to use for re-alignment", default=defaults["mapper"],
              type=click.Choice(['bwamem', 'last']), show_default=True)
@click.option('--map-script', help="""External shell script for mappping fastq files. \
Overrides --mapper argument. Script must take positional arguments as: $1 reference genome; \
$2 read1.fastq file $3 read2.fastq if paired-end \
$4 threads to use""", default=None, type=click.Path(exists=True))
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--dest', help="Destination folder to use/create for saving results. Defaults to current directory",
              default=None, type=click.Path())
@click.option('-I', help="Insert size and stdev as 'FLOAT,FLOAT'. If not provided, automatically inferred",
              default=defaults["I"], type=str)
@click.pass_context
def run_command(ctx, **kwargs):
    """Run the fusion-finder pipeline."""
    ctx = apply_ctx(ctx, kwargs)
    pipeline(ctx.obj)


@cli.command("find-reads")
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.option('--post-fix', help="Post fix to tag temp files with. Default is to use 'fnfi'", default='fnfi', type=str)
@click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept", default=defaults["clip_length"], type=int,
              show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=defaults["procs"], show_default=True)
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
              default=None, type=click.Path(exists=True))
@click.option('--dest', help="Destination folder to use/create for saving results. Defaults to current directory",
              default=None, type=click.Path())
@click.pass_context
def find_reads(ctx, **kwargs):
    """Filters input .bam for read-pairs that are discordant or have a soft-clip of length >= '--clip-length'"""
    # Add arguments to context insert_median, insert_stdev, read_length, out_name
    ctx = apply_ctx(ctx, kwargs)
    return find_pairs.process(ctx.obj)


@cli.command("align")
@click.argument("sam", type=click.File('r', encoding="ascii"), required=True)
@click.argument("output", required=False, type=click.Path())
@click.option("--paired", help="Paired end reads or single", default=defaults["paired"],
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-I', help="Insert size and stdev as 'FLOAT,FLOAT'",
              default=defaults["I"], type=str, show_default=True)
@click.option("--replace-hardclips",  help="Replace hard-clips with soft-clips when possible", default=defaults["replace_hardclips"], type=click.Choice(["True", "False"]), show_default=True)
@click.option("--fq1",  help="Fastq reads 1, used to add soft-clips to all hard-clipped read 1 alignments",
              default=defaults["fq1"], type=click.Path(), show_default=True)
@click.option("--fq2",  help="Fastq reads 2, used to add soft-clips to all hard-clipped read 2 alignments",
              default=defaults["fq2"], type=click.Path(), show_default=True)
@click.option("--max-insertion", help="Maximum insertion within read", default=defaults["max_insertion"], type=float, show_default=True)
@click.option("--min-aln", help="Minimum alignment length", default=defaults["min_aln"], type=float, show_default=True)
@click.option("--max-overlap", help="Maximum overlap between successive alignments", default=defaults["max_overlap"], type=float, show_default=True)
@click.option("--ins-cost", help="Insertion cost", default=defaults["ins_cost"], type=float, show_default=True)
@click.option("--ol-cost", help="Overlapping alignment cost", default=defaults["ol_cost"], type=float)
@click.option("--inter-cost", help="Cost of inter-chromosomal jump", default=defaults["inter_cost"], type=float, show_default=True)
@click.option("--u", help="Pairing heuristic cost", default=defaults["u"], type=float, show_default=True)
@click.option("--match-score", help="Matched base score used for input sam reads", default=defaults["match_score"], type=float, show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--include', help=".bed file, elevate alignment scores in these regions. Determined by '--bias'",
              default=None, type=click.Path(exists=True))
@click.option("--bias", help="""Multiply match score by bias if alignment falls within regions .bed file
Unused if .bed not provided""", default=defaults["bias"], type=float, show_default=True)
@click.pass_context
def fnfi_aligner(ctx, **kwargs):
    """Choose an optimal set of alignments from from a collection of candidate alignments.
    If reads are paired, alignments must be sorted by read-name with the bit flag
    designating read_1 vs read_2."""
    ctx = apply_ctx(ctx, kwargs)
    input_stream_alignments.process_reads(ctx.obj)


@cli.command("call-events")
# @click.argument('raw-aligns', required=True, type=click.Path(exists=True))
@click.argument('sv-aligns', required=True, type=click.Path(exists=True))
@click.argument("svs-out", required=False, type=click.Path())
@click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept.", default=defaults["clip_length"], type=int,
              show_default=True)
@click.option('--max-cov', help="Regions with > max-cov that do no overlap 'include' are discarded.", default=defaults["max_cov"], type=float,
              show_default=True)
@click.option('-I', help="Insert size and stdev as 'FLOAT,FLOAT'",
              default=defaults["I"], type=str, show_default=True)
@click.option('--verbose', type=click.Choice(["True", "False"]), default="True", help="If set to 'True' output is directed to a folder with the same name as the input file. Otherwise a .vcf file is generated")
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=defaults["procs"], show_default=True)
@click.option('--include', help=".bed file, limit calls to regions", default=None, type=click.Path(exists=True))
@click.option('--dest', help="Folder to use/create for saving results. Defaults to current directory",
              default=None, type=click.Path())
@click.option("--buffer-size", help="Number of alignments to load into buffer", default=defaults["buffer_size"],
              type=int, show_default=True)
@click.pass_context
def call_events(ctx, **kwargs):
    """Clusters reads into SV-events. Takes as input the original .bam file, and a .bam file with only sv-like reads."""
    # Create dest in not done so already
    ctx = apply_ctx(ctx, kwargs)
    cluster.cluster_reads(ctx.obj)


@cli.command("run-tests", context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option('--dest', required=True,  help="Folder to use/create for saving results",
              default=None, type=click.Path())
@click.option('--reference', required=True, type=click.Path(),  help="Path to reference for mapper to use")
@click.pass_context
def test_run_command(ctx, **kwargs):
    """Runs fnfi commands on test data using default options.
    fnfi defaults for the run command can be modified by adding the respective option, e.g. add --procs 4"""

    def update_ctx(kwargs, ctx):
        for k, v in defaults.items() + kwargs.items() + {ctx.args[i][2:]: ctx.args[i + 1] for i in
                                                         xrange(0, len(ctx.args), 2)}.items():
            if isinstance(v, str):
                if v.isdigit():
                    if int(v) == float(v):
                        v = int(v)
                    else:
                        v = float(v)
            ctx.obj[k] = v

    ctx.ensure_object(dict)

    data_io.mk_dest(kwargs["dest"])
    tests_path = os.path.dirname(__file__) + "/tests"
    click.echo("Input data for tests: {}".format(tests_path), err=True)

    # Test aligner using single reads
    f1 = "{}/long_contigs.fa".format(tests_path)
    success = False
    try:
        com = "lastal -D1000 -K3 -C3 -r1 -q4 -a6 -b1 -P1 {ref} {f} \
         | maf-convert -f {tp}/hg38.dict sam > {dest}/long_contigs.last.sam".format(tp=tests_path, f=f1,
                                              ref=kwargs["reference"],
                                              dest=kwargs["dest"],
                                              )
        call(com, shell=True)
        success = True

    except OSError as e:
        click.echo("Skipping align test using longer contigs", err=True)
        if e.errno == os.errno.ENOENT:
            click.echo("No such file or directory", err=True)
        else:
            raise SystemError("Something went wrong")
    if success:
        sam = "{dest}/long_contigs.last.sam".format(dest=kwargs["dest"])
        output = "{dest}/long_contigs.fnfi.sam".format(dest=kwargs["dest"])

        update_ctx(kwargs, ctx)

        call("fnfi align --paired False {sam} > {o}".format(sam=sam, o=output), shell=True)


    quit()

    # Test run command using paired-end data
    b1 = "{}/bwa.0.5.srt.bam".format(tests_path)

    click.echo(b1)
    inc = "{}/include_tels.bed".format(tests_path)

    kwargs["bam"] = b1
    kwargs["include"] = inc

    update_ctx(kwargs, ctx)
    click.echo(ctx.obj)
    ctx.invoke(run_command)


if __name__ == "__main__":

    # df = defaults
    # df["sam"] = open("/Users/kezcleal/Documents/Data/fusion_finder_development/Split_read_simulator/paired_end_with_one_split/pairs.bwamem_allf.sam", "r")
    # df["output"] = "/Users/kezcleal/Documents/Data/fusion_finder_development/Split_read_simulator/sv_gen_mapped/test.fnfi.unsrt.sam"

    #input_stream_alignments.process_reads(df)
    k = defaults
    k["sv_aligns"] = "/Users/kezcleal/Documents/Data/fusion_finder_development/kates_benchmarked_data/fnfi2_out/output/DB120.hq_all.fnfi.srt.bam"
    # k["raw_aligns"] = "/Users/kezcleal/Documents/Data/fusion_finder_development/Event_simulator/Events/bwa.0.2.srt.bam"
    k["include"] = "/Users/kezcleal/Documents/Data/fusion_finder_development/test/include_tels.bed"
    k["svs_out"] = "/Users/kezcleal/Documents/Data/fusion_finder_development/kates_benchmarked_data/fnfi2_out/DB120.test.csv"
    k["procs"] = 1
    k["I"] = "142,187"
    #name = "fufi2_id{}".format("0.2")

    cluster.cluster_reads(k)
    #call("diff /Users/kezcleal/Documents/Data/fusion_finder_development/Split_read_simulator/paired_end_with_one_split/pairs.fnfi.unsrt.sam /Users/kezcleal/Documents/Data/fusion_finder_development/Split_read_simulator/sv_gen_mapped/test.fnfi.unsrt.sam", shell=True)
