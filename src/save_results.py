

"""
Save results in folder
"""

import os
from shutil import copyfile
import glob
from subprocess import call
import pysam


def to_folder(args):

    if args.default_outpath == "":
        # Use current directory
        args.default_outpath = os.getcwd()

    out_folder = args.default_outpath
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    print("OUTPUT=", out_folder)

    # if args.merge:
    #     call("samtools merge {d}/{v}.splitmap.bam \
    #     {t}/{s}.na_sr.bam \
    #     {t}/{s}.pairs.sc.srt.fixed.srt.bam".format(d=out_folder, v=args.samp[:-4], s=args.samp, t=args.tmp), shell=True)
    #     call("samtools index {d}/{s}.splitmap.bam".format(d=out_folder, s=args.samp[:-4]), shell=True)
    # else:
    to_move = [#(".na_sr.bam", ".na_sr.bam"),
               (".pairs.sc.srt.fixed.srt.bam", ".sr.bam"),
               (".pairs.sc.srt.fixed.srt.bam.bai", ".sr.bam.bai")]

    for old, new in to_move:
        copyfile("{t}/{s}".format(t=args.tmp, s=args.samp) + old,
                 "{d}/{s}".format(d=out_folder, s=args.samp) + new)
    call("samtools index {d}/{s}.sr.bam".format(d=out_folder, s=args.samp), shell=True)
    print ("samtools index {d}/{s}.sr.bam".format(d=out_folder, s=args.samp))
    # cleanup
    if not args.debug:
        fs = glob.glob("{t}/{s}.pairs.*".format(t=args.tmp, s=args.samp)) + \
             ["{t}/{s}.sizes".format(t=args.tmp, s=args.samp)]
        for item in fs:
            print("Deleting temp file", item)
            #os.remove(item)


def splitters_only(args):
    final_outname = os.path.dirname(args.bam) + "/" + os.path.splitext(args.samp)[0] + ".splitmap.splitters.bam"
    b = os.path.dirname(args.bam) + "/" + os.path.splitext(args.samp)[0] + ".splitmap.bam"

    inbam = pysam.AlignmentFile(b, "rb")
    outbam = pysam.AlignmentFile(final_outname, "wb", template=inbam)

    count = 0
    for r in inbam:
        #if r.flag & 2048:
        if any((i == "SP" and j == "1") for i, j in r.get_tags()):
            outbam.write(r)
            count += 1

    print("Total splitters: ", count)

    inbam.close()
    outbam.close()

    call("samtools index " + final_outname, shell=True)

