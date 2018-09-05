from subprocess import call
import os

inputf = "/Volumes/Kez8T/breastCancer88/Analysis/structural_variants/splitmap_development/DB143.last_out_small.sam"

inputf = "/Volumes/Kez8T/breastCancer88/Analysis/structural_variants/splitmap_development/DB143.small1.bwa-mem-all.sam"

inputf2 = "/Volumes/Kez8T/breastCancer88/Analysis/structural_variants/splitmap_development/DB143.big.bayName.fastq"
outf = "/Volumes/Kez8T/breastCancer88/Analysis/structural_variants/splitmap_development/DB143.big.bwa-mem-all.splitmap.sam"
test_bam = os.path.dirname(__file__) + "/data/test.bam"


call("bwa mem -t8 -P -a -p ~/Documents/Data/db/hg38/hg38.fa {} | python __main__.py -k sam > {}".format(inputf2, outf), shell=True)

#call("time bwa mem -t8 -p ~/Documents/Data/db/hg38/hg38.fa {} > ~/test.sam".format(inputf2, outf), shell=True)
#call("less {} | python __main__.py -k sam > {}".format(inputf, outf), shell=True)

