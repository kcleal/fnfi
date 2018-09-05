fufi_

fufi is a structural variant detection tool that specialises
in identifying rearrangements involving telomere fusions.


Example Usage:
$ pip install fufi
$ fufi run -i your.bam -o output_dir
$ fufi locate-telomeres -i your_telomeres.fa; fufi run -i your.bam -r your_telomeres.bed -ooutput_dir
$ bma mem ref.fa your.fq | fufi align -o output_dir
