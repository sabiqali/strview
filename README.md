# Introduction

A bunch of tools to help process and count STRs.

# strview
A tool to elaborate upon alignments between the read and reference around the STR region <br />

usage: ` strview.py [-h] [--bam BAM] [--read READ] --ref REF --config CONFIG --output OUTPUT [--parasail PARASAIL] [--pysam PYSAM] [--verbose VERBOSE] ` <br />

optional arguments: <br />
  -h, --help           show this help message and exit <br />
  --bam BAM            the bam file for the alignments between the reads and the reference genome(Required)<br />
  --read READ          the read file<br />
  --ref REF            the reference genome(Required)<br />
  --config CONFIG      the config file which contains information regarding the <br />
  --output OUTPUT      the output file where the output will be written in .tsv format<br />
  --parasail PARASAIL  Obtain the alignments using parasail. Must be used with the verbose flag<br />
  --pysam PYSAM        Obtain the alignments from the BAM File. Must be used with the verbose flag<br />
  --verbose VERBOSE    display the alignments of the different regions when enabled but will default to displaying the value in a tsv format <br />

  The config file should be in the following format:<br />

  | chr | begin | end | name | repeat | prefix | suffix | 
  | --- | --- | --- | --- | --- | --- | --- | 
  | chr9 | 27573527 | 27573544 | c9orf72 | GGCCCC | CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC | TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC | 

# Repeat Counting 
A tool to find and count the number of times a motif is located in the provided genome. Only the motif and the reference genome is required. 

usage: `repeat_counting_kmp.py [-h] --ref REF --output OUTPUT --motif MOTIF ` <br />

optional arguments: <br />
  -h, --help       show this help message and exit <br />
  --ref REF        the reference genome <br />
  --output OUTPUT  the output file <br />
  --motif MOTIF    the motif to be searched <br />


