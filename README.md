# strview
A small tool to elaborate upon alignments between the read and reference around the STR region

usage: strview.py [-h] [--bam BAM] [--read READ] --ref REF --config CONFIG --output OUTPUT [--parasail PARASAIL] [--pysam PYSAM] [--verbose VERBOSE]

optional arguments:
  -h, --help           show this help message and exit
  --bam BAM            the bam file
  --read READ          the read file
  --ref REF            the ref file
  --config CONFIG      the config file
  --output OUTPUT      the config file
  --parasail PARASAIL  Obtain the alignments using parasail. Must be used with the verbose flag
  --pysam PYSAM        Obtain the alignments from the BAM File. Must be used with the verbose flag
  --verbose VERBOSE    display the alignments of the different regions
