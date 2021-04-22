import parasail
import sys
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref', help='the reference file', required=False)
parser.add_argument('--motif', help='the motif', required=False)
parser.add_argument('--start', help='the start', required=False)
parser.add_argument('--end', help='the end', required=False)
parser.add_argument('--chr', help='the chromosome', required=False)
parser.add_argument('--name', help='the name of the motif', required=False)

args = parser.parse_args()

control_file = args.ref

control_file = args.ref
motif = args.motif
start = args.start
end = args.end
chromosome = args.chr
name_motif = args.name

print("chr\tbegin\tend\tname\trepeat\tprefix\tsuffix")

for read in pysam.FastxFile(control_file):
    if read.name == chromosome:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chromosome,begin,end,name_motif,motif,read.sequence[int(begin)-101:int(begin)-1],read.sequence[int(end)+1:int(end)+101]))