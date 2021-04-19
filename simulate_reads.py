import random
import argparse
import pysam
import random


parser = argparse.ArgumentParser()
parser.add_argument('--indel-error-probability', help='the repeat motif, which will have the repeat expansion', dest="prob", default="ATG",required=True)
parser.add_argument('--num-reads', help='the number of repeat counts to be expanded in the genome', dest="num", default=100,required=False)
parser.add_argument('--genome', help='the length of the flanking sequence on either side of the repeat', dest="ref",required=True)
parser.add_argument('--output-reads', help='the STRique config file that will be generated', dest="output",required=False)

args = parser.parse_args()

probability = float(args.prob)
num_reads = int(args.num)
genome_file = args.ref

for read in pysam.FastxFile(genome_file):
    ref_seq = read.sequence

random.seed(1)

for i in range(0,num_reads):
    read_sequence = ""
    for base in ref_seq:
        if random.random() < probability:  # an error has occurred at this base
            if random.random() < 0.5:  # the error is an insertion
                inserted_base = random.choice('CGTA')
                read_sequence += base + inserted_base
            #else: # the error is a deletion, do nothing
        else: # no error, append base
            read_sequence += base
    print(">Simulated_Read_%d"%(i))
    print("%s"%(read_sequence))
    