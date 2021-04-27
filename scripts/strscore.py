#! /usr/env/python

import parasail
import sys
import pysam
import argparse
import math

def roundup(x):
    return int(math.ceil(x / 100000.0)) * 100000

parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the bam file', required=False)
parser.add_argument('--read', help='the read file', required=False)
parser.add_argument('--ref', help='the ref file', required=True)
parser.add_argument('--config', help='the config file', required=True)
parser.add_argument('--output', help='the config file', required=True)
parser.add_argument('--verbose', help='display the alignments of the different regions', type=int, required=False, default=0)

args = parser.parse_args()

reads_file = args.read
control_file = args.ref
config = args.config
in_bam = args.bam

align_data_file = open(args.output,'w')

# read in the control sequences
configs = list()
config_fh = open(config)
header = config_fh.readline()

for line in config_fh:
    configs.append(line.rstrip().split()[0:7])

for (chromosome,begin,end,name,repeat,prefix,suffix) in configs:
    print(chromosome)
    break

ref_seq = []
local_seq = []
for read in pysam.FastxFile(control_file):
    print(read.name)
    if(read.name == chromosome):
        ref_seq = read.sequence

bamfile = pysam.AlignmentFile(in_bam)
if args.verbose == 0:
    print("read_name\tchromosome\tcount\tposition\taligned_query\taligned_ref")
    align_data_file.write("read_name\tchromosome\tcount\tposition\taligned_query\taligned_ref\n")

upper_limit = roundup(int(begin))
lower_limit = upper_limit - 100000    
idx = 0
scoring_matrix = parasail.matrix_create("ACGT", 5, -1)

with pysam.FastaFile(reads_file) as fh:
    for entry in fh.fetch(chromosome,lower_limit,upper_limit):
        print(entry.name)

        read_seq = entry.sequence
        prev_score = 0 
        ideal_read = prefix + repeat + suffix
        result = parasail.sw_trace_scan_32(read_seq, ideal_read, 5, 4, scoring_matrix)
        prev_result_ref = result.traceback.ref
        result_ref = result.traceback.ref
        prev_result_query = result.traceback.query
        result_query = result.traceback.query
        score = result.score
        c = 1
        while score > prev_score:
            c = c + 1
            prev_score = score
            prev_result_ref = result_ref
            prev_result_query = result_query
            ref_seq = prefix + ( repeat * c ) + suffix
            result = parasail.sw_trace_scan_32(read_seq, ideal_read, 5, 4, scoring_matrix)
            score = result.score
            result_ref = result.traceback.ref
            result_query = result.traceback.query
        max_score = prev_score
        count = c - 1

        print("%s\t%s\t%d\t%s\t%s\t%s\n" % (alignment.qname,chromosome,tmp_count,alignment.pos,prev_result_query,prev_result_ref))
                
        align_data_file.write("%s\t%s\t%d\t%s\t%s\t%s\n" % (alignment.qname,chromosome,tmp_count,alignment.pos,prev_result_query,prev_result_ref))

