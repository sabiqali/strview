#! /usr/env/python

import parasail
import sys
import pysam
import argparse
import math

class ReadAlignment:
    def __init__(self, name):
        self.read_name: name
        self.has_prefix_match = False
        self.has_suffix_match = False
    
    count = 0

def percentage_identity(cigar_exp):
    m = 0
    nm = 0
    for chr in cigar_exp:
        if chr == '|':
            m = m + 1
        else:
            nm = nm + 1
    pi = float(m/(nm + m))
    return pi

def roundup(x):
    return int(math.ceil(x / 100000.0)) * 100000

def alignment_contains_str_flank(alignment, flank, matrix):
    result = parasail.sw_trace_scan_32(flank, alignment.query_sequence, 5, 4, scoring_matrix)
    if(percentage_identity(result.traceback.comp) > 0.6):
        return True
    else: 
        return False

parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the bam file', required=False)
parser.add_argument('--read', help='the read file', required=False)
parser.add_argument('--ref', help='the ref file', required=True)
parser.add_argument('--config', help='the config file', required=True)
parser.add_argument('--output', help='the config file', required=False)
parser.add_argument('--verbose', help='display the alignments of the different regions', type=int, required=False, default=0)

args = parser.parse_args()

reads_file = args.read
reference_file = args.ref
config = args.config
in_bam = args.bam

# read in the configs and store them in the respective variables
configs = list()
config_fh = open(config)
header = config_fh.readline()

for line in config_fh:
    configs.append(line.rstrip().split()[0:7])

chromosome,begin,end,name,repeat,prefix,suffix = configs[0]
    
#First pass through the alignment file to determine what are the matches.
bamfile = pysam.AlignmentFile(in_bam)
upper_limit = roundup(int(begin))
lower_limit = upper_limit - 100000
idx = 0
scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
reads = dict()
for alignment in bamfile.fetch(chromosome,lower_limit,upper_limit):
    if alignment.qname not in reads:
        reads[alignment.qname] = ReadAlignment(alignment.qname)
    if alignment_contains_str_flank( alignment , prefix, scoring_matrix):
        reads[alignment.qname].has_prefix_match = True
    if alignment_contains_str_flank( alignment, suffix, scoring_matrix):
        reads[alignment.qname].has_suffix_match = True

#Once we have all the matches, we can iterate through them to get the count
print("\t".join(["read_name","chromosome","repeat_name","count","aligned_query","aligned_ref"]))
for read_name , alignment in reads.items():
    if not reads[read_name].has_prefix_match or not reads[read_name].has_suffix_match:
        continue
    with pysam.FastaFile(reads_file) as fh:
        entry = fh.fetch(read_name)
        read_seq = entry
        prev_score = 0 
        ideal_read = prefix + repeat + suffix
        result = parasail.sw_trace_scan_32(read_seq, ideal_read, 5, 4, scoring_matrix)
        prev_result_ref = result.traceback.ref
        result_ref = result.traceback.ref
        result_comp = result.traceback.comp
        prev_result_query = result.traceback.query
        result_query = result.traceback.query
        score = result.score
        c = 1
        while (score > prev_score) and (percentage_identity(result_comp) > 0.5):
            c = c + 1
            prev_score = score
            prev_result_ref = result_ref
            prev_result_query = result_query
            ideal_read = prefix + ( repeat * c ) + suffix
            result = parasail.sw_trace_scan_32(read_seq, ideal_read, 5, 4, scoring_matrix)
            score = result.score
            result_comp = result.traceback.comp
            result_ref = result.traceback.ref
            result_query = result.traceback.query
        max_score = prev_score
        reads[read_name].count = c - 1
    
    print( "\t".join([read_name,chromosome,name,str(reads[read_name].count),prev_result_query,prev_result_ref]))
