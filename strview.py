#! /usr/env/python

import parasail
import sys
import pysam
import argparse

#parser = argparse.ArgumentParser(prog='strview', usage= 'python %(strview)s [options]',description='A small tool to elaborate upon alignments between the read and reference around the STR region',)
#parser.add_argument('--bam',dest='bam_file', help='the bam file that contains the alignment of the read to the reference')
#parser.add_argument('--ref',dest='ref_file', help='the reference file')
#parser.add_argument('--reads',dest='reads_file' help='the reads file')
#parser.add_argument('--config',dest='config_file', help='the STR config file that contains the repeat configuration with the coordinates, motif, prefix, and suffix')
#args = parser.parse_args()
parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the bam file')
parser.add_argument('--read', help='the read file')
parser.add_argument('--ref', help='the ref file')
parser.add_argument('--config', help='the config file')

args = parser.parse_args()

reads_file = args.read
control_file = args.ref
config = args.config

# read in the control sequences
configs = list()
config_fh = open(config)
header = config_fh.readline()

for line in config_fh:
    configs.append(line.rstrip().split()[0:7])

for (chromosome,begin,end,name,repeat,prefix,suffix) in configs:
    print(chromosome)
    break

scoring_matrix = parasail.matrix_create("ACGT", 5, -1)

chr9_seq = []
local_seq = []
read_seq = []

align_data_file = open('alignment_data.txt','w')

for read in pysam.FastxFile(control_file):
    print(read.name)
    if(read.name == chromosome):
        chr9_seq = read.sequence
        #print(len(chr9_seq))
        local_seq = chr9_seq[27572000 : 27575000]
        #print(local_seq)
        #print(read.sequence)

for read in pysam.FastxFile(reads_file):
    read_seq = read.sequence
    break

#prefix = "CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC"

result = parasail.sw_trace_scan_32(local_seq, prefix, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

result = parasail.sw_trace_scan_32(read_seq, prefix, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

align_data_file.write("Read Name: %s \n" % (read.name))
align_data_file.write("Prefix Score: %d \n" % (result.score))
align_data_file.write(result.traceback.query)
align_data_file.write("\n")
align_data_file.write(result.traceback.comp)
align_data_file.write("\n")
align_data_file.write(result.traceback.ref)
align_data_file.write("\n") 

#suffix = "TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC"

result = parasail.sw_trace_scan_32(local_seq, suffix, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

result = parasail.sw_trace_scan_32(read_seq, suffix, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

align_data_file.write("Read Name: %s \n" % (read.name))
align_data_file.write("Suffix Score: %d \n" % (result.score))
align_data_file.write(result.traceback.query)
align_data_file.write("\n")
align_data_file.write(result.traceback.comp)
align_data_file.write("\n")
align_data_file.write(result.traceback.ref)
align_data_file.write("\n") 

#repeat = "GGCCCC"

result = parasail.sw_trace_scan_32(local_seq, repeat, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

result = parasail.sw_trace_scan_32(read_seq, repeat, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

total_analysis = prefix+repeat+repeat+repeat+suffix

result = parasail.sw_trace_scan_32(local_seq, total_analysis, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)

result = parasail.sw_trace_scan_32(read_seq, total_analysis, 5, 4, scoring_matrix)

print("Score: %d " % (result.score))
print(read.name)
print(result.traceback.query)
print(result.traceback.comp)
print(result.traceback.ref)
align_data_file.write("Read Name: %s \n" % (read.name))
align_data_file.write("Score: %d \n" % (result.score))
align_data_file.write(result.traceback.query)
align_data_file.write("\n")
align_data_file.write(result.traceback.comp)
align_data_file.write("\n")
align_data_file.write(result.traceback.ref)
align_data_file.write("\n") 

align_data_file.close()
#print(read_seq)
