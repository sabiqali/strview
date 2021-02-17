#! /usr/env/python

import parasail
import sys
import pysam
import argparse
import math

def find_read(read_file_name, read_from_cigar):
    read_seq = []
    for query_read in pysam.FastxFile(read_file_name):
        #print("%s\t%s"%(query_read.name,read.qname))
        if query_read.name == read_from_cigar: 
            read_seq = query_read.sequence
            #print(read_from_cigar)
            #print(query_read.name)
            return (read_seq,query_read.name)
    return ('not found','not found')

#def round_up(n, decimals=0):
#    multiplier = 10 ** decimals
#    return int(math.ceil(n * multiplier) / multiplier)

def roundup(x):
    return int(math.ceil(x / 100000.0)) * 100000

#parser = argparse.ArgumentParser(prog='strview', usage= 'python %(strview)s [options]',description='A small tool to elaborate upon alignments between the read and reference around the STR region',)
#parser.add_argument('--bam',dest='bam_file', help='the bam file that contains the alignment of the read to the reference')
#parser.add_argument('--ref',dest='ref_file', help='the reference file')
#parser.add_argument('--reads',dest='reads_file' help='the reads file')
#parser.add_argument('--config',dest='config_file', help='the STR config file that contains the repeat configuration with the coordinates, motif, prefix, and suffix')
#args = parser.parse_args()
parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the bam file', required=False)
parser.add_argument('--read', help='the read file', required=False)
parser.add_argument('--ref', help='the ref file', required=True)
parser.add_argument('--config', help='the config file', required=True)
parser.add_argument('--output', help='the config file', required=True)
parser.add_argument('--parasail', help='Obtain the alignments using parasail. Must be used with the verbose flag', required=False, type=int, default=0)
parser.add_argument('--pysam', help='Obtain the alignments from the BAM File. Must be used with the verbose flag', required=False, type=int, default=0)
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
        #print(len(ref_seq))
        #local_seq = ref_seq[27572000 : 27575000]
        #print(local_seq)
        #print(read.sequence)

#bam_file = pysam.AlignmentFile(in_bam)

#for read in bamfile:
#    pair_out = read.get_aligned_pairs(True)
    #for (x,y) in pair_out:
    #    if type(y) == int and type(x) == int :
    #        print(x,y)
    #        final_pair_out.extend(x,y)
#    print(pair_out)
#    break

bamfile = pysam.AlignmentFile(in_bam)
if args.verbose == 0:
    print("read_name\tchromosome\tcount\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence")
    align_data_file.write("read_name\tchromosome\tcount\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence\n")

upper_limit = roundup(int(begin))
lower_limit = upper_limit - 100000    
for alignment in bamfile.fetch(chromosome,lower_limit,upper_limit):
    aligned_prefix = ""
    aligned_ref_prefix = ""
    prefix_cigar = ""
    aligned_ref_suffix = ""
    aligned_suffix = ""
    suffix_cigar = ""
    aligned_ref_repeat = ""
    aligned_repeat = ""
    repeat_cigar = ""
    count = 0 
    idx = 0
    pair_out_identity = alignment.get_aligned_pairs()

    pair_out = alignment.get_aligned_pairs(True)
    
    tmp = alignment.query_sequence
    read_seq = tmp
    
    for tmp_pair in pair_out:                                                         #reconstruct the aligned segments from the cigar
        if tmp_pair[1] > (int(begin) - len(prefix)) and tmp_pair[1] < int(begin):
            aligned_prefix = aligned_prefix + read_seq[tmp_pair[0]]
            aligned_ref_prefix = aligned_ref_prefix + ref_seq[tmp_pair[1]]
            prefix_cigar = prefix_cigar + "|"
        if tmp_pair[1] > int(end) and tmp_pair[1] < (int(end) + len(suffix)):
            aligned_suffix = aligned_suffix + read_seq[tmp_pair[0]]
            aligned_ref_suffix = aligned_ref_suffix + ref_seq[tmp_pair[1]]
            suffix_cigar = suffix_cigar + "|"
        if tmp_pair[1] > int(begin) and tmp_pair[1] < int(end):
            aligned_repeat = aligned_repeat + read_seq[tmp_pair[0]]
            aligned_ref_repeat = aligned_ref_repeat + ref_seq[tmp_pair[1]]
            repeat_cigar = repeat_cigar + "|"
            count = aligned_repeat.count("GGCCCC")
    if aligned_prefix and aligned_repeat and aligned_suffix and args.verbose == 0:
        print("%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,count,alignment.pos,aligned_prefix,aligned_repeat,aligned_suffix))
        #print(alignment.rname)
        align_data_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,alignment.pos,aligned_prefix,aligned_repeat,aligned_suffix))
    if aligned_prefix and aligned_repeat and aligned_suffix and args.verbose == 1 and args.pysam == 1:
        print(aligned_ref_prefix)
        print(prefix_cigar)
        print(aligned_prefix)
        print(aligned_ref_repeat)
        print(repeat_cigar)
        print(aligned_repeat)
        print(aligned_ref_suffix)
        print(suffix_cigar)
        print(aligned_suffix)


#TODO::CONVERT SAM CIGAR TO ALIGNMENT

idx = 0

#PARASAIL ALIGNMENT SEGMENT
if args.verbose == 1 and args.parasail == 1 :
    scoring_matrix = parasail.matrix_create("ACGT", 5, -1)

    read_seq = []

    align_data_file = open('alignment_data.txt','w')

    for read in pysam.FastxFile(reads_file):
        read_seq = read.sequence
        #break

        #align prefix to read

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

        #align suffix to read

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

        #align repeat to read

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

        #align the final analysis to the read
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

        print("\n\nProcessed %d bam reads"%(idx))
        idx = idx + 1


#if args.verbose == 1 and args.pysam == 1 :
#    read_from_cigar = ""
#    reference_from_cigar = ""
#    expanded_cigar = ""
#    hit = 0
#    for read in bamfile:
#        pair_out = read.get_aligned_pairs(True)
#        tmp = find_read(reads_file,read.qname)
#        read_seq = tmp[0]
#        if tmp[1] != read.qname:
#            #print("%s\t%s"%(tmp[1],read.qname))
#            continue
#        for tmp_pair in pair_out:
#            #print(read_seq[tmp_pair[0]])
#            if tmp_pair[1] == 27573528:
#                hit = hit + 1
#            read_from_cigar = read_from_cigar + read_seq[tmp_pair[0]]
#            reference_from_cigar = reference_from_cigar + ref_seq[tmp_pair[1]]
#            if tmp_pair[0] and tmp_pair[1]:
#                expanded_cigar = expanded_cigar + '|'
#            else:
#                expanded_cigar = expanded_cigar + ' '
#        print(reference_from_cigar)
#        print(expanded_cigar)
#        print(read_from_cigar)
#        if hit != 0 :
#            print('FOUND IT!!!!!!!!!!!!!!!!!!!!!')
#            hit = 0
#        print("\n\nProcessed %d bam reads"%(idx))
#        idx = idx + 1

align_data_file.close()
#print(read_seq)
