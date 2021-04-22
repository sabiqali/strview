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
parser.add_argument('--score', help='Have a look at the read, with an ideal reference, to see what count would maximise the alignment score', required=False, type=int, default=0)
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
    if args.score == 1:
        print("read_name\tchromosome\tcount\talt_count\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence")
        align_data_file.write("read_name\tchromosome\tcount\talt_count\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence\n")
    else:
        print("read_name\tchromosome\tcount\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence")
        align_data_file.write("read_name\tchromosome\tcount\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence\n")

upper_limit = roundup(int(begin))
lower_limit = upper_limit - 100000    
idx = 0
scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
if args.pysam == 1:
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
        pair_out_identity = alignment.get_aligned_pairs()

        pair_out = alignment.get_aligned_pairs(True)
        
        tmp = alignment.query_sequence
        read_seq = tmp

        tmp_num = 0
        tmp_repeat_low = 0
        tmp_repeat_upper = 0
        for tmp_pairs in  pair_out:
            if tmp_pairs[1] == int(begin) - 1:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 1
            if tmp_pairs[1] == int(begin) - 2:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 2
            if tmp_pairs[1] == int(begin) - 3:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 3
            if tmp_pairs[1] == int(begin) - 4:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 4
            if tmp_pairs[1] == int(begin) - 5:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 5
            if tmp_pairs[1] == int(begin) - 6:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 6
            if tmp_pairs[1] == int(begin) - 7:
                tmp_num = tmp_pairs[0]
                tmp_repeat_low = tmp_num + 7

            if tmp_pairs[1] == int(end) + 1:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 1
            if tmp_pairs[1] == int(end) + 2:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 2
            if tmp_pairs[1] == int(end) + 3:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 3
            if tmp_pairs[1] == int(end) + 4:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 4
            if tmp_pairs[1] == int(end) + 5:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 5
            if tmp_pairs[1] == int(end) + 6:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 6
            if tmp_pairs[1] == int(end) + 7:
                tmp_num = tmp_pairs[0]
                tmp_repeat_upper = tmp_num - 7
        
        idx = idx + 1
        #print(idx)
        tmp_count = int(len(read_seq[tmp_repeat_low:tmp_repeat_upper]) / len(repeat))
        aligned_prefix = read_seq[(tmp_repeat_low - 100) if tmp_repeat_low - 100 > 0 else 0 : tmp_repeat_low]
        aligned_repeat = read_seq[tmp_repeat_low:tmp_repeat_upper]
        aligned_suffix = read_seq[tmp_repeat_upper: (tmp_repeat_upper + 100) if tmp_repeat_upper + 100 < len(read_seq) else len(read_seq)]

        #for tmp_pair in pair_out:                                                         #reconstruct the aligned segments from the cigar
        #    if tmp_pair[1] > (int(begin) - len(prefix)) and tmp_pair[1] < int(begin):
        #        aligned_prefix = aligned_prefix + read_seq[tmp_pair[0]]
        #        aligned_ref_prefix = aligned_ref_prefix + ref_seq[tmp_pair[1]]
        #        prefix_cigar = prefix_cigar + "|"
        #    if tmp_pair[1] > int(end) and tmp_pair[1] < (int(end) + len(suffix)):
        #        aligned_suffix = aligned_suffix + read_seq[tmp_pair[0]]
        #        aligned_ref_suffix = aligned_ref_suffix + ref_seq[tmp_pair[1]]
        #        suffix_cigar = suffix_cigar + "|"
        #    if tmp_pair[1] > int(begin) and tmp_pair[1] < int(end):
        #        aligned_repeat = aligned_repeat + read_seq[tmp_pair[0]]
        #        aligned_ref_repeat = aligned_ref_repeat + ref_seq[tmp_pair[1]]
        #        repeat_cigar = repeat_cigar + "|"
        #        count = aligned_repeat.count(repeat)
        if args.score == 1:
            prev_score = 0 
            ref_seq = prefix + repeat + suffix
            result = parasail.sw_trace_scan_32(read_seq, ref_seq, 5, 4, scoring_matrix)
            score = result.score
            c = 1
            while score > prev_score:
                c = c + 1
                prev_score = score
                ref_seq = prefix + ( repeat * c ) + suffix
                result = parasail.sw_trace_scan_32(read_seq, ref_seq, 5, 4, scoring_matrix)
                score = result.score
            max_score = prev_score
            alt_count = c - 1

        if aligned_prefix and aligned_repeat and aligned_suffix and args.verbose == 0 and tmp_count != 0:
            if args.score == 1:
                print("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,tmp_count,alt_count,alignment.pos,aligned_prefix,aligned_repeat,aligned_suffix))
                #print(alignment.rname)
                align_data_file.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,tmp_count,alt_count,alignment.pos,aligned_prefix,aligned_repeat,aligned_suffix))
            else:
                print("%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,tmp_count,alignment.pos,aligned_prefix,aligned_repeat,aligned_suffix))
                #print(alignment.rname)
                align_data_file.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,tmp_count,alignment.pos,aligned_prefix,aligned_repeat,aligned_suffix))
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


#idx = 0

#if args.score == 1:
#    scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
#    for alignment in bamfile.fetch(chromosome,lower_limit,upper_limit):
#        read_seq = alignment.query_sequence
#        prev_score = 0 
#        ref_seq = prefix + repeat + suffix
#        result = parasail.sw_trace_scan_32(read_seq, ref_seq, 5, 4, scoring_matrix)
#        score = result.score
#        c = 1
#        while score > prev_score:
#            c = c + 1
#            prev_score = score
#            ref_seq = prefix + ( repeat * c ) + suffix
#            result = parasail.sw_trace_scan_32(read_seq, ref_seq, 5, 4, scoring_matrix)
#            score = result.score
#        max_score = prev_score
#        count = c - 1
#        print("%s\t%s\t%d\t%s\n" % (alignment.qname,chromosome,count,alignment.pos))

idx = 0

#PARASAIL ALIGNMENT SEGMENT
if args.parasail == 1 :
    #scoring_matrix = parasail.matrix_create("ACGT", 5, -1)

    align_data_file = open('alignment_data.txt','w')

    for alignment in bamfile.fetch(chromosome,lower_limit,upper_limit):
        read_seq = alignment.query_sequence

        result_pre_traceback = parasail.sw_trace_scan_32(read_seq, prefix, 5, 4, scoring_matrix)
        result_pre = parasail.ssw(read_seq, prefix, 5, 4, scoring_matrix)

        #print("Score: %d " % (result_pre_traceback.score))
        #print(alignment.qname)
        #print(result_pre.ref_begin1)
        #print(result_pre.ref_end1)
        #print(percentage_identity(result_pre_traceback.traceback.comp))

        result_suf_traceback = parasail.sw_trace_scan_32(read_seq, suffix, 5, 4, scoring_matrix)
        result_suf = parasail.ssw(read_seq, suffix, 5, 4, scoring_matrix)

        #print("Score: %d " % (result_suf_traceback.score))
        #print(alignment.qname)
        #print(result_suf.ref_begin1)
        #print(result_suf.ref_end1)
        #print(percentage_identity(result_suf_traceback.traceback.comp))

        if result_pre.ref_end1 < result_suf.ref_begin1:
            count = len(read_seq[result_pre.ref_end1+1:result_suf.ref_begin1-1]) / len(repeat)
            aligned_repeat = read_seq[result_pre.ref_end1+1:result_suf.ref_begin1-1]
        else: 
            count = 0
            aligned_repeat = ""

        print("%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (alignment.qname,chromosome,count,alignment.pos,result_pre_traceback.traceback.query,aligned_repeat,result_suf_traceback.traceback.query))

        #print("\n\nProcessed %d bam reads"%(idx))
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
