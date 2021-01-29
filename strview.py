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
    print("read_name\tchromosome\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence")
    align_data_file.write("read_name\tchromosome\tposition\tprefix_sequence\trepeat_sequence\tsuffix_sequence\n")

    for alignment in bamfile:
        chr_rname = []
        if alignment.rname == 0:
            chr_rname = "chr1"
        if alignment.rname == 1:
            chr_rname = "chr2"
        if alignment.rname == 2:
            chr_rname = "chr3"
        if alignment.rname == 3:
            chr_rname = "chr4"
        if alignment.rname == 4:
            chr_rname = "chr5"
        if alignment.rname == 5:
            chr_rname = "chr6"
        if alignment.rname == 6:
            chr_rname = "chr7"
        if alignment.rname == 7:
            chr_rname = "chr8"
        if alignment.rname == 8:
            chr_rname = "chr9"
        if alignment.rname == 9:
            chr_rname = "chr10"
        if alignment.rname == 10:
            chr_rname = "chr11"
        if alignment.rname == 11:
            chr_rname = "chr12"
        if alignment.rname == 12:
            chr_rname = "chr13"
        if alignment.rname == 13:
            chr_rname = "chr14"
        if alignment.rname == 14:
            chr_rname = "chr15"
        if alignment.rname == 15:
            chr_rname = "chr16"
        if alignment.rname == 16:
            chr_rname = "chr17"
        if alignment.rname == 17:
            chr_rname = "chr18"
        if alignment.rname == 18:
            chr_rname = "chr19"
        if alignment.rname == 19:
            chr_rname = "chr20"
        if alignment.rname == 20:
            chr_rname = "chr21"
        if alignment.rname == 21:
            chr_rname = "chr22"
        if alignment.rname == 22:
            chr_rname = "chrX"
        print("%s\t%s\t%s\t%s\t%s\t%s\n" % (alignment.qname,chr_rname,alignment.pos,prefix,repeat,suffix))
        print(alignment.rname)
        align_data_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (alignment.qname,chr_rname,alignment.pos,prefix,repeat,suffix))
        break

#TODO::CONVERT SAM CIGAR TO ALIGNMENT


#PARASAIL ALIGNMENT SEGMENT
if args.verbose == 1:
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

align_data_file.close()
#print(read_seq)
