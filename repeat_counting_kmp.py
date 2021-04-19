import sys
import pysam
import argparse

# Python program for KMP Algorithm 
def KMPSearch(pat, txt): 
	M = len(pat) 
	N = len(txt) 

	locations = list()
	# create lps[] that will hold the longest prefix suffix 
	# values for pattern 
	lps = [0]*M 
	j = 0 # index for pat[] 

	# Preprocess the pattern (calculate lps[] array) 
	computeLPSArray(pat, M, lps) 

	i = 0 # index for txt[] 
	while i < N: 
		if pat[j] == txt[i]: 
			i += 1
			j += 1

		if j == M: 
			#print ("Found pattern at index " + str(i-j)) 
			locations.append((i-j))
			j = lps[j-1] 

		# mismatch after j matches 
		elif i < N and pat[j] != txt[i]: 
			# Do not match lps[0..lps[j-1]] characters, 
			# they will match anyway 
			if j != 0: 
				j = lps[j-1] 
			else: 
				i += 1
	return locations

def computeLPSArray(pat, M, lps): 
	len = 0 # length of the previous longest prefix suffix 

	lps[0] # lps[0] is always 0 
	i = 1

	# the loop calculates lps[i] for i = 1 to M-1 
	while i < M: 
		if pat[i]== pat[len]: 
			len += 1
			lps[i] = len
			i += 1
		else: 
			# This is tricky. Consider the example. 
			# AAACAAAA and i = 7. The idea is similar 
			# to search step. 
			if len != 0: 
				len = lps[len-1] 

				# Also, note that we do not increment i here 
			else: 
				lps[i] = 0
				i += 1

parser = argparse.ArgumentParser()
parser.add_argument('--ref', help='the reference file', required=True)
parser.add_argument('--output', help='the output file', required=True)
parser.add_argument('--motif', help='the motif to be searched', required=True)

args = parser.parse_args()

control_file = args.ref
output_file = args.output
pat = args.motif
output_fd = open(output_file,'w')

output_fd.write("Chromosome\tStart\tEnd\tCount\n")
for read in pysam.FastxFile(control_file):
	print("Chromosome: %s"%(read.name))
	seq = read.sequence
	locations = KMPSearch(pat, seq) 
	
	init = 0
	count = 1
	start = 0
	end = 0
	for pos in locations:
		if pos - init == len(pat):
			if count == 1:
				start = init
			count = count + 1
			init = pos
			end = pos + len(pat)
		else:
			#end = init + len(pat)
			if count >= 10: 
				output_fd.write("%s\t%d\t%d\t%d\n"%(read.name,start,end,count))
				print("Found %d %s consecutively from position %d till %d in %s"%(count,pat,start,end,read.name))
			count = 1 
			init = pos
			start = pos
	

