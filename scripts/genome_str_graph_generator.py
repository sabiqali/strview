import sys
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref', help='the ref file', required=True)
parser.add_argument('--config', help='the config file', required=True)

args = parser.parse_args()

reference_file = args.ref
config = args.config

# read in the configs and store them in the respective variables
configs = list()
config_fh = open(config)
header = config_fh.readline()

for line in config_fh:
    configs.append(line.rstrip().split()[0:7])

#In this case I am extracting only 1 config but in case we have more than 1, we can have a loop here to extract the other configs
chromosome,begin,end,name,repeat,prefix,suffix = configs[0]

print("H\tVN:Z:a.0")

#appraoch it chromosome wise and then have a variable for sides and one for links for the chromosome and then in the end iterate over all and print. 

sides = dict()
links = dict()

for chr in pysam.FastxFile(reference_file):
    if(chr.name == chromosome):
        print("\t".join(["S", chr.name + "_prefix", chr.sequence[:int(begin)-1]]))
	print("\t".join(["S", "repeat" , chr.sequence[int(begin):int(end)]]))
	print("\t".join(["S", chr.name + "_suffix", chr.sequence[int(end)+1:]]))
	#links 
