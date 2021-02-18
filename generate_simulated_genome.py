import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--repeat-motif', help='the repeat motif, which will have the repeat expansion', dest="repeat", default="ATG",required=False)
parser.add_argument('--repeat-count', help='the number of repeat counts to be expanded in the genome', dest="count", default=100,required=False)
parser.add_argument('--flanking-sequence-length', help='the length of the flanking sequence on either side of the repeat', dest="length",required=True)
parser.add_argument('--output-config', help='the STRique config file that will be generated', dest="config",required=True)

args = parser.parse_args()

repeat = args.repeat
count = int(args.count)
length = int(args.length)
config_file_name = args.config

config_file = open(config_file_name,'w')

def DNA(length_of_flanking_sequence):
    return ''.join(random.choice('CGTA') for _ in range(length_of_flanking_sequence))

def repeat_expansion(repeat_motif, number_of_repeats):
    return repeat_motif * number_of_repeats

config_file.write("chr\tbegin\tend\tname\trepeat\tprefix\tsuffix\n")

left_flank = DNA(length)
right_flank = DNA(length)
RE = repeat_expansion(repeat,count)
sequence_name = "Synthetic_Sequence"

config_file.write("%s\t%d\t%d\t%s\t%s\t%s\t%s"%(sequence_name,length,(length+150),"synthetic",repeat,left_flank[length-150:],right_flank[:150]))

print(">%s"%(sequence_name))
print (left_flank+RE+right_flank)
