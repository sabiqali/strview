import matplotlib.pyplot as plt
import pandas as pd
import argparse
import seaborn as sns

import statistics
from statistics import mode

import os

#def get_file_name(file_path):
#    return dir.rsplit('.', 1)[0].split('/')[-1]
    #return os.path.basename(file_path).split('.')[0]
 
def most_common(List):
    return(mode(List))

parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from strview', required=True)
#parser.add_argument('--name', help='the graph name', required=True)

args = parser.parse_args()

input_file = open(args.input)
name = os.path.basename(args.input).split('.')[0]

output_file_name = name+".png"

outputs = list()
for line in input_file:
    outputs.append(line.rstrip().split()[0:21])

count_list=[]
for line in outputs:
    c = 0
    character = line[11].split(">")
    if(len(character) == 1):
        character = line[11].split("<")
    for word in character:
        if (word == "repeat"):
            c = c + 1
    #if(c != 0):
    count_list.append(c)
    print("\t".join([line[0],str(c),line[11]]))
    #print(character)

count_list.sort(reverse=True)
#n, bins, patches = plt.hist(count_list)
n, bins, patches = plt.hist(count_list, bins='auto')
plt.xlabel('count')
plt.ylabel('number of reads')
plt.title('distribution')
#plt.show()
plt.savefig(output_file_name)

#print(count_list)
#print("Mean %d"%(statistics.mean(count_list)))
#print("Median %d"%(statistics.median(count_list)))
#print("Mode %d"%(most_common(count_list)))
#print("SD %d"%(statistics.stdev(count_list)))
