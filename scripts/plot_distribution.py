import matplotlib.pyplot as plt
import pandas as pd
import argparse
import seaborn as sns

#def most_frequent(List):
#    counter = 0
#    num = List[0]
#     
#    for i in List:
#        curr_frequency = List.count(i)
#        if(curr_frequency> counter):
#            counter = curr_frequency
#            num = i
# 
#    return num


import statistics
from statistics import mode
 
def most_common(List):
    return(mode(List))

parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from strview', required=True)
parser.add_argument('--name', help='the graph name', required=True)

args = parser.parse_args()

input_file = open(args.input)
name = args.name

output_file_name = name+".png"

header = input_file.readline()
outputs = list()
for line in input_file:
    outputs.append(line.rstrip().split()[0:6])

count_list=[]
c = 0
d = 0
for line in outputs:
    if(int(line[3]) != 3):
        count_list.append(int(line[3]))

count_list.sort(reverse=True)
#n, bins, patches = plt.hist(count_list)
n, bins, patches = plt.hist(count_list, bins='auto')
plt.xlabel('count')
plt.ylabel('number of reads')
plt.title('distribution')
#plt.show()
plt.savefig(output_file_name)

print(count_list)
print("%d %d"%(c,most_common(count_list)))

