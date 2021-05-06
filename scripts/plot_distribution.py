import matplotlib.pyplot as plt
import pandas as pd
import argparse
import seaborn as sns

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
    outputs.append(line.rstrip().split()[0:7])

count_list=[]
c = 0
d = 0
for line in outputs:
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
print("%d %d"%(c,d))