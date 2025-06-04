#! /usr/bin/env python

import os, sys
from collections import OrderedDict

current_chrom = None
current_data = None
N = 0

def process_chrom(chrom, data):
    for pos, (pos_count, neg_count, base) in data.items():
        if pos_count > 0 or neg_count > 0:
            score = pos_count/(pos_count + neg_count)
        else:
            score = "NA"
        sys.stdout.write("\t".join([chrom, str(pos), str(pos+1), str(pos_count), str(neg_count), str(score), base])+"\n")

line = sys.stdin.readline().strip()
while line:
    line = line.split('\t')
    ch = line[0]
    st = int(line[1])
    ed = int(line[2])
    
    # If we're moving to a new chromosome, process and output the previous one
    if ch != current_chrom:
        if current_chrom is not None:
            process_chrom(current_chrom, current_data)
        current_chrom = ch
        current_data = OrderedDict()
    
    sites = [int(i) for i in line[4].split(',')]
    scores = [int(i) for i in line[5].split(',')]

    C_or_G = line[6]
    if C_or_G == "CT":
        C_or_G = "C"
    elif C_or_G == "GA":
        C_or_G = "G"
    else:
        raise ValueError

    for i in range(len(sites)):
        site = sites[i]
        score = scores[i]
        pos = st + site
        
        if pos not in current_data:
            current_data[pos] = [0, 0, C_or_G]
        
        if score == 1:
            current_data[pos][0] += 1
        elif score == 0:
            current_data[pos][1] += 1
        else:
            raise ValueError

    N += 1
    if (N % 1000000 == 0):
        sys.stderr.write("Processed: "+str(N)+"\n")
    line = sys.stdin.readline().strip()

# Process the last chromosome
if current_chrom is not None:
    process_chrom(current_chrom, current_data)