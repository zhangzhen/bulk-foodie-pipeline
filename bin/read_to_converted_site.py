import os, sys
import numpy as np
from scipy.stats import binom_test
import pandas as pd

line = sys.stdin.readline().strip()
chrom_sizes = pd.read_csv(sys.argv[1], header=None, index_col=0, sep="\t")

regions = {}
N = 0

while line:
    line = line.split('\t')
    ch = line[0]
    st = int(line[1])
    ed = int(line[2])
    
    if ch not in list(chrom_sizes.index):
        line = sys.stdin.readline().strip()
        continue    
    
    region = ch
    st_r, ed_r = 0, chrom_sizes.loc[ch,1]

    if region not in regions:
        length = ed_r - st_r
        regions[region] = [region, [0]*length, [0]*length, ["N"]*length]
    
    # sites = [int(i) for i in line[4].split(',')[1:-1]]
    # scores = [int(i) for i in line[5].split(',')[1:-1]]
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
        pos = st+site
        if pos >= st_r and pos < ed_r:
            regions[region][3][pos-st_r] = C_or_G
            if score == 1:
                regions[region][1][pos-st_r] += 1
            elif score == 0:
                regions[region][2][pos-st_r] += 1
            else:
                raise ValueError
    N += 1
    if (N%1000000 == 0):
        sys.stderr.write("Processed: "+str(N)+"\n")
    line = sys.stdin.readline().strip()

for i in regions:
    region = regions[i]
    length = chrom_sizes.loc[region[0],1]
    pos = region[1]
    neg = region[2]
    for j in range(length):
        if (pos[j] > 0 or neg[j] > 0):
            score = pos[j]/(pos[j]+neg[j])
        else:
            score = "NA"
#            score = binom_test(pos[j], n=pos[j]+neg[j], p=0.085,alternative="greater")
#            if score < 1e-20:
#                score = 20
#            else:
#                score = -np.log10(score)
        sys.stdout.write("\t".join([region[0],str(j),str(j+1),str(pos[j]),str(neg[j]),str(score),region[3][j]])+"\n")

