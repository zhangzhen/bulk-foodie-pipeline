import os, sys
import numpy as np
from scipy.stats import binom_test
import pandas as pd

# Read chromosome sizes once and convert to dictionary for faster lookup
chrom_sizes = pd.read_csv(sys.argv[1], header=None, index_col=0, sep="\t")
chrom_sizes_dict = chrom_sizes[1].to_dict()

# Initialize data structures
regions = {}
N = 0

# Process input line by line
for line in sys.stdin:
    line = line.strip().split('\t')
    if not line:
        continue
        
    ch = line[0]
    st = int(line[1])
    ed = int(line[2])
    
    if ch not in chrom_sizes_dict:
        continue
    
    region = ch
    st_r, ed_r = 0, chrom_sizes_dict[ch]

    # Initialize region data if not exists
    if region not in regions:
        length = ed_r - st_r
        # Use numpy arrays instead of lists
        regions[region] = {
            'name': region,
            'pos': np.zeros(length, dtype=np.int32),
            'neg': np.zeros(length, dtype=np.int32),
            'bases': np.full(length, 'N', dtype='U1')
        }
    r
    # Parse sites and scores
    sites = np.array([int(i) for i in line[4].split(',')], dtype=np.int32)
    scores = np.array([int(i) for i in line[5].split(',')], dtype=np.int32)

    C_or_G = 'C' if line[6] == "CT" else 'G' if line[6] == "GA" else None
    if C_or_G is None:
        raise ValueError(f"Invalid C_or_G value: {line[6]}")

    # Calculate positions and update arrays
    positions = st + sites
    valid_mask = (positions >= st_r) & (positions < ed_r)
    
    if np.any(valid_mask):
        valid_positions = positions[valid_mask]
        valid_scores = scores[valid_mask]
        
        # Update bases
        regions[region]['bases'][valid_positions - st_r] = C_or_Gr
        
        # Update scores using boolean indexing
        pos_mask = valid_scores == 1
        neg_mask = valid_scores == 0
        
        if np.any(pos_mask):
            regions[region]['pos'][valid_positions[pos_mask] - st_r] += 1
        if np.any(neg_mask):
            regions[region]['neg'][valid_positions[neg_mask] - st_r] += 1

    N += 1
    if N % 1000000 == 0:
        sys.stderr.write(f"Processed: {N}\n")

# Output results
for region_data in regions.values():
    region = region_data['name']
    pos = region_data['pos']
    neg = region_data['neg']
    bases = region_data['bases']
    
    # Calculate scores only for positions with data
    data_mask = (pos > 0) | (neg > 0)
    scores = np.full(len(pos), "NA", dtype='U10')
    
    if np.any(data_mask):
        valid_pos = pos[data_mask]
        valid_neg = neg[data_mask]
        scores[data_mask] = (valid_pos / (valid_pos + valid_neg)).astype(str)
    
    # Output results
    for j in range(len(pos)):
        if data_mask[j]:
            sys.stdout.write(f"{region}\t{j}\t{j+1}\t{pos[j]}\t{neg[j]}\t{scores[j]}\t{bases[j]}\n")

