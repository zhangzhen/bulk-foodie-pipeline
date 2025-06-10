#! /usr/bin/env python

import argparse
import tabix
import pandas as pd

def is_open_region(site_handle, win_size, CG_num, ratio, high_num, chrom, scan_start):
    try:
        sites = site_handle.query(chrom, scan_start, scan_start+win_size)
    except:
        return False
    
    C_num, G_num, high_ratio_num = 0, 0, 0
    for site in sites:
        if site[3] == 'C':
            C_num += 1
        else:
            G_num += 1
        
        if float(site[5]) >= ratio:
            high_ratio_num += 1
            
    if (C_num >= CG_num) and (G_num >= CG_num) and (high_ratio_num >= high_num):
        return [scan_start, scan_start+win_size]
    else:
        return False

def trim_open_region(site_handle, win_size, ratio, border_num, chrom, open_region):
    sites = site_handle.query(chrom, open_region[0], open_region[1])
    sites = pd.DataFrame(sites, columns = ['chrom', 'start', 'end', 'site_type', 'depth', 'ratio'])
    sites = sites[['end', 'ratio']]
    sites = sites.astype({'end': int, 'ratio': float})
    
    sites_ris, border_start_num, border_end_num = list(range(sites.shape[0])), 0, 0
    for sites_ri in sites_ris:
        if sites.iloc[sites_ri, 1] >= ratio:
            border_start_num += 1
            if border_start_num > border_num:
                trimmed_start = sites.iloc[sites_ri, 0] - 1
                break
                
    for sites_ri in sites_ris[::-1]:
        if sites.iloc[sites_ri, 1] >= ratio:
            border_end_num += 1
            if border_end_num > border_num:
                trimmed_end = sites.iloc[sites_ri, 0]
                break

    if trimmed_end - trimmed_start >= win_size:
        return [trimmed_start, trimmed_end]
    else:
        return False
    
def output_open_regions(chrom, open_regions):
    last_open_region = None
    for open_region in open_regions:
        if last_open_region:
            if open_region[0] <= last_open_region[1]:
                last_open_region[1] = open_region[1]
            else:
                print('%s\t%d\t%d' % (chrom, last_open_region[0], last_open_region[1]))
                last_open_region = open_region
        else:
            last_open_region = open_region
    if last_open_region:
        print('%s\t%d\t%d' % (chrom, last_open_region[0], last_open_region[1]))

def get_open_region(site_handle, peak_handle, win_size, CG_num, ratio, high_num, border_num):
    for line in peak_handle:
        chrom, start, end, name = line.split()
        start, end, last_open_region, trimmed_open_regions = int(start), int(end), None, []
        for scan_start in range(start, end-win_size):
            open_region = is_open_region(site_handle, win_size, CG_num, ratio, high_num, chrom, scan_start)
            if open_region:
                if last_open_region:
                    if last_open_region[1] == scan_start+win_size-1:
                        last_open_region[1] += 1
                    else:
                        trimmed_open_region = trim_open_region(site_handle, win_size, ratio, border_num, chrom, last_open_region)
                        if trimmed_open_region:
                            trimmed_open_regions.append(trimmed_open_region)
                        last_open_region = open_region
                else:
                    last_open_region = open_region
        if last_open_region:
            trimmed_open_region = trim_open_region(site_handle, win_size, ratio, border_num, chrom, last_open_region)
            if trimmed_open_region:
                trimmed_open_regions.append(trimmed_open_region)
        output_open_regions(chrom, trimmed_open_regions)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sitefile', required = True, type = str, help = 'site file')
    parser.add_argument('-p', '--peakfile', required = True, type = str, help = 'peak file')
    parser.add_argument('-w', '--winsize', type = int, default = 100, help = 'window size, default 100 bp')
    parser.add_argument('-c', '--CGnum', type = int, default = 10, help = 'min. C/G site number, default 10')
    parser.add_argument('-r', '--ratio', type = float, default = 0.8, help = 'min. ratio of conv. rate to pos. control, default 0.8')
    parser.add_argument('-n', '--highnum', type = int, default = 10, help = 'min. high-ratio site number, default 10')
    parser.add_argument('-b', '--bordernum', type = int, default = 2, help = 'min. number of high-ratio sites outside border, default 2')
    args = parser.parse_args()

    site_handle = tabix.open(args.sitefile)
    peak_handle = open(args.peakfile)
    get_open_region(site_handle, peak_handle, args.winsize, args.CGnum, args.ratio, args.highnum, args.bordernum)

