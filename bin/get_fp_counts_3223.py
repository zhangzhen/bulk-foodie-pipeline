#! /usr/bin/env python

import argparse

def output(sites):
    total_depth = 0
    for site in sites:
        total_depth += site[3]
    ave_depth = total_depth / len(sites)

    for site in sites:
        print('%s\t%s\t%s\tid-%d\t%d' % (site[0], site[1], site[2], site[5], round(site[4] * ave_depth)))

def get_fp_counts(handle):
    last_dnase, i = '', 1
    for line in handle:
        chrom, start, end, depth, c_rate, dnase = line.split()
        depth, c_rate = int(depth), float(c_rate)

        if dnase == last_dnase:
            sites.append([chrom, start, end, depth, c_rate, i])
        else:
            if last_dnase != '':
                output(sites)
            last_dnase = dnase
            sites = [[chrom, start, end, depth, c_rate, i]]

        i += 1

    if last_dnase != '':
        output(sites)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required = True, type = str, help = 'input file')
    args = parser.parse_args()
    get_fp_counts(open(args.infile))

