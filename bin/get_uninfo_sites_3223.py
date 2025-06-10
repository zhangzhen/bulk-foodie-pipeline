#! /usr/bin/env python

import argparse

def output(dnase, sites):
    chrom, start, end = dnase
    while start < end:
        if start not in sites:
            print('%s\t%d\t%d\t.\t100' % (chrom, start, start + 1))
        start += 1

def get_uninfo_sites(handle, dnase_handle):
    dnases = {}
    for line in dnase_handle:
        chrom, start, end, dnase_id = line.split()
        start, end = int(start), int(end)
        dnases[dnase_id] = [chrom, start, end]

    last_dnase_id = ''
    for line in handle:
        chrom, start, end, depth, c_rate, dnase_id = line.split()
        start = int(start)

        if dnase_id == last_dnase_id:
            sites.append(start)
        else:
            if last_dnase_id != '':
                output(dnases[last_dnase_id], sites)
            last_dnase_id = dnase_id
            sites = [start]

    if last_dnase_id != '':
        output(dnases[last_dnase_id], sites)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required = True, type = str, help = 'input file')
    parser.add_argument('-d', '--dnase', required = True, type = str, help = 'DNase file')
    args = parser.parse_args()
    get_uninfo_sites(open(args.infile), open(args.dnase))

