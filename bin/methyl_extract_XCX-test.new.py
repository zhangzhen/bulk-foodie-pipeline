#!/usr/bin/env python

'''
methylation extracting for single molecule footprinting

Created by Jian Fanchong, jfc@pku.edu.cn, May 25 2021
Modified by Shen Ke

tag meaning:
ACA/TGT A/a
ACT/AGT B/b
ACC/GGT C/c
ACG/CGT D/d
TCA/TGA E/e
TCT/AGA F/f
TCC/GGA G/g
TCG/CGA H/h
CCA/TGG I/i
CCT/AGG J/j
CCC/GGG K/k
CCG/CGG L/l
GCA/TGC M/m
GCT/AGC N/n
GCC/GGC O/o
GCG/CGC P/p
.: others

'''

import sys, os
import pysam
import numpy as np
import multiprocessing as mp
import argparse
import time  
import functools  

def timer(func):  
    @functools.wraps(func)  
    def wrapper(*args, **kwargs):  
        start_time = time.time()  
        result = func(*args, **kwargs)  
        end_time = time.time()  
        print(f"{func.__name__} ran in {end_time - start_time:.4f} seconds", file=sys.stderr)  
        return result  
    return wrapper

parser = argparse.ArgumentParser()
parser.add_argument("-fo", "--output", type = str)
parser.add_argument("-a", "--alignment", type = str)
parser.add_argument("-c", "--cpu", type = int, default = 1)
parser.add_argument("-t", "--trimming_Tn5", type = int, default = 15)
parser.add_argument("-q", "--map_quality_filter", type = int, default = 5)
parser.add_argument("--ATAC", type = str, default = 'Y')
parser.add_argument("--remove_chimera", type = str, default = 'Y')
parser.add_argument("--cell_prefix", type=str, required=True)

tag_dict = {
    'ACA':'A', 'ACT':'B', 'ACC':'C', 'ACG':'D', 
    'TCA':'E', 'TCT':'F', 'TCC':'G', 'TCG':'H',
    'CCA':'I', 'CCT':'J', 'CCC':'K', 'CCG':'L',
    'GCA':'M', 'GCT':'N', 'GCC':'O', 'GCG':'P'
}
tags = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPp"
tag2id = {}
baseComplement = {'A':'T','C':'G','T':'A','G':'C'}

def reverseComplement(seq):
    ans = ''
    for i in range(len(seq)-1, -1, -1):
        ans += baseComplement[seq[i]]
    return ans

def get_context(seq, start, end):
    if start >= 0 and end < len(seq):
        con = seq[start:end+1].upper()
        if 'N' in con:
            return 'U'
        if con[1] == 'G':
            con = reverseComplement(con)
        return tag_dict[con]
    else:
        return '.'

def call_methyl(readSeq, refSeq, refPos, refStart, XR, XG):
    is_methylated = ""
    GA_converted_ref_pos,CT_converted_ref_pos = [],[]
    j = refStart
    if XG == "GA":
        for i in range(len(readSeq)):
            if refPos[i] is None:
                if i==len(readSeq)-1:
                    is_methylated += '.'
                continue
            is_methylated += '.'*(refPos[i]-j-1)
            j = refPos[i]
            if refSeq[refPos[i]-refStart] == "g" and readSeq[i] == "A":
                is_methylated += get_context(refSeq, refPos[i]-refStart-1, refPos[i]-refStart+1).lower()
                GA_converted_ref_pos.append(j)
            elif refSeq[refPos[i]-refStart] == "G" and readSeq[i] == "G":
                is_methylated += get_context(refSeq, refPos[i]-refStart-1, refPos[i]-refStart+1).upper()
            elif refSeq[refPos[i]-refStart] == "c" and readSeq[i] == "T":
                is_methylated += "."
                CT_converted_ref_pos.append(j)
            else:
                is_methylated += "."
            
    elif XG == "CT":
        for i in range(len(readSeq)):
            if refPos[i] is None:
                if i==len(readSeq)-1:
                    is_methylated += '.'
                continue
            is_methylated += '.'*(refPos[i]-j-1)
            j = refPos[i]
            if refSeq[refPos[i]-refStart] == "c" and readSeq[i] == "T":
                is_methylated += get_context(refSeq, refPos[i]-refStart-1, refPos[i]-refStart+1).lower()
                CT_converted_ref_pos.append(j)
            elif refSeq[refPos[i]-refStart] == "C" and readSeq[i] == "C":
                is_methylated += get_context(refSeq, refPos[i]-refStart-1, refPos[i]-refStart+1).upper()
            elif refSeq[refPos[i]-refStart] == "g" and readSeq[i] == "A":
                is_methylated += '.'
                GA_converted_ref_pos.append(j)
            else:
                is_methylated += "."
    return is_methylated, GA_converted_ref_pos,CT_converted_ref_pos

def get_counts(methyl):
    ans = np.zeros(len(tags), dtype = int)
    for i in methyl:
        if i in tag2id:
            ans[tag2id[i]] += 1
    return ans

def read2Bed(read1, read2, trim, ATAC,remove_chimera, cell_prefix):
    isRev1 = '+'
    if read1.is_reverse:
        isRev1 = '-'
    # Reference sequence:    A  T  C  G  A  C  T  G  C  A
    # Aligned read:          A  T  -  -  -  -  T  G  C  A
    # [0, 1, None, None, None, None, 6, 7, 8, 9]
    refPos1 = np.array(read1.get_reference_positions(True))
    XR1 = read1.get_tag("XR")
    XG1 = read1.get_tag("XG")
    CB1 = read1.get_tag("CB")
    refStart1 = read1.reference_start
    refEnd1 = read1.reference_end
    refSeq1 = read1.get_reference_sequence()
    readSeq1 = read1.query_sequence
    ch1 = read1.reference_name
    methyl1, GA_converted_ref_pos1,CT_converted_ref_pos1 = call_methyl(readSeq1, refSeq1, refPos1, refStart1, XR1, XG1)
    readName1 = read1.query_name.replace('/', '_')
    assert(len(methyl1) == refEnd1-refStart1)

    isRev2 = '+'
    if read2.is_reverse:
        isRev2 = '-'
    refPos2 = np.array(read2.get_reference_positions(True))
    XR2 = read2.get_tag("XR")
    XG2 = read2.get_tag("XG")
    refStart2 = read2.reference_start
    refEnd2 = read2.reference_end
    refSeq2 = read2.get_reference_sequence()
    readSeq2 = read2.query_sequence
    ch2 = read2.reference_name
    methyl2, GA_converted_ref_pos2,CT_converted_ref_pos2 = call_methyl(readSeq2, refSeq2, refPos2, refStart2, XR2, XG2)
    readName2 = read2.query_name.replace('/', '_')
    assert(len(methyl2) == refEnd2-refStart2)

    assert(ch1 == ch2)
    assert(XG1 == XG2)
    assert(XR1 != XR2)

    l1 = 0
    l2 = 0
    r1 = len(methyl1)
    r2 = len(methyl2)

    if isRev1 == '+':
        refStart1 += trim
        l1 += trim
        refEnd2 -= trim
        r2 -= trim
        if refStart1 > refStart2:
            l2 += (refStart1-refStart2)
            refStart2 = refStart1
        if refEnd1 > refEnd2:
            r1 -= (refEnd1-refEnd2)
            refEnd1 = refEnd2
    else:
        refEnd1 -= trim
        r1 -= trim
        refStart2 += trim
        l2 += trim
        if refEnd2 > refEnd1:
            r2 -= (refEnd2-refEnd1)
            refEnd2 = refEnd1
        if refStart2 > refStart1:
            l1 += (refStart2-refStart1)
            refStart1 = refStart2
    # if reads is too short
    if refStart1>=refEnd1 or refStart2>=refEnd2:
        return None,None,None
    
    # constrain all pos of read1 within refStart1~refEnd1, pos of read2 within refStart2~refEnd2
    GA_converted_ref_pos1,GA_converted_ref_pos2,CT_converted_ref_pos1,CT_converted_ref_pos2 = [np.array(i) for i in [GA_converted_ref_pos1,GA_converted_ref_pos2,CT_converted_ref_pos1,CT_converted_ref_pos2]]
    GA_converted_ref_pos1,GA_converted_ref_pos2,CT_converted_ref_pos1,CT_converted_ref_pos2 = \
    GA_converted_ref_pos1[(GA_converted_ref_pos1>=refStart1)&(GA_converted_ref_pos1<refEnd1)],\
    GA_converted_ref_pos2[(GA_converted_ref_pos2>=refStart2)&(GA_converted_ref_pos2<refEnd2)],\
    CT_converted_ref_pos1[(CT_converted_ref_pos1>=refStart1)&(CT_converted_ref_pos1<refEnd1)],\
    CT_converted_ref_pos2[(CT_converted_ref_pos2>=refStart2)&(CT_converted_ref_pos2<refEnd2)]
    GA_converted_ref_pos,CT_converted_ref_pos = \
        list(GA_converted_ref_pos1)+list(GA_converted_ref_pos2),list(CT_converted_ref_pos1)+list(CT_converted_ref_pos2)
    if len(GA_converted_ref_pos)>0 and len(CT_converted_ref_pos)>0:
        if np.min(GA_converted_ref_pos)>np.max(CT_converted_ref_pos) or np.min(CT_converted_ref_pos)>np.max(GA_converted_ref_pos):
            if remove_chimera == 'Y':
                # remove chimera reads: GA and CT conversion are at 2 sides of read pair
                return None,None,None
    # print('\t',GA_converted_ref_pos,CT_converted_ref_pos)
    if len(GA_converted_ref_pos)+len(CT_converted_ref_pos)==0:
        if ATAC == 'N':
            # remove unconverted reads
            return None,None,None

    methyl1 = methyl1[l1:r1]
    methyl2 = methyl2[l2:r2]
    # print(l1,r1,refStart1,refEnd1,len(methyl1))
    # print('\t',l2,r2,refStart2,refEnd2,refEnd2-refStart2,len(methyl2))
    assert(len(methyl1) == (refEnd1-refStart1))
    assert(len(methyl2) == (refEnd2-refStart2))

    full_cell_barcode = f"{cell_prefix}_{CB1}"

    bedstring1 = '\t'.join(
        [ch1, str(refStart1), str(refEnd1), readName1+'/'+methyl1, '.', isRev1, XG1, full_cell_barcode]
    )
    bedstring2 = '\t'.join(
        [ch2, str(refStart2), str(refEnd2), readName2+'/'+methyl2, '.', isRev2, XG1, full_cell_barcode]
    )
    cnt = get_counts(methyl1) + get_counts(methyl2)
    return cnt, bedstring1, bedstring2

@timer
def main():
    for i in range(len(tags)):
        tag2id[tags[i]] = i
    args = parser.parse_args()
    bamin = args.alignment
    bedout = args.output
    ncpu = args.cpu
    trim = args.trimming_Tn5
    lowq = args.map_quality_filter
    ATAC = args.ATAC
    remove_chimera = args.remove_chimera
    cell_prefix = args.cell_prefix
    # print(ATAC,remove_chimera)
    counts = np.zeros(len(tags), dtype = int)
    refresh_count = 0
    mates = {}
    with open(bedout, "w") as output:
        with pysam.AlignmentFile(bamin, 'rb') as bamfile:
            for read in bamfile: 
                if not read.is_proper_pair:
                    continue
                if read.query_name not in mates:
                    mates[read.query_name] = read
                    continue
                read_ = mates.pop(read.query_name)
                read1, read2 = (read, read_) if read.is_read1 else (read_, read)

                if not (read1.mapping_quality > lowq and read2.mapping_quality > lowq):
                    continue
                # try:
                #     cnt, bedstring1, bedstring2 = read2Bed(read1, read2, trim, ATAC,remove_chimera, cell_prefix)
                # except:
                #     continue
                
                cnt, bedstring1, bedstring2 = read2Bed(read1, read2, trim, ATAC,remove_chimera, cell_prefix)

                if bedstring1==None: # too short reads or chimera reads or unconverted reads
                    continue
                counts += cnt
                output.write(bedstring1+'\n'+bedstring2+'\n')
                
                refresh_count += 1
                if refresh_count % 1000 == 0:
                    output.flush()
                
                # print(refresh_count,end='\r')
                # if refresh_count > 100:
                #     raise ValueError
        
    sys.stdout.write('\t'.join(list(counts.astype(str))))

if __name__ == "__main__":
    main()
