#! /usr/bin/env python
import os,sys

sub = sys.argv[2]
bkgd = sys.argv[3]

def getsum(string, sub, bkgd):
    pos = []
    score = []
    sites_num = 0
    converted_sites_num = 0

    for i in range(len(string)):
        if string[i] in sub:
            score.append('1')
            pos.append(str(i))
            sites_num += 1
            converted_sites_num += 1
        elif string[i] in bkgd:
            score.append('0')
            pos.append(str(i))
            sites_num += 1
    # # delete unconverted reads
    # return (','.join(pos), ','.join(score), sites_num) if converted_sites_num>0 else ('.','.',0)
    return (','.join(pos), ','.join(score), sites_num) # methyl_extract_XCX-test.py has already control whether to delete unconverted reads or not, in --ATAC==Y/N parameter

with open(sys.argv[1], 'r') as bedin:
    line1 = bedin.readline()
    line2 = bedin.readline()
    while line2:
        r1 = line1.strip().split('\t')
        r2 = line2.strip().split('\t')
        
        ch1,st1,ed1 = r1[:3]
        ch2,st2,ed2 = r2[:3]
        
        names1,meth1 = r1[3].split('/')
        names2,meth2 = r2[3].split('/')
        
        XG1, XG2 = r1[6], r2[6]

        cell_barcode = r1[7]
        
        st1,ed1,st2,ed2 = list(map(int,[st1,ed1,st2,ed2]))
        st,ed = min(st1,st2),max(ed1,ed2)

        assert(names1==names2 and ch1==ch2 and XG1==XG2)
        
        meth = ['.']*(ed-st)
        
        # methylation of read1 is preferred
        meth[st2-st:ed2-st] = meth2
        meth[st1-st:ed1-st] = meth1
        
        pos, score, sites_num = getsum(meth, sub, bkgd)
        
        st1_ed1_st2_ed2 = (f'{st1},{ed1},{st2},{ed2}' if st1<=st2 else f'{st2},{ed2},{st1},{ed1}')
    
        if sites_num > 0:
            sys.stdout.write("\t".join([ch1, str(st), str(ed), names1, pos, score, XG1, str(sites_num), st1_ed1_st2_ed2, cell_barcode]) + "\n")

        line1 = bedin.readline()
        line2 = bedin.readline()
