#! /usr/bin/env python
import os,sys

sub = sys.argv[2]

def getsum(string, sub):
    ans = 2
    tot = 0
    pos = ['0']
    full = sub.lower() + sub.upper()
    for i in range(len(string)):
        if string[i] in sub:
            ans += 1
            tot += 1
            pos.append(str(i+1))
        elif string[i] in full:
            tot += 1
    pos.append(str(len(string)+1))
    return ans, tot, ','.join(pos)

with open(sys.argv[1], 'r') as bedin:
    line = bedin.readline()
    while line:
        tmp = line.strip().split('\t')
        try:
            tmp[4] = int(tmp[4])
        except:
            tmp[4] = 0
        tmp[4] = str(tmp[4])
        names = tmp[3].split('/')
        meth = names[-1]
        # print(len(meth), int(tmp[2])-int(tmp[1]))
        assert(len(meth) == int(tmp[2])-int(tmp[1]))
        meth = meth[1:-1]
        tmp.append(str(int(tmp[1])+1))
        tmp.append(str(int(tmp[2])-1))
        tmp.append('0')
        s,tot,pos = getsum(meth, sub)
        names[-1] = str(s-2)+','+str(tot)
        tmp[3] = '/'.join(names)
        tmp.append(str(s))
        tmp.append(','.join(['1']*s))
        tmp.append(pos)
        sys.stdout.write("\t".join(tmp))
        line = bedin.readline()
        if line:
            sys.stdout.write("\n")
