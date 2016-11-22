#!/usr/bin/python
"""
Xuenan Pi
7-8-2015
Script to correct the homologs insert locations from blast by ITIS.
"""

from sys import argv
import copy

def read_file(infile_name):
    """read all results from list"""
    chr_list = [0]*13
    for i in range(len(chr_list)):
        chr_list[i] = [] 
    infile = open(infile_name)
    for line in infile:
        if line.startswith('SL2.40'):
            chr = int(line.strip().split()[0][-2:])
            loci = map(int,line.strip().split('\t')[1:-1])+[line.strip().split('\t')[-1]]
            chr_list[chr] += [loci]
        else:
            pass
    infile.close()
    return chr_list

def merge(i,homo_list,distan):
    """merge the blast result should belong to one result"""
    result = []
    group = [homo_list[0]]
    for inx,loci in enumerate(homo_list[1:]):
        if loci[3]=='COPY 1'and loci[2]<=4000:#considered
            if group[0][2]>=4000:#remove the first elem
                result += group
                group = [loci]
            elif loci[1]-group[0][0]<=distan:
                group += [loci]
            else:
                if len(group)>1:
                    result += [[group[0][0],group[-1][1],group[-1][1]-group[0][0],'COPY 2']]
                else:
                    result += group
                group = [loci]
        else:
            if len(group)>1:
                result += [[group[0][0],group[-1][1],group[-1][1]-group[0][0],'COPY 2']]
            else:
                result += group
            group = [loci]
    if len(group)>1:
        result += [[group[0][0],group[-1][1],group[-1][1]-group[0][0],'COPY 2']]
    else:
        result += group
    chr_name_list = ['SL2.40ch00','SL2.40ch01','SL2.40ch02','SL2.40ch03','SL2.40ch04','SL2.40ch05','SL2.40ch06','SL2.40ch07','SL2.40ch08','SL2.40ch09','SL2.40ch10','SL2.40ch11','SL2.40ch12']
    for loci in result:
        chr_name = chr_name_list[i]
        print '%s\t%s' % (chr_name,'\t'.join(map(str,loci)))

if  __name__=="__main__":
    file1 = argv[1]
    homo_list = read_file(file1)
    for i in range(13):
        merge(i,homo_list[i],30000)
    #mer_list = merge(list,30000)
 #   cut_list = cut(mer_list,800)
    #write(mer_list,outfilename)
