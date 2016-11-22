#!/usr/bin/python
"""
Xuenan Pi
13-8-2015
Script to extract the paired reads around an insertion
"""

from sys import argv
import copy
import re

def read_file(infile_name,dis):
    left = []
    right = []
    infile = open(infile_name)
    for line in infile:
        if line.startswith('@'):
            pass
        else:
            info_list = line.strip().split()
            chr = info_list[2]
            loci = int(info_list[3])
            if chr=='SL2.40ch01' and info_list[6]=='Rider':
                if loci+5-dis < 500 and loci+5-dis > 0:
                    right += [info_list[:2]]
                elif dis-5-loci < 500 and dis-5-loci > 0:
                    left += [info_list[:2]]
    infile.close()
    return right,left

def extract_fastq(infile_name,list):
    infile = open(infile_name)
    for line in infile:
        if line.startswith('@'):
            pass
        else:
            info_list = line.strip().split()
            pair = info_list[1]
            read_n = info_list[0]
            for loci in list:
                if loci[0]==read_n:
                    if int(loci[1][-1])+int(pair[-1])==3:
                        print '>'+read_n
                        print info_list[9]
    infile.close()
    return


if  __name__=="__main__":
    right,left = read_file(argv[1],15882500)
    extract_fastq(argv[1],left)
