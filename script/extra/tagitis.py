#!/usr/bin/python
"""
Xuenan Pi
27-8-2015
"""

from sys import argv
import copy
import re

def read_file(infile_name):
    """read all results from one software into a list of list, every index
    of sub list indicate chr
    """
    chr_list = [0]*13 
    for i in range(len(chr_list)):
        chr_list[i] = [] 
    infile = open(infile_name)
    for line in infile:
        if line.startswith('SL2.40'):
            chr = int(line.strip().split()[0][-2:])
            start = int(line.strip().split()[1])
            end = int(line.strip().split()[2])
            chr_list[chr].append([start,end])
        else:
            pass
    infile.close()
    return chr_list

def write_file(infile_name,homo_list,filter_list,gap_list,distan):
    """tag the result from itis with insertion
    """
    infile = open(infile_name)
    flag = 0
    for line in infile:
        if line.startswith('SL2.40'):
            line_list = line.strip().split()
            chr = int(line_list[0][-2:])
            start = int(line_list[1])
            end = int(line_list[2])
            homolist = homo_list[chr]
            filterlist = filter_list[chr]
            gaplist = gap_list[chr]
            homo_start_distan = min([abs(start-homo_loci[0]) for homo_loci in homolist])
            homo_end_distan = min([abs(end-homo_loci[1]) for homo_loci in homolist])
            gap_start_distan = min([abs(start-gap_loci[0]) for gap_loci in gaplist])
            gap_end_distan = min([abs(end-gap_loci[1]) for gap_loci in gaplist])
            if homo_start_distan<=distan or homo_end_distan<=distan:
                print "%s\tblast" % line.strip()
            elif gap_start_distan<=distan or gap_end_distan<=distan:
                print "%s\tgap" % line.strip()
            elif [start,end] in filterlist:
                print "%s\tfilter" % line.strip()
            else:
                map_qua = int(re.search(r'MQ=([0-9]*)',line_list[3]).group(1))
                sr = int(re.search(r'SR=([0-9]*)',line_list[3]).group(1))
                if map_qua >= 30 and sr >= 36:
                    print "%s\tmaybe" % line.strip()
                else:
                    print line.strip()

if  __name__=="__main__":
    filterlist = read_file(argv[2])#read the filter file 
    homolist = read_file(argv[3])#read the blast file
    gaplist = read_file(argv[4])#read the gap file
    write_file(argv[1],homolist,filterlist,gaplist,400)
