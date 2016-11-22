#!/usr/bin/python
"""
Xuenan Pi
7-8-2015
Script to correct the homologs insert locations from blast by ITIS.
"""

from sys import argv
from functions import read_file_1
import copy

def merge(chr_dict):
    """merge the blast results that overlap with each other"""
    for chr_name,chr_list in chr_dict.items():
        inx = 1
        while inx < len(chr_list):
            if chr_list[inx][0] - chr_list[inx-1][1] <= 0:
                loci_end = chr_list[inx][1]
                chr_dict[chr_name][inx-1][1]= loci_end
                chr_dict[chr_name].remove(chr_list[inx])
            else:
                inx += 1
    return chr_dict

def write(chr_dict,filename):
    """write the result to a new file"""
    chr_name_list = sorted(chr_dict.keys())
    file = open(filename,'w+')
    for chr_name in chr_name_list:
        chr_list = chr_dict[chr_name]
        for i in range(len(chr_list)):
            file.write('%s\t%s\t%d\n' % (chr_name,'\t'.join(map(str,\
                        chr_list[i])),chr_list[i][1]-chr_list[i][0]))
    file.close()

if  __name__=="__main__":
    infile = argv[1]
    outfilename = argv[2]
    chr_dict = read_file_1(infile)
    merge_dict = merge(chr_dict)
    write(merge_dict,outfilename)
