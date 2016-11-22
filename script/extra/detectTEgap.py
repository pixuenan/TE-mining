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

def write_file(infile_name,homo_list,gap_list,distan):
    """tag the result from itis with insertion
    """
    infile = open(infile_name)
    flag = 0
    gap_group = [0]*13 
    for i in range(len(gap_group)):
        gap_group[i] = [] 
    for line in infile:
        if line.startswith('SL2.40'):
            line_list = line.strip().split()
            chr = int(line_list[0][-2:])
            start = int(line_list[1])
            end = int(line_list[2])
            homolist = homo_list[chr]
            gaplist = gap_list[chr]
            homo_start_distan = min([abs(start-homo_loci[0]) for homo_loci in homolist])
            homo_end_distan = min([abs(end-homo_loci[1]) for homo_loci in homolist])
            gap_start_distan = min([abs(start-gap_loci[0]) for gap_loci in gaplist])
            gap_end_distan = min([abs(end-gap_loci[1]) for gap_loci in gaplist])
            if homo_start_distan<=distan or homo_end_distan<=distan:
                pass
            elif gap_start_distan<=distan or gap_end_distan<=distan:
                gap_group[chr] += [[start,end]]
            else:
                pass
    print '###############'
    for i in range(13):
        result = group_gap(gap_group[i],6000)
        for j in result:
            print 'SL2.40ch%s\t%d\t%d\t%d' % (i,j[0],j[1],j[1]-j[0])

def group_gap(chr_list,distan):
    """put all putative in a group"""
    result = []
    for ind,loci in enumerate(chr_list[:-1]):
        if chr_list[ind+1][1]-loci[0]<= distan:
            result += [[loci[0],chr_list[ind+1][1]]]
    return result

if  __name__=="__main__":
    homolist = read_file(argv[2])#read the blast file
    gaplist = read_file(argv[3])#read the gap file
    write_file(argv[1],homolist,gaplist,400)
