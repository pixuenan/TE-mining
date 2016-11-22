#!/usr/bin/python
"""
Xuenan Pi
12-8-2015
Script to print the common
"""

from sys import argv
import copy

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
            loci = int(line.strip().split()[1])
            chr_list[chr] += [loci]
        else:
            pass
    infile.close()
    return chr_list

def print_file(chr_list,filename):
    """read all results from one software into a list of list, every index
    of sub list indicate chr
    """
    infile = open(filename)
    for line in infile:
        if line.startswith('SL2.40'):
            chr = int(line.strip().split()[0][-2:])
            loci = int(line.strip().split()[1])
            for chr_i,chr_l in enumerate(chr_list):
                for loc in chr_l:
                    if chr==chr_i and loci==loc:
                        print line
    return 

if  __name__=="__main__":
    chr_list = read_file(argv[1])
    print_file(chr_list,argv[2])

