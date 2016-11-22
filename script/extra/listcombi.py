#!/usr/bin/python
"""
Xuenan Pi
28-8-2015
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

def comp_list(list1,num1=None,list2=None,num2=None):
    counter = 0
    if not list2:
        for chr in list1:
            for loci in chr:
                loci.append([num1])
    else:
        for inx,chr in enumerate(list1):
            for loci in chr:
                if loci[0] in [list2loci[0] for list2loci in list2[inx]]\
                        and loci[1] in [list2loci[1] for list2loci in list2[inx]]:
                    loci[-1].append(num2)
                    list2[inx].remove([loci[0],loci[1]])
                else:
                    counter += 1
        if counter:
            for inx,chr in enumerate(list2):
                for loci in chr:
                    list1[inx].append([loci[0],loci[1],[num2]])
    return list1

def print_file(list):
    """Write the insertions to a file.
    """
    chr_name_list = ['SL2.40ch00','SL2.40ch01','SL2.40ch02','SL2.40ch03','SL2.40ch04','SL2.40ch05','SL2.40ch06','SL2.40ch07','SL2.40ch08','SL2.40ch09','SL2.40ch10','SL2.40ch11','SL2.40ch12']
    for index,chr_list in enumerate(list):
        if chr_list:
            chr = chr_name_list[index]
            for loci in chr_list:
                print "%s\t%d\t%d\t%s\t%d" % (chr,loci[0],loci[1],loci[2],len(loci[2]))



if  __name__=="__main__":
    chrlist1 = read_file(argv[1])
    chrlist2 = read_file(argv[2])
    chrlist3 = read_file(argv[3])
    list1 = comp_list(chrlist1,"01")
    list2 = comp_list(list1,None,chrlist2,"11")
    list3 = comp_list(list2,None,chrlist3,"19")
    print_file(list3)
