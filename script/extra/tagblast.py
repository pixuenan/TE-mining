#!/usr/bin/python
"""
Xuenan Pi
3-9-2015
Script to collapse the fragment result from blast to a single copy of TE.
Another file come from run ITIS on simulated reference data
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

def group_frag(chr_list,distan):
    """put all putative in a group"""
    result_list = []
    p_list = [chr_list[0]] 
    i = 1
    while i < len(chr_list):
        loci = chr_list[i]
        result = check_line(p_list,loci,distan)
        if result[0]:
            p_list += [loci]
            i += 1
        else:
            result_list += [p_list]
            p_list = [loci]
            i += 1
    return result_list

def check_line(p_list,loci,distan):
    """check if the fragments may belong one TE"""
    if loci[0]-p_list[-1][1] <= distan:
        return True,loci
    else:
        return False,loci

def uniq(lst):
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item

def write_file(i,itis_list,group_list,distan):
    """tag the result from blast with copy
    """
    twoside_list = []
    oneside_list = []
    no_list = []
    x = 0
    while x < len(group_list):
        loci_list = group_list[x]
        min_start_distan = min([abs(loci_list[0][0]-itis_loci[0]) for itis_loci in itis_list])
        min_end_distan = min([abs(loci_list[-1][-1]-itis_loci[1]) for itis_loci in itis_list])
        if min_start_distan<=distan and min_end_distan<=distan:
            twoside_list.append([loci_list[0][0],loci_list[-1][-1],"COPY"])
            group_list.remove(loci_list)
        elif min_start_distan<=distan or min_end_distan<=distan:
            oneside_list.append([loci_list[0][0],loci_list[-1][-1],"COPY 1"])
            group_list.remove(loci_list)
        else:
            no_list.append([loci_list[0][0],loci_list[-1][-1],""])
            x += 1
    outlist = list(uniq(sorted(no_list+twoside_list+oneside_list)))
    for j in outlist:
        if i < 10:
            print "SL2.40ch0%d\t%d\t%d\t%d\t%s" % (i,j[0],j[1],j[1]-j[0],j[2])
        else:
            print "SL2.40ch%d\t%d\t%d\t%d\t%s" % (i,j[0],j[1],j[1]-j[0],j[2])

if  __name__=="__main__":
    itislist = read_file(argv[1])
    homolist = read_file(argv[2])
    rawlist = []
    homogroup = []
    copylist = []
    num = 0
    for i in range(13):
        homogroup.append(group_frag(homolist[i],5000))
    for i in range(13):
        write_file(i,itislist[i],homogroup[i],400)
