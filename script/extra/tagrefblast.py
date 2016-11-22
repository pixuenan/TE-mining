#!/usr/bin/python
"""
Xuenan Pi
15-9-2015
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
            strand = line.strip().split()[-1]
            chr_list[chr].append([start,end,strand])
        else:
            pass
    infile.close()
    return chr_list

def group_frag(given_list,distan):
    """put all putative in a group"""
    p_list = [given_list[0]]
    if len(given_list)==1:
        return given_list
    else:
        i = 1
        while i < len(given_list):
            loci = given_list[i]
            result = check_line(p_list,loci,distan)
            if result[0]:
                p_list += [loci]
                if i+1==len(given_list):
                    return p_list
                else:
                    i += 1
            else:
                return p_list

def check_line(p_list,loci,distan):
    """check if the fragments may belong one TE"""
    if loci[1]-p_list[0][0] <= distan:
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

def itis_loci_finder(itis_list,loci_num,loci,min_distan):
    for ind,itis_loci in enumerate(itis_list):
        if abs(itis_loci[loci_num]-loci)==min_distan:
            loci = itis_loci
            strand = itis_loci[-1]
            break
    return ind,loci,strand

def find_copy(i,itis_list,chr_list,distan,max_len):
    """tag the result from blast with copy
    """
    twoside_list = []
    oneside_list = []
    no_list = []
    x = 0
    while x < len(chr_list):
        loci_list = group_frag(chr_list[x:],max_len)
        min_start_distan = min([abs(loci_list[0][0]-itis_loci[0]) for itis_loci in itis_list])
        start_ind,start_loci,start_strand = itis_loci_finder(itis_list,0,loci_list[0][0],min_start_distan)
        flag = False
        if min_start_distan<=distan:
            if start_ind == len(itis_list)-1:
                start_ind = start_ind-1
            for i,loci in enumerate(loci_list):
                end_distan = abs(loci[1]-itis_list[start_ind+1][1])
                if end_distan <=distan and start_strand==itis_list[start_ind+1][2]:
                    flag = True
                    break
            if flag:#one copy has been confirmed from both end
                twoside_list = process_list(loci_list[:i+1],'COPY',start_strand,twoside_list)
                x += i+1
            elif itis_list[start_ind+1][0] > loci_list[-1][1]+distan:#one copy is only confirmed at end side, there is no sr disrupt the group
                oneside_list = process_list(loci_list,'COPY 1',start_strand,oneside_list)
                x += len(loci_list)
            else:#the group cannot be confirmed at either side
                no_list.append([loci_list[0][0],loci_list[-1][1],"",""])
                x += 1
        else:
            min_end_distan = min([abs(loci_list[-1][1]-itis_loci[1]) for itis_loci in itis_list])
            end_ind,end_loci,end_strand = itis_loci_finder(itis_list,1,loci_list[-1][1],min_end_distan)
            if min_end_distan<=distan:
                if itis_list[end_ind-1][1] < loci_list[0][0]-distan:#one copy is only confirmed at end side
                    oneside_list = process_list(loci_list,'COPY 1',end_strand,oneside_list)
                    x += len(loci_list)
                else:#the group cannot be confirmed at either side
                    no_list.append([loci_list[0][0],loci_list[-1][1],"",""])
                    x += 1
            else:#the group cannot be confirmed at either side
                no_list.append([loci_list[0][0],loci_list[-1][1],"",""])
                x += 1
    outlist = list(uniq(sorted(twoside_list+oneside_list)))
    return outlist

def process_list(loci_list,tag,strand,copy_list):
    copy_list.append([loci_list[0][0],loci_list[-1][1],loci_list[-1][1]-loci_list[0][0],tag,strand])
    return copy_list

def write_file(i,outlist,min_len):
    chr_name_list = ['SL2.40ch00','SL2.40ch01','SL2.40ch02','SL2.40ch03','SL2.40ch04','SL2.40ch05','SL2.40ch06','SL2.40ch07','SL2.40ch08','SL2.40ch09','SL2.40ch10','SL2.40ch11','SL2.40ch12']
    for loci in outlist:
        chr_name = chr_name_list[i]
        if loci[2]>=min_len:
            print '%s\t%s' % (chr_name,'\t'.join(map(str,loci)))

if  __name__=="__main__":
    itislist = read_file(argv[1])
    homolist = read_file(argv[2])
    max_len = int(argv[3])
    min_len = int(argv[4])
    homo_copy_list = []
    num = 0
    for i in range(13):
        if itislist[i]:
            homo_copy_list.append(find_copy(i,itislist[i],homolist[i],400,max_len))
        else:
            homo_copy_list.append([])
    for i in range(13):
        write_file(i,homo_copy_list[i],min_len)
