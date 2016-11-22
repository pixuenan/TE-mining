#!/usr/bin/python
"""
Xuenan Pi
28-8-2015
Script to detect non-reference TIPs
"""

from sys import argv
from functions import read_file_2,initial_list
import copy
import re

def combi_list(dict1,dict2=None,i=0,distan=0):
    """combine results from different accessions to detect TIPs"""
    counter = 0
    if not dict2:#initial the first list
        for chr_name,chr_list in dict1.items():
            for loci in chr_list:
                if loci[2]=='introgression':
                    loci.append(['0','i'])
                else:
                    loci.append(['0','1'])
    else:
        for chr_name,chr_list in dict1.items():
            if (not dict1[chr_name] and dict2[chr_name]):
                counter += 1
            for loci in chr_list:
                if dict2[chr_name]:
                    min_distan = min([abs(loci[1]-list2loci[0])\
                                      for list2loci in dict2[chr_name]]) 
                    if min_distan < distan:
                        for loci2 in dict2[chr_name]:
                            if abs(loci[1]-loci2[0])==min_distan:
                                break
                        if loci2[2]=='introgression':
                            loci[-1].append('i')
                        else:
                            loci[-1].append('1')
                        dict2[chr_name].remove(loci2)
                    else:
                        loci[-1].append('0')
                        counter += 1
                else:
                    loci[-1].append('0')
        if counter:#the location that is only present in the second list
            for chr_name,chr_list in dict2.items():
                for loci in chr_list:
                    if loci[2]=='introgression':
                        dict1[chr_name].append(loci[:3]+[list((i+1)*'0')+\
                                               ['i']])
                    else:
                        dict1[chr_name].append(loci[:3]+[list((i+1)*'0')+\
                                               ['1']])
    for chr_list in dict1.values():
        chr_list.sort()
    return dict1

def combi_insertion_list(insertion_list,distan):
    """detect TIPs by accession and accession"""
    result_list = initial_list(len(insertion_list))
    result_list[0] = combi_list(insertion_list[0])
    for inx in range(1,len(result_list)):
        result_list[inx] = combi_list(result_list[inx-1],insertion_list[inx],\
                                      inx,distan)
    return result_list[-1]

def print_file(chr_dict):
    """Write the output.
    """
    chr_name_list = sorted(chr_dict.keys())
    for chr_name in chr_name_list:
        chr_list = chr_dict[chr_name]
        if chr_list:
            for loci in chr_list:
                print "%s\t%d\t%d\tinsertion\t%d\t%s\t%d" % (chr_name,loci[0],\
                      loci[1],loci[1]-loci[0],'\t'.join(loci[3]),\
                      len(loci[3])-loci[3].count('0'))

if  __name__=="__main__":
    te_name = argv[1]
    access_num_ls = argv[2].split(',')#001,011,019
    important_data_path = argv[3]
    insert_size = int(argv[4])
    insertion_list = initial_list(len(access_num_ls))
    for i,access_num in enumerate(access_num_ls):
        insertion_file = '%s/result/%s/rf%s/itisrf%s.insertion.anno.bed'\
        % (important_data_path,te_name,access_num,access_num)
        insertion_list[i] = read_file_2(insertion_file)
    max_distan = insert_size
    print_file(combi_insertion_list(insertion_list,max_distan))
