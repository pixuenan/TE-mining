#!/usr/bin/python
"""
Xuenan Pi
14-9-2015
Script to detect reference TIPs
"""

from sys import argv
from functions import read_file_2,initial_list
import copy
import re
import subprocess

def combi_list(dict1,distan,del_dict=None,dict2=None):
    """combine results from different accessions to detect TIPs"""
    if not dict2:#load the copy loci from ref
        for chr_name,chr_list in dict1.items():
            for loci in chr_list:
                if loci[2]=='introgression':
                    loci.append(['i'])
                else:
                    loci.append(['1'])
    else:
        for chr_name,chr_list in dict1.items():
            for loci in chr_list:
                flag = True
                if dict2[chr_name]:
                    for loci2 in dict2[chr_name]:
                        if loci[:2] == loci2[:2]:
                            flag = False
                            if loci2[2]=='introgression':
                                loci[-1].append('i')
                            else:
                                loci[-1].append('1')
                        dict2[chr_name].remove(loci2)
                if flag:#deletion
                    dellist = del_dict[chr_name]
                    #check if the deletion can be confirmed
                    del_start_distan = min([abs(loci[0]-del_loci[0]) \
                                       for del_loci in dellist])
                    del_end_distan = min([abs(loci[1]-del_loci[1]) \
                                     for del_loci in dellist])
                    if del_start_distan<=distan and del_end_distan<=distan:
                        loci[-1].append('confirmed')
                    else:
                        loci[-1].append('unconfirmed')
    for chr_list in dict1.values():
        chr_list.sort()
    return dict1

def combi_copy_list(copy_list,del_list,distan):
    """detect TIPs by accession and accession"""
    result_list = initial_list(len(copy_list))
    result_list[0] = combi_list(copy_list[0],distan)
    for inx in range(1,len(result_list)):
        result_list[inx] = combi_list(result_list[inx-1],distan,del_list[inx],\
                                      copy_list[inx])
    return result_list[-1]

def print_file(chr_dict,flag):
    """Write the output.
    """
    chr_name_list = sorted(chr_dict.keys())
    for chr_name in chr_name_list:
        chr_list = chr_dict[chr_name]
        if chr_list:
            for loci in chr_list:
                if flag=='gap':
                    print "%s\t%d\t%d\tNgap\t%d\t%s\t%d" % (chr_name,loci[0],\
                            loci[1],loci[1]-loci[0],'\t'.join(loci[3]),\
                            loci[3].count('1')+loci[3].count('unconfirmed'))
                else:
                    print "%s\t%d\t%d\tcopy\t%d\t%s\t%d" % (chr_name,loci[0],\
                            loci[1],loci[1]-loci[0],'\t'.join(loci[3]),\
                            loci[3].count('1')+loci[3].count('unconfirmed'))

if  __name__=="__main__":
    te_name = argv[1]
    access_num_ls = argv[2].split(',')#000,001,011,019
    flag = argv[3]
    important_data_path = argv[4]
    insert_size = int(argv[5])
    copy_list = initial_list(len(access_num_ls))
    del_list = initial_list(len(access_num_ls))
    max_distan = insert_size
    for i,access_num in enumerate(access_num_ls):
        deletion_file = '%s/store/deletion/rf%s_deletion.info' % \
                (important_data_path,access_num)
        if flag=='gap':
            copy_file = '%s/result/%s/rf%s/itisrf%s.gap.copy.corr.bed' % \
                    (important_data_path,te_name,access_num,access_num)
        else:
            copy_file ='%s/result/%s/rf%s/itisrf%s.copy.corr.bed' % \
                    (important_data_path,te_name,access_num,access_num)
        copy_list[i] = read_file_2(copy_file)
        if access_num!='000':
            del_list[i] = read_file_2(deletion_file)
    print_file(combi_copy_list(copy_list,del_list,max_distan),flag)
