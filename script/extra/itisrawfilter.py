#!/usr/bin/python
"""
Xuenan Pi
13-8-2015
Script to filter the raw result from ITIS 
"""

from sys import argv
import copy
import re

def read_file(infile_name):
    """read all results from one software into a list of list, every index
    of sub list indicate chr
    """
    infile = open(infile_name)
    for line in infile:
        if line.startswith('SL2.40'):
            info_list = line.strip().split()[3].split(';')
            map_qua = re.search(r'MQ=([0-9]*)',info_list[1]).group(1)
            te_length = re.search(r'TE=([0-9]*)',info_list[-1]).group(1)
            sr = re.search(r'SR=(.*)',info_list[0]).group(1)
            sr_total = int(sr.split(',')[0])
            sr_start = int(sr.split(',')[-2])
            sr_end = int(sr.split(',')[-1])
            info_start = int(sr.split(',')[2])
            info_end = int(sr.split(',')[3])
            if sr_total>=3 and (sr_start>=1 and sr_end>=1) or (info_start>=1 and info_end>=1):
                print line
        else:
            pass
    infile.close()
    return 

if  __name__=="__main__":
    read_file(argv[1])
