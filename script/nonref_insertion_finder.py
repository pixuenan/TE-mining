#!/usr/bin/python
"""
Xuenan Pi
27-8-2015
Script to detect non-reference insertions from the .raw file from ITIS result
"""

from sys import argv
from functions import read_file_1
import copy
import re

def write_file(infile_name,homo_dict,filter_dict,gap_dict,distan,\
               outfile_ins):
    """tag the result from itis with insertion
    """
    infile = open(infile_name)
    out_ins = open(outfile_ins,'w+')
    flag = 0
    for line in infile:
        if line.strip() != '\n':
            line_list = line.strip().split()
            chr_name = line_list[0]
            start = int(line_list[1])
            end = int(line_list[2])
            homolist = homo_dict.get(chr_name,[])
            filterlist = filter_dict.get(chr_name,[])
            gaplist = gap_dict.get(chr_name,[])
            #filter out the results around blast results
            if check_loci(start,end,homolist,distan):
                pass 
            #filter out the results around assembly gaps
            elif check_loci(start,end,gaplist,distan):
                pass
            elif [start,end] in filterlist:
                split_reads_num = sum(map(int,line_list[3].split(';')[0].\
                                          split(',')[2:4]))
                sr = int(re.search(r'SR=([0-9]*)',line_list[3]).group(1))
                map_qua = int(re.search(r'MQ=([0-9]*)',line_list[3]).\
                              group(1))
                if map_qua >=30 and (split_reads_num>=1 or sr >5):
                    out_ins.write(line)
            else:
                map_qua = int(re.search(r'MQ=([0-9]*)',line_list[3]).\
                              group(1))
                sr = int(re.search(r'SR=([0-9]*)',line_list[3]).group(1))
                if map_qua >= 30 and sr >= 30:
                    tsd = re.search(r'TS=([0-9]*|NA)',line_list[3]).\
                                    group(1)
                    if sum(map(int,line_list[3].split(';')[0].split(',')\
                               [2:4])):#the number of split reads >=1
                        out_ins.write(line)
                else:
                    pass
    infile.close()
    out_ins.close()
    return 

def check_loci(start,end,lst,distan):
    """check if the location satify the criteria"""
    flag = False
    if lst:
        min_start_distan = min([abs(start-loci[0]) for loci in lst])
        min_end_distan = min([abs(end-loci[1]) for loci in lst])
        if min_start_distan<=distan or min_end_distan<=distan:
            flag = True
    return flag

if  __name__=="__main__":
    rawfile = argv[1]#the .raw file
    filterlist = read_file_1(argv[2])#read the .filter file 
    homolist = read_file_1(argv[3])#read the blast file
    gaplist = read_file_1(argv[4])#read the gap file
    out_ins = argv[5]
    insert_size = int(argv[6])
    max_distance = 2*insert_size
    write_file(rawfile,homolist,filterlist,gaplist,max_distance,out_ins)
