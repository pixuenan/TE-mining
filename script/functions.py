#!/usr/bin/python
"""
Xuenan Pi
31-12-2015
Script includes the frequently used functions in the pipeline.
"""
def read_file_1(infile_name):
    """read all results from file into a dictionary of lists"""
    chr_dict = {}
    infile = open(infile_name)
    for line in infile:
        if line.strip() != '\n':
            chr_name = line.strip().split()[0]
            loci = map(int,line.strip().split()[1:3])
            if chr_name in chr_dict:
                chr_dict[chr_name] += [loci]
            else:
                chr_dict[chr_name] = [loci]
        else:
            pass
    infile.close()
    return chr_dict

def read_file_2(infile_name):
    """read all results into a list of list, every index
    of sub list indicate chr
    """
    chr_dict = {}
    infile = open(infile_name)
    for line in infile:
        if line.strip() != '\n':
            info = line.strip().split('\t')
            chr_name = info[0]
            if chr_name in chr_dict:
                chr_dict[chr_name].append(map(int,info[1:3])+[info[3]])
            else:
                chr_dict[chr_name] = [map(int,info[1:3])+[info[3]]]
        else:
            pass
    infile.close()
    return chr_dict

def initial_list(length):
    ini_list = [0]*length
    for i in range(length):
        ini_list[i] = [] 
    return ini_list
