#!/usr/bin/python
"""
Xuenan Pi
18-11-2015
Script to annotate TE location with introgression 
"""

from sys import argv
import copy

def read_file(infile_name):
    """read all results from one software into a list of list, every index
    of sub list indicate chr
    """
    chr_dict = {}
    infile = open(infile_name)
    for line in infile:
        if line.strip() != '\n':
            chr_name = line.strip().split()[0]
            loci = map(int,line.strip().split()[1:3])
            if chr_name in chr_dict:
                chr_dict[chr_name] += [loci+['']]
            else:
                chr_dict[chr_name] = [loci+['']]
        else:
            pass
    infile.close()
    return chr_dict

def write_file(chr_dict):
    """Write the solo insertions to a file.
    """
    chr_name_list = sorted(chr_dict.keys())
    for chr_name in chr_name_list:
        chr_list = chr_dict[chr_name]
        if chr_list:
            for loc in chr_list:
                print "%s\t%d\t%d\t%s" % (chr_name,loc[0],loc[1],loc[2])
        
def compare_2_dicts(dict1,dict2):
    """Extract the overlap between two locations.
    """
    for chr_name,chr_list in dict1.items():
        for k,loc1 in enumerate(chr_list):
            flag = True
            if chr_name in dict2:
                for loc2 in dict2[chr_name]:
                    #check whether the insertion is inside of an introgression
                    #segement
                    if loc1[0]>=loc2[0] and loc1[1]<=loc2[1]:
                        loc1[2] = 'introgression'
                        flag = False
                        break
            if flag:
                loc1[2] = 'nonIG'
    return dict1

if  __name__=="__main__":
    file1 = argv[1]#the location file
    file2 = argv[2]#the introgression annotation file
    dict1 = read_file(argv[1])
    dict2 = read_file(argv[2])
    common = compare_2_dicts(dict1,dict2)
    write_file(common)
   

    
