#!/usr/bin/python
"""
Xuenan Pi
27-7-2015
Script to compare TE-locate,Tif,ITIS result for one accession
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

def write_common_file(outfile_name,list):
    """Write the common insertions to a file.
    """
    chr_name_list = ['SL2.40ch00','SL2.40ch01','SL2.40ch02','SL2.40ch03','SL2.40ch04','SL2.40ch05','SL2.40ch06','SL2.40ch07','SL2.40ch08','SL2.40ch09','SL2.40ch10','SL2.40ch11','SL2.40ch12']
    outfile = open(outfile_name,'w+')
    for index,chr_list in enumerate(list):
        if chr_list:
            chr = chr_name_list[index]
            for loci in chr_list:
                outfile.write("%s\t%d\t%d\n" % (chr,loci[0],loci[1]))
    outfile.close()

def write_solo_file(outfile_name,list):
    """Write the solo insertions to a file.
    """
    chr_name_list = ['SL2.40ch00','SL2.40ch01','SL2.40ch02','SL2.40ch03','SL2.40ch04','SL2.40ch05','SL2.40ch06','SL2.40ch07','SL2.40ch08','SL2.40ch09','SL2.40ch10','SL2.40ch11','SL2.40ch12']
    outfile = open(outfile_name,'w+')
    for index,chr_list in enumerate(list):
        if chr_list:
            chr = chr_name_list[index]
            for loci in chr_list:
                outfile.write("%s\t%d\n" % (chr,loci))
    outfile.close()


def check_smallest(list):
    """Check the smallest distance between the insertion
    """
    distan = 20000000
    for i in range(len(list)):
        if list[i]:
            for j in range(len(list[i])-1):
                if list[i][j+1]-list[i][j] < distan:
                    distan = list[i][j+1]-list[i][j]
    return distan
        
def compare_2_lists(list1,list2,distance):
    """Extract the overlape between two programs.
    """
    #for each chr, compare the overlap
    chr_list = [0]*13 
    for i in range(len(chr_list)):
        chr_list[i] = [] 
    list1_solo = copy.deepcopy(list1)
    list2_solo = copy.deepcopy(list2)
    overlap = copy.deepcopy(chr_list)
    for i in range(len(list1)):
        j = 0
        if list1[i] and list2[i]:
            for loci1 in list1[i]:
                min_distan = min([abs(loci1-loci2) for loci2 in list2[i]])
                flag = False
                if min_distan <= distance:
                    k = list1[i].index(loci1)
                    for j in range(len(list2[i])):
                        if abs(loci1 - list2[i][j]) == min_distan and min_distan < abs(list2[i][j]-list1[i][k-1]):
                            if k < len(list1[i])-1:
                                if min_distan < abs(list2[i][j]-list1[i][k+1]):
                                    flag = True
                            else:
                                flag = True
                        if flag:
                            loci2 = list2[i][j]
                            overlap[i] += [loci1]
                            list1_solo[i].remove(loci1)
                            list2_solo[i].remove(loci2)
                            break
    return list1_solo,list2_solo,overlap
   
def number(list):
    num = 0
    for i in list:
        if i:
            num += len(i)
    return num

if  __name__=="__main__":
    file1 = argv[1]
    file2 = argv[2]
    distance = int(argv[3])
    list1 = read_file(argv[1])
    list2 = read_file(argv[2])
    solo1,solo2,common = compare_2_lists(list1,list2,distance)
    flag = argv[4]
    if flag == 'C':
        outfile_name = argv[5]
        write_solo_file(outfile_name,common)
    elif flag == 'S':
        outfile_name1 = argv[5]
        outfile_name2 = argv[6]
        write_solo_file(outfile_name1,solo1)
        write_solo_file(outfile_name2,solo2)
    elif flag == 'PC':
        print common
    else:
        print number(list1),number(list2),"common number: %d" % number(common)
        #print common
    #itis_file = argv[3]
    #itis_list = read_file(argv[3])
    #print itis_list
    list1 = [[1,4,7],[4,9,13],[]]
    list2 = [[1,2],[],[4,5,6]]
    list3 = [[10,40,60,90]]
    list4 = [[10,20,30,40,50,60,70,80,90]]
    #solo1,solo2,common = compare_2_lists(te_list,tif_list)
    #solo1,solo2,common = compare_2_lists(list4,list3)
    #print common
   

    
