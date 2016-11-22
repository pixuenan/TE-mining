#!/usr/bin/python
"""
Xuenan Pi
1-9-2015
Script to detect the location of the assembly gap (NNNs) in a fasta file. 
"""

from sys import argv

def Ngap_in_line(lines):
    """Detect Ngap for each line in the fasta file, and add their locations
    to every chromosome. Every line in the fasta file is assumed to have same
    sequence length, except the header lines and last line of each chromosome.
    
    Result
    chr_lists: a list contains sublists, every sublist is corrseponding to a
    chromosome and contains Ngap location in that chromosome
    """
    chr_dict = {}
    for line in lines:
        if line.strip()=='\n':
            pass
        #creat a new list for every chromosome
        elif line.startswith('>'):
            chr_idx = 1
            chr_name = line.strip()[1:]
            chr_dict[chr_name] = []
        #skip the line without NNNs
        elif not 'N' in line.strip():
            chr_idx += len(line.strip())
        else:
            #find the Ngap locations in the line
            ngap_locs = findN(translate(line.strip()),[])
            #convert the Ngap locations in the lines to their location in the
            #chromosome
            for loc in ngap_locs:
                chr_dict[chr_name] += [[chr_idx+loc[0],chr_idx+loc[1]]]
            chr_idx += len(line.strip())
    return chr_dict

def translate(dna_seq):
    """translate the nucleotide sequence to binary sequence, 
    N will be turned into 1, ATCG will be turned into 0"""
    bi_seq = ['0']*len(dna_seq)
    for i,ch in enumerate(dna_seq):
        if ch=='N':
            bi_seq[i] = '1'
    return ''.join(bi_seq)

def findN(bi_seq,result=[],begin=0):
    """find the location of the 1 within a binary sequence recursively"""
    #if there is all 0 in the remained seq
    if not '1' in bi_seq:
        return result
    #if there are all 1 in the remained seq
    elif not '0' in bi_seq:
        result += [[begin+1,begin+len(bi_seq)]]
        return result
    else:
        #find the first 1
        start = bi_seq.index('1')
        if '0' in bi_seq[start:]:
            #end until find 0
            end = bi_seq[start:].index('0')+start
            result += [[start+begin+1,end+begin]]
            remained_seq = bi_seq[end:]
            return findN(remained_seq,result,end)
        else:
            end = len(bi_seq)
            result += [[start+begin+1,end+begin]]
            return result

def link(ngap_list):
    """join the Ngap stop at the end of one line in fasta file and start at
    the begining of the next line together"""
    idx = 0
    while idx+1 < len(ngap_list):
        #if the location of the Ngap is adjacent to the next one, extend
        #this location
        start,end = ngap_list[idx][0],ngap_list[idx][1]
        start_next,end_next = ngap_list[idx+1][0],ngap_list[idx+1][1]
        if end+1==start_next:
            ngap_list[idx] = [start,end_next]
            ngap_list.remove(ngap_list[idx+1])
        else:
            idx += 1
    return ngap_list

if __name__=='__main__':
    infile = open(argv[1])#the input is the reference.fasta
    Ngap_dict = Ngap_in_line(infile)
    infile.close()
    chr_name_list = sorted(Ngap_dict.keys())
    for chr_name in chr_name_list:
        for ngap_loc in link(Ngap_dict[chr_name]):
            #print: name of chromosome, start location of Ngap, end location 
            #of N gap, length of Ngap
            print "%s\t%d\t%d\t%d" % (chr_name,ngap_loc[0],\
                    ngap_loc[1],ngap_loc[1]-ngap_loc[0])

