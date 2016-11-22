#!/usr/bin/python
"""
Xuenan Pi
31-8-2015
Script to run pipeline for transposon mining in tomato genome project
"""

from sys import argv
import subprocess
import os
import os.path

def split_genome(genome_list):
    """split the genome into groups of 5 members"""
    #write the command line into a bash file to execute them at the same time on screen
    #ITIS will run on 5 genomes every time
    genome_num = len(genome_list)
    for i in range(genome_num/5+1):
        if 5*(i+1) > genome_num: 
            if genome_num%5 == 0:
                continue
            else:
                infile_list = genome_list[i*5:]
        else:
            infile_list = genome_list[i*5:(i+1)*5]
        yield infile_list

if  __name__=="__main__":
    genome_list = argv[1].split(',')#input should be a list of genomes, like 01,11,19
    for i in split_genome(genome_list):
        print i
