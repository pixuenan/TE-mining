#!/usr/bin/python
"""
Xuenan Pi
Script to run ITIS
"""

from sys import argv
import subprocess
import os
import time

class Run_itis(object):
    def __init__(self,access_num,te_name,important_data_path,raw_data_path,\
                 refseq_name):
        self.te_seq = "%s/transposon/%s.fasta" % (raw_data_path,te_name)
        self.run_path = "%s/itis/%s/" % (raw_data_path,te_name)
        self.output = "itisrf%s" % access_num
        self.refseq = "%s/ref/%s" % (raw_data_path,refseq_name)
        self.reads1 = '%s/fastq/RF_%s_1.fastq' % (raw_data_path,access_num)
        self.reads2 = '%s/fastq/RF_%s_2.fastq' % (raw_data_path,access_num)
        self.bamfile = '%s/BAM/RF_%s.bam' % (raw_data_path,access_num)
        self.st_raw_file = '%s/store/%s/rf%s/itisrf%s.%s.raw.bed' % \
                (important_data_path,te_name,access_num,access_num,te_name)

    def itis_run(self,access_num,insert_size,refseq_insert_size,ITIS_path):
        #make directory for the TE
        if not os.path.isdir(self.run_path):#/mnt/scratch/pi002/itis/Rider/
            cmd = 'mkdir %s' % self.run_path
            subprocess.check_call(cmd,shell=True)
        #cd /mnt/scratch/pi002/itis to run the program, 
        #otherwise the result cannot be put in a proper place
        os.chdir(self.run_path)#cd /mnt/scratch/pi002/itis/Rider/
        folder_lst = self.check_dir_folder()
        #check if there is result in important data path or raw data path
        if not (os.path.isfile(self.st_raw_file) or self.output in folder_lst):
            if access_num=='000':#reference 
                cmd = ("screen -d -m perl %s/ITIS_v1/itis.pl -g %s -t %s -l"+\
                      " %d -N %s -1 %s -2 %s -e Y") % (ITIS_path,self.refseq,\
                      self.te_seq,refseq_insert_size,self.output,\
                      self.reads1,self.reads2)
            else:
                cmd = ("screen -d -m perl %s/ITIS_v1/itis.pl -g %s -t %s -l"+\
                      " %d -N %s -1 %s -2 %s -F Y -B %s -e Y") % (ITIS_path,\
                      self.refseq,self.te_seq,insert_size,self.output,\
                      self.reads1,self.reads2,self.bamfile)
            subprocess.check_call(cmd,shell=True)

    def check_dir_folder(self):
        """get the name of folders in ITIS results"""
        result = []
        lst = os.listdir(self.run_path)
        for name in lst:
            result += name.split('.')
        return result

if __name__=='__main__':
    important_data_path = "/mnt/nexenta/pi002"
    raw_data_path = "/mnt/scratch/pi002"
    ITIS_path = "/home/pi002/progs_nobackup"
    refseq_name = "S_lycopersicum_chromosomes.2.40.fa"
    insert_size = 500
    refseq_insert_size = 300
    te_name = argv[1]
    access_num_ls = argv[2]#001,011,019
    for access_num in access_num_ls.split(','):
        run = Run_itis(access_num,te_name,important_data_path,raw_data_path,\
                       refseq_name)
        run.itis_run(access_num,insert_size,refseq_insert_size,ITIS_path)
