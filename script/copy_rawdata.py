#!/usr/bin/python
"""
17-9
Xuenan Pi
Copy the important results from ITIS in raw_data_path to 
important_data_path for storage
"""

from sys import argv
import subprocess
import os
import time

def copy_file(access_num,te_name,important_data_path,raw_data_path):
    """copy neccessary data from raw_data_path to important_data_path"""
    te_seq = te_name+'.fasta'
    run_path = "%s/itis/%s/" % (raw_data_path,te_name)
    store_path = '%s/store/%s/' % (important_data_path,te_name)
    output = 'itisrf%s' % access_num
    st_raw_file = '%srf%s/%s.%s.raw.bed' % (store_path,access_num,output,\
                                            te_name)
    if not os.path.isfile(st_raw_file):
        # cd ITIS result folder
        os.chdir(run_path)
        folder_lst = check_dir_folder(run_path)
        # make all directories if they are not existed yet
        if output in folder_lst:
            folder_inx = folder_lst.index(output)
            folder_nam = '.'.join(folder_lst[folder_inx:folder_inx+2])
            if not os.path.isdir(store_path):
                cmd = 'mkdir %s' % store_path
                subprocess.check_call(cmd,shell=True)
            cmd = 'mkdir %srf%s' % (store_path,access_num)
            subprocess.check_call(cmd,shell=True)
            # copy blast result of reference genome
            if access_num=='000':
                blast_result = '%s/te_homo_in_ref.lst' % folder_nam
                cmd = 'cp %s %srf%s' % (blast_result,store_path,access_num)
                subprocess.check_call(cmd,shell=True)
            # copy the .raw, .filtered and shrinked bam file from ITIS
            rn_raw_file = '%s/%s.%s.raw.bed' % (folder_nam,output,te_name)  
            rn_fil_file = '%s/%s.%s.filtered.bed' % (folder_nam,output,\
                                                     te_name)
            ref_aln_bam = '%s/%s.ref_and_te.sorted.bam' % (folder_nam,\
                                                           output)
            cmd = 'cp %s %s %s %srf%s' % (rn_raw_file,rn_fil_file,\
                                          ref_aln_bam,store_path,access_num)
            subprocess.check_call(cmd,shell=True)
            #rm the itis result
            #cmd = 'rm %s' % output
            #subprocess.check_call(cmd,shell=True)

def check_dir_folder(run_path):
    """get the name of folders in ITIS results"""
    result = []
    lst = os.listdir(run_path)
    for name in lst:
        result += name.split('.')
    return result

if __name__=='__main__':
    te_name = argv[1]
    access_num_ls = argv[2]#001,011,019
    important_data_path = argv[3]
    raw_data_path = argv[4] 
    for access_num in access_num_ls.split(','):
        copy_file(access_num,te_name,important_data_path,raw_data_path)
