#!/usr/bin/python
"""
6-9-2015
Xuenan Pi
Script to detect all kinds of TE locations from ITIS results
"""
from sys import argv
import os.path
import os
import subprocess
import time

class Detection(object):
    def __init__(self,access_num,te_name,important_data_path,raw_data_path):
        #e.g. /mnt/nexenta/pi002/store/Rider/rf001
        self.store_path = '%s/store/%s/rf%s' % (important_data_path,te_name,\
                           access_num)
        #e.g. /mnt/nexenta/pi002/result/Rider        
        self.result_path_te = '%s/result/%s' % (important_data_path,te_name)
        #e.g. /mnt/nexenta/pi002/result/Rider/rf001    
        self.result_path = '%s/result/%s/rf%s' % \
                           (important_data_path,te_name,access_num)
        self.itis_raw = '%s/itisrf%s.%s.raw.bed' % \
                        (self.store_path,access_num,te_name)
        self.itis_raw_corr = '%s/itisrf%s.raw.corr.bed' % \
                             (self.result_path,access_num)
        self.filter_file = '%s/itisrf%s.%s.filtered.bed' % \
                           (self.store_path,access_num,te_name)
        self.ref_aln_bam = '%s/itisrf%s.ref_and_te.sorted.bam' % \
                           (self.store_path,access_num)
        self.filter_nbn = '%s/itisrf%s.nbn.bed' % \
                          (self.store_path,access_num)
        self.copy_ori_file = '%s/itisrf%s.copy.bed' % \
                             (self.result_path,access_num)
        self.copy_corr_file = '%s/itisrf%s.copy.corr.bed' % \
                              (self.result_path,access_num)
        self.copy_anno_file = '%s/itisrf%s.copy.anno.bed' % \
                              (self.result_path,access_num)
        self.insertion_ori_file = '%s/itisrf%s.insertion.bed' % \
                                  (self.result_path,access_num)
        self.insertion_anno_file = '%s/itisrf%s.insertion.anno.bed' % \
                                   (self.result_path,access_num)
        self.gap_copy_ori_file = '%s/itisrf%s.gap.copy.bed' % \
                                 (self.result_path,access_num)
        self.gap_copy_corr_file = '%s/itisrf%s.gap.copy.corr.bed' % \
                                  (self.result_path,access_num)
        self.gap_copy_anno_file = '%s/itisrf%s.gap.copy.anno.bed' % \
                                  (self.result_path,access_num)
        self.introgression_file = ('%s/store/introgression/'+\
                'rf%s.median_polish.bed') % (important_data_path,access_num)
        self.deletion_file = '%s/store/deletion/rf%s_deletion.info' % \
                             (important_data_path,access_num)
        self.homo_file = '%s/store/%s/rf000/corred_homo.lst' % \
                         (important_data_path,te_name)
        self.gap_file = '%s/store/deletion/reffa.Ngap.info' % \
                        important_data_path

##detect non-reference insertions##

    def filter_raw(self,tag_script,te_name,blast_corr_script,insert_size):
        """find non-reference insertion from itis raw result"""
        if not os.path.isfile(self.insertion_ori_file):
            #if no /mnt/nexenta/pi002/result/Rider    
            if not os.path.isdir(self.result_path_te):
                cmd = 'mkdir %s' % self.result_path_te
                subprocess.check_call(cmd,shell=True)
            #mkdir /mnt/nexenta/pi002/store/Rider/rf001        
            cmd = 'mkdir %s' % self.result_path
            subprocess.check_call(cmd,shell=True)
            if not os.path.isfile(self.itis_raw_corr):
                cmd = 'grep -v MQ=0 %s > %s' % (self.itis_raw,\
                    self.itis_raw_corr)
                subprocess.check_call(cmd,shell=True)
            #new insertions from itis filtered result
            if not os.path.isfile(self.filter_nbn):
                self.run_nbn()
            if not os.path.isfile(self.homo_file):#modified blast result
                self.modify_blast_result(te_name,blast_corr_script)
            cmd = 'python %s %s %s %s %s %s %d' % \
                (tag_script,self.itis_raw_corr,self.filter_nbn,self.homo_file,\
                self.gap_file,self.insertion_ori_file,insert_size)
            subprocess.check_call(cmd,shell=True)

    def run_nbn(self):
        """extract new insertions from itis filtered result"""
        cmd = "grep NB=N %s > %s" % (self.filter_file,self.filter_nbn)
        subprocess.call(cmd,shell=True)

    def modify_blast_result(self,te_name,blast_corr_script):
        """modify the blast result"""
        blast_result = '%s/rf000/te_homo_in_ref.lst' % self.result_path_te
        sorted_file = '%s/rf000/sorted_homo.lst' % self.result_path_te
        cmd = 'sort -k1,1 -k2n %s > %s' % (blast_result,sorted_file)
        subprocess.check_call(cmd,shell=True)
        cmd = 'python %s %s %s' % (blast_corr_script,sorted_file,\
                                   self.homo_file)
        subprocess.check_call(cmd,shell=True)

##detect non-reference insertions##

    def detect_copy(self,tag_script,te_name,min_len,max_len,insert_size,\
                    read_length,wait_time):
        """output copies by combine blast result and ITIS result"""
        if not os.path.isfile(self.copy_ori_file):
            cmd = 'screen -m -d python %s %s %s %s %s %s %s %d %d'  % \
                  (tag_script,self.ref_aln_bam,self.homo_file,self.gap_file,\
                  self.copy_ori_file,self.gap_copy_ori_file,te_name,\
                  insert_size,read_length)
            subprocess.check_call(cmd,shell=True)
            time.sleep(wait_time)
        if os.path.isfile(self.copy_ori_file) and \
           os.path.isfile(self.gap_copy_ori_file):
            cmd_corr = "awk '$4>=%d && $4<=%d ' %s > %s" % \
            (min_len,max_len,self.copy_ori_file,self.copy_corr_file)
            subprocess.check_call(cmd_corr,shell=True)
            cmd_corr = "awk '$4>=%d && $4<=%d ' %s > %s" % \
            (min_len,max_len,self.gap_copy_ori_file,self.gap_copy_corr_file)
            subprocess.check_call(cmd_corr,shell=True)

    def deletion(self,deletion_script,access_num,insert_size,read_length,\
                 raw_data_path,wait_time):
        if not os.path.isfile(self.deletion_file):
            print 'Creating deletion file for RF%s...' % access_num
            bam_index = '%s/BAM/RF_%s.bam.bai' % (raw_data_path,access_num)
            if not os.path.isfile(bam_index):
                index_cmd = 'screen -d -m samtools index %s/BAM/RF_%s.bam' % \
                            (raw_data_path,access_num)
                subprocess.check_call(index_cmd,shell=True)
                time.sleep(wait_time)
                cmd = 'screen -d -m python %s %s/BAM/RF_%s.bam %s %d %d' % \
                      (deletion_script,raw_data_path,access_num,\
                      self.deletion_file,insert_size,read_length)
                subprocess.check_call(cmd,shell=True)
            else:
                cmd = 'screen -d -m python %s %s/BAM/RF_%s.bam %s %d %d' % \
                      (deletion_script,raw_data_path,access_num,\
                      self.deletion_file,insert_size,read_length)
                subprocess.check_call(cmd,shell=True)

    def annotate_introgression(self,introgression_anno_script,access_num):
        if os.path.isfile(self.copy_corr_file) and \
           os.path.isfile(self.gap_copy_corr_file) and \
           os.path.isfile(self.insertion_ori_file):
            #annotate all accessions except reference genome and wild species 
            #with introgression
            if not access_num in ['000','073']:
                cmd = 'python %s %s %s > %s' % \
                (introgression_anno_script,self.copy_corr_file,\
                self.introgression_file,self.copy_anno_file)
                subprocess.check_call(cmd,shell=True)
                cmd = 'python %s %s %s > %s' % \
                (introgression_anno_script,self.gap_copy_corr_file,\
                self.introgression_file,self.gap_copy_anno_file)
                subprocess.check_call(cmd,shell=True)
                cmd = 'python %s %s %s > %s' % \
                (introgression_anno_script,self.insertion_ori_file,\
                self.introgression_file,self.insertion_anno_file)
                subprocess.check_call(cmd,shell=True)
            else:
                awk_cmd = r'''awk {'printf("%s\t%s\t%s\t%s\n",$1,$2,$3,"ref")'}'''
                cmd =  '%s %s > %s' % \
                (awk_cmd,self.copy_corr_file,self.copy_anno_file)
                subprocess.check_call(cmd,shell=True)
                cmd =  '%s %s > %s' % \
                (awk_cmd,self.gap_copy_corr_file,self.gap_copy_anno_file)
                subprocess.check_call(cmd,shell=True)
                cmd =  '%s %s > %s' % \
                (awk_cmd,self.insertion_ori_file,self.insertion_anno_file)
                subprocess.check_call(cmd,shell=True)


if __name__=='__main__':
    important_data_path = "/mnt/nexenta/pi002"
    raw_data_path = "/mnt/scratch/pi002"
    insert_size = 500
    read_length = 100
    wait_time = 1800
    ref_insertion_script = important_data_path+\
                           '/script/ref_insertion_finder.py'
    nonref_insertion_script = important_data_path+\
                              '/script/nonref_insertion_finder.py'
    copy_data_script = important_data_path+'/script/copy_rawdata.py'
    deletion_script = important_data_path+'/script/large_deletion_finder.py'
    introgression_anno_script = important_data_path+\
                                '/script/introgression_annotation.py'
    blast_corr_script = important_data_path+'/script/blastcorr.py'
    te_name = argv[1]
    access_num_ls = argv[2]#000,001,011,019
    min_len,max_len = int(argv[3]),int(argv[4])
    min_len = min_len + min_len*0.1
    cmd = 'python %s %s %s %s %s' % (copy_data_script,te_name,access_num_ls,\
                                     important_data_path,raw_data_path)
    subprocess.check_call(cmd,shell=True)
    for access_num in access_num_ls.split(','):
        detect = Detection(access_num,te_name,important_data_path,\
                           raw_data_path)
        detect.filter_raw(nonref_insertion_script,te_name,blast_corr_script,\
                          insert_size)
        detect.detect_copy(ref_insertion_script,te_name,min_len,max_len,\
                           insert_size,read_length,wait_time/10)
        detect.annotate_introgression(introgression_anno_script,access_num)
        if access_num!='000':
            detect.deletion(deletion_script,access_num,insert_size,\
                            read_length,raw_data_path,wait_time)
    
