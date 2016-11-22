#!/usr/bin/python
"""
6-10-2015
Xuenan Pi
Script to extract the region around TE locations in the BAM file and upload
it into JBrowse.
"""
from sys import argv
import subprocess
import os.path
import time

def bam_extracter(bamfile,loci_list,extract_bam,access_num,wait_time,\
                  region_length):
    """parse the BAM file to a truncated file that only contains the 
    region around the TE locations"""
    #index the BAM file e.g. RF_001.bam
    if not os.path.isfile(bamfile+".bam.bai"):
        sort_cmd = "screen -m -d samtools index %s.bam" % bamfile
        subprocess.check_call(sort_cmd,shell=True)
        print "index of %s is running" % bamfile
        #wait this step to finish
        time.sleep(wait_time)
    #creat the truncated bam file only contain the region around the
    #TE locations
    if not os.path.isfile(extract_bam+".bam.bai"):
        bam_cmd = ("cat %s | awk \'{print $1\":\"$2-%d\"-\"$3+%d}\' | "+\
            "xargs samtools view -F 4 -b %s.bam -o %s.bam") % (loci_list,\
                region_length,region_length,bamfile,access_num)
        subprocess.check_call(bam_cmd,shell=True)
        sort_cmd = "samtools sort %s.bam %s" % (access_num,extract_bam)
        subprocess.check_call(sort_cmd,shell=True)
        index_cmd = "samtools index %s.bam" % extract_bam
        subprocess.check_call(index_cmd,shell=True)
        rm_cmd = "rm %s.bam" % access_num 
        subprocess.check_call(rm_cmd,shell=True)

def bam_copy(te_name,access_num,extract_bam,jbrowse_path):
    """copy the truncated bam file to jbrowse folder and write the 
    'tracks.conf' file"""
    jbrowse_BAM_path = '%s/data/BAM/%s' % (jbrowse_path,te_name)
    if not os.path.isdir(jbrowse_path):
        cmd = 'mkdir %s' % jbrowse_path
        subprocess.check_call(cmd,shell=True)
    copied_bam = '%s/rf%s.bam' % (jbrowse_path,access_num)
    copy_cmd = 'cp %s.bam %s.bam.bai %s' % (extract_bam,extract_bam,\
                                            jbrowse_BAM_path)
    subprocess.check_call(copy_cmd,shell=True)
    track_file_name = '%s/data/tracks.conf' % jbrowse_path
    #caution, the same information may write again and again
    track_file = open(track_file_name,'a')
    track_file.write('[tracks.%s_rf%s]\n' % (te_name,access_num))
    track_file.write(('storeClass     = JBrowse/Store/SeqFeature/BAM\n'+\
            'urlTemplate    = BAM/%s/rf%s.bam\nbaiUrlTemplate = BAM/%s/'+\
            'rf%s.bam.bai\n') % (te_name,access_num,te_name,access_num)) 
    track_file.write(('category = NGS/%s\nmaxHeight = 150\ntype = JBrowse'+\
            '/View/Track/Alignments2\nkey = truncated bam for RF_%s\n\n') \
            % (te_name,access_num))
    track_file.close()

if __name__=='__main__':
    important_data_path = "/mnt/nexenta/pi002"
    raw_data_path = "/mnt/scratch/pi002"
    jbrowse_path = '/mnt/geninf15/prog/www/htdocs/transposon/jbrowse'
    wait_time = 1800
    insert_size = 500
    read_length = 100
    region_length = insert_size - read_length
    te_name = argv[1]
    access_num_ls = argv[2]#000,001,011,019
    loci_list = argv[3]
    for access_num in access_num_ls.split(','):
        bamfile = '%s/BAM/RF_%s' % (raw_data_path,access_num)
        extract_bam = '%s/result/%s/rf%s/rf%s' % (important_data_path,\
                te_name,access_num,access_num)
        bam_extracter(bamfile,loci_list,extract_bam,access_num,wait_time,\
                      region_length)
        #wait until the parse step finished
        time.sleep(wait_time/10)
        bam_copy(te_name,access_num,extract_bam,jbrowse_path)
