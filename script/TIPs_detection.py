#!/usr/bin/python
"""
Xuenan Pi
14-9-2015
Script to detect TIPs.
"""

from sys import argv
import subprocess
import operator
import csv

def call_TIPs(important_data_path,te_name,access_num_ls,insert_size):
    cmd_ins = "python %s/script/nonref_TIPs_detection.py %s %s %s %d" % \
            (important_data_path,te_name,access_num_ls[4:],important_data_path\
            ,insert_size)
    p_ins = subprocess.Popen([cmd_ins],shell=True,stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE)
    outfile_ins,err = p_ins.communicate()
    cmd_copy = "python %s/script/ref_TIPs_detection.py %s %s copy %s %d" % \
            (important_data_path,te_name,access_num_ls,important_data_path,\
            insert_size)
    p_copy = subprocess.Popen([cmd_copy],shell=True,stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE)
    outfile_copy,err = p_copy.communicate()
    cmd_gap = "python %s/script/ref_TIPs_detection.py %s %s gap %s %d" % \
            (important_data_path,te_name,access_num_ls,important_data_path,\
            insert_size)
    p_gap = subprocess.Popen([cmd_gap],shell=True,stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE)
    outfile_gap,err = p_gap.communicate()
    outfile_combine = outfile_ins+outfile_copy+outfile_gap
    return outfile_combine

def gene_annotation(outfile_content,important_data_path,access_num_ls):
    out = open("tmpout",'w+')
    out.write(outfile_content)
    out.close()
    cmd = ("sort -k1,1 -k2n tmpout | python %s/script/gene_annotation.py"+\
          " stdin %s %s extra > %s") % (important_data_path,access_num_ls,\
                important_data_path,outfile)
    subprocess.check_call(cmd,shell=True)
    cmd_rm = 'rm tmpout' 
    subprocess.check_call(cmd_rm,shell=True)

if __name__=='__main__':
    important_data_path = "/mnt/nexenta/pi002"
    te_name = argv[1]
    access_num_ls = argv[2]#000,001,011,019
    insert_size = 500
    outfile = "%s_%dgn.tsv" % (argv[1],argv[2].count(',')+1)
    out_content = call_TIPs(important_data_path,te_name,access_num_ls,\
                  insert_size)
    gene_annotation(out_content,important_data_path,access_num_ls)
