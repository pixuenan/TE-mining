#!/usr/bin/python
"""
9-9-2015
Xuenan Pi
Script to make gene annotation based on ITAG annotation. The script will 
firstly annotate whether the insertion is inside of a gene or not, if so which
part, its distance to nearest gene and gene id. Then the function of the gene 
will be added corresponde to gene id. The extra annotation file means the file
has different annotation content to the same gene id.
"""

from sys import argv
import os
import subprocess

def call_closest(infile,bedtools2_path,annotation_file):
    """call BEDTools closest"""
    os.chdir(bedtools2_path)
    cmd = './bin/closestBed -D \'b\' -t last -a %s -b %s' % \
          (infile,annotation_file)
    p = subprocess.Popen([cmd],shell=True,stdout=subprocess.PIPE,\
        stderr=subprocess.PIPE)
    outfile,err = p.communicate()
    return outfile

def read_file(infile):
    """read the output file from BEDTools closest and extract 
    needed information"""
    line_list = []
    for line in infile.strip().split('\n'):
        if line:
            info_list = line.strip().split('\t')
            insert_part = info_list[-8]
            strand = info_list[-4]
            gene_info = info_list[-2].split(';')
            distan = int(info_list[-1])
            if distan != 0:
                insert_part = 'intergenic'
                ID = gene_info[1].split(':')[1]
                line_list += [info_list[:-10]+[insert_part,strand,distan,ID]]
            else:
                ID = gene_info[1].split(':')[1]
                line_list += [info_list[:-10]+[insert_part,strand,distan,ID]]
    return line_list

def read_function(infile_name):
    """read the truncated ITAG2.3 gene model file"""
    result = {}
    infile = open(infile_name)
    for line in infile:
        if line.strip() != '\n':
            info = line.strip().split('\t')[-1]
            info_list = info.split(';')
            ID = info_list[0].split(':')[1]
            func = info_list[2].split('=')[1]
            result[ID] = func
    infile.close()
    return result

def read_function2(infile_name):
    """read extra gene annotation file"""
    result = {}
    infile = open(infile_name)
    for line in infile:
        if line.strip() != '\n':
            info = line.strip().split(',')
            ID = info[0]
            func = info[1:]
            result[ID] = func
    infile.close()
    return result

def anno_func(lines,func_dict2,access_ls,header,func_dict1=None):
    """add gene function annotation"""
    print header
    for line in lines:
        ID1 = line[-1].strip().split('.')[0]#e.g. Solyc00g005000
        ID = line[-1]#e.g. Solyc00g005000.2.1
        if func_dict1:
            if ID1 in func_dict1.keys():
                print '\t'.join(map(str,line+func_dict1[ID1]+[func_dict2[ID]]))
            else:
                print '\t'.join(map(str,line+[func_dict2[ID]]))
        else:
            print '\t'.join(map(str,line+[func_dict2[ID]]))


if __name__=='__main__':
    infile = argv[1]
    access_ls = argv[2]#000,001,011,019
    important_data_path = argv[3]
    flag = argv[4]
    bedtools2_path = '%s/store/annotation/bedtools2/' % important_data_path
    annotation_file = '%s/store/annotation/ITAG2.3_gene_models_corr.gff3' % \
                      important_data_path
    outfile = call_closest(infile,bedtools2_path,annotation_file)
    lines = read_file(outfile)
    #extra gene annotation
    if flag == 'extra':
        #func_file1 is the extra function annotation file
        func_file1 = '%s/store/annotation/Ruud_gene_list_annotated.csv' % \
                 important_data_path
        func_dict1 = read_function2(func_file1)
        header = ("chr\tstart\tend\ttype\tlength\trf%s\tshared_num\tinsert\t"+\
                "strand\tdistance\tgene_id\tRuud\tlocus_symbol\t"+\
                "locus_name/description\tRemarks\tITAG_2.30_gene_model") % \
                '\trf'.join(access_ls.split(','))
    else:
        func_dict1 = None
        header = ("chr\tstart\tend\ttype\tlength\trf%s\tshared_num\tinsert\t"+\
                "strand\tdistance\tgene_id\t"+\
                "ITAG_2.30_gene_model") % '\trf'.join(access_ls.split(','))
    #ITAG2.30 gene annotation
    #func_file2 is the standard function annotation file
    func_file2 = '%s/store/annotation/ITAG2.3_gene_function.gff3' % \
                 important_data_path
    func_dict2 = read_function(func_file2)
    anno_func(lines,func_dict2,access_ls,header,func_dict1)
    
