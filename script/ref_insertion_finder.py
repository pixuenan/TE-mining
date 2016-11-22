#!/usr/bin/python
"""
Xuenan Pi
29-9-2015
Script to detect reference insertions
"""

from sys import argv
from functions import read_file_1
import os.path
import os
import subprocess

def extract_reads(infile_name,te_name,mq=30):
    """use samtool to extract reads mapped around the end of the transposon 
    sequence
    """ 
    outfile = None
    if os.path.isfile(infile_name):
        cmd_f = "samtools view -F 16 %s | awk '$3!=\"%s\" && $7!=\"=\" && \
                $5 >= %d'" % (infile_name,te_name,mq)
        p_f = subprocess.Popen([cmd_f],shell=True,stdout=subprocess.PIPE,\
                               stderr=subprocess.PIPE)
        forward,err = p_f.communicate()
        cmd_r = "samtools view -f 16 %s | awk '$3!=\"%s\" && $7!=\"=\" && \
                $5 >= %d'" % (infile_name,te_name,mq)
        p_r = subprocess.Popen([cmd_r],shell=True,stdout=subprocess.PIPE,\
                               stderr=subprocess.PIPE)
        reverse,err = p_r.communicate()
    return forward,reverse

def group(chr_prefix,infile,insert_size,read_length,flag):
    """put all supportive reads in a group"""
    final_result = []
    fileline = infile.strip().split('\n')
    p_list = [[0,0,0,0]] 
    i = 0
    while i < len(fileline):
        line = fileline[i].strip().split('\t')[:10]
        if flag:
            line += ['-']
        else:
            line += ['+']
        if chr_prefix in line[2]:
            result = check_line(p_list,line,insert_size,read_length)
            if result[0]:
                p_list += [result[1:]]
                i += 1
            else:
                if not i+1 == len(fileline):
                    final_result += [p_list]
                    p_list = [result[1:]]
                    i += 1
                else:#only for the last read
                    final_result += [p_list]
                    p_list = [result[1:]]
                    i += 1
        else:
            i += 1
            pass
    return final_result

def check_line(p_list,line,insert_size,read_length):
    """check if the reads belong to the group"""
    line_list = line[:]
    chr = line_list[2]
    loci1 = int(line_list[3])#start of alignment
    if chr == p_list[0][0] and abs(loci1-p_list[0][1]) <= insert_size:
        return True,chr,loci1,len(line[-2]),line[-1]
    else:
        return False,chr,loci1,len(line[-2]),line[-1]


def info(raw_list,read_length):
    """"generate the start and end location for each group. Like chr start 
    end num_sr length: SL2.40ch00 1200 1300 10 100."""
    for i in raw_list:
        if len(i) >= 3:#filter out sr<=3  
            chr = i[0][0]
            if i[0][-1]=='-':#left alignment
                end = i[-1][1]+read_length
                leftsr = len(i)
                yield [chr,end,leftsr,"left"]
            else:
                start = i[0][1]
                rightsr = len(i)
                yield [chr,start,rightsr,"right"]

def write_file(result_list,homo_dict,gap_dict,distan,outfile_copy,outfile_gap):
    """confirm the group of face each other reads are matched with blast result
    """
    i = 0
    out_gap = open(outfile_gap,'w+')
    out_copy = open(outfile_copy,'w+')
    while i<len(result_list)-1:
        if result[i][-1]=='left' and result[i+1][-1]=='right':
            left = result[i]
            right = result[i+1]
            chr_name = left[0]
            start = left[1]
            end = right[1]
            homolist = homo_dict[chr_name]
            gaplist = gap_dict[chr_name]
            homo_flag,homo_start_distan,homo_end_distan = check_loci(start,\
                                                        end,homolist,distan)
            #location is confirmed, extract the exact length by blast result
            if homo_flag:
                for homo_loci in homolist:
                    if abs(start-homo_loci[0])==homo_start_distan:
                        homo_start=homo_loci[0]
                    if abs(end-homo_loci[1])==homo_end_distan:
                        homo_end=homo_loci[1]
                length = homo_end-homo_start
                out_copy.write('\t'.join([result_list[i][0]]+map(str,\
                               [homo_start,homo_end,length])+['\n']))
                i += 2
            else:
                gap_flag,gap_start_distan,gap_end_distan = check_loci(start,\
                                                           end,gaplist,distan)
                length = end-start
                flag = False
                if gap_flag:
                    for gap_loci in gaplist:
                        if abs(start-gap_loci[0])==gap_start_distan and \
                                abs(end-gap_loci[1])==gap_end_distan:
                            gap_start=gap_loci[0]
                            gap_end=gap_loci[1]
                            out_gap.write('\t'.join([result_list[i][0]]+map(\
                            str,[gap_start,gap_end,gap_end-gap_start])+['\n']))
                            flag = True
                            break
                if flag:
                    i += 2
                else:
                    i += 1
        else:
            i += 1
    out_gap.close()
    out_copy.close()

def check_loci(start,end,lst,distan):
    """check if the location has minimal distance with a location in a list"""
    flag,min_start_distan,min_end_distan = False,None,None
    if lst:
        min_start_distan = min([abs(start-loci[0]) for loci in lst])
        min_end_distan = min([abs(end-loci[1]) for loci in lst])
        if min_start_distan<=distan and min_end_distan<=distan:
            flag = True
    return flag,min_start_distan,min_end_distan

if  __name__=="__main__":
    bam_file,homofile,gapfile = argv[1],argv[2],argv[3]
    out_copy,out_gap = argv[4],argv[5]
    te_name = argv[6]
    insert_size = int(argv[7])
    read_length = int(argv[8])
    max_distance = insert_size - read_length
    homodict = read_file_1(homofile)
    gapdict = read_file_1(gapfile)
    chr_prefix = os.path.commonprefix(homodict.keys())
    f_read,r_read = extract_reads(bam_file,te_name)
    group_list = group(chr_prefix,f_read,insert_size,read_length,flag=1)+group\
                      (chr_prefix,r_read,insert_size,read_length,flag=0)
    group_list.sort()
    result = []
    for i in info(group_list,read_length):
        result += [i]
    result.sort()
    write_file(result,homodict,gapdict,max_distance,out_copy,out_gap)

