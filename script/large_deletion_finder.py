#!/usr/bin/python
"""
Xuenan Pi
2-9-2015
Script to detect deletions, the deletion has been filtered with max length
of 20kb and reads mapping quality above 30
"""

from sys import argv
import os.path
import os
import subprocess

def extract_reads(infile_name,pair_reads_distance):
    """use samtool to extract paired reads with abnormal insert size, and
    use pipeline command to choose the paired reads on same chromosome, 
    with insert size larger than threshold and mapping quality above 30"""
    outfile = None
    if os.path.isfile(infile_name):
        cmd = ("samtools view -f1 -F14 %s | awk '$7 == \"=\" && $5 >= 30'"+\
               "|awk '$9 >= %d || $9 <= -%d'") % (infile_name,\
               pair_reads_distance,pair_reads_distance)
        p = subprocess.Popen([cmd],shell=True,stdout=subprocess.PIPE,\
                             stderr=subprocess.PIPE)
        outfile,err = p.communicate()
    return outfile

def read_file(infile,insert_size,read_length):
    """put all supportive reads in a group"""
    fileline = infile.strip().split('\n')
    p_list = [[0,0,0,0]] 
    i = 0
    while i < len(fileline):
        line = fileline[i].strip().split('\t')[:9]
        result = check_line(p_list,line,insert_size,read_length)
        if result[0]:
            p_list += [result[1:]]
            i += 1
        else:
            if not i+1 == len(fileline):
                second_line = fileline[i+1].strip().split('\t')[:9]
                #check the second followed read
                result = check_line(p_list,second_line,insert_size,\
                                    read_length)
                if result[0]:
                    p_list += [result[1:]]
                    i += 2
                else:
                    yield p_list
                    p_list = [result[1:]]
                    i += 1
            else:#only for the last read
                yield p_list
                p_list = [result[1:]]
                i += 1
    return

def check_line(p_list,line,insert_size,read_length):
    """check if the reads belong to the group"""
    line_list = line[:]
    chr = line_list[2]
    loci1 = int(line_list[3])#start of alignment
    loci2 = int(line_list[7])#end of alignment
    distan = int(line_list[-1])#alignment distance of paired reads
    if chr == p_list[0][0] and abs(loci1-p_list[0][1]) <= insert_size\
            and abs(loci2-p_list[0][2]) <= insert_size and abs(distan-\
            p_list[0][3]) <= read_length:
        return True,chr,loci1,loci2,distan
    else:
        return False,chr,loci1,loci2,distan


def info(raw_list,insert_size,read_length):
    """"generate the start and end location for each group. Like chr start
    end num_sr length: SL2.40ch00 1200 1300 10 100.
    The list will be sorted"""
    result_list = []
    for i in raw_list:
        #filter out sr<=3 or length of deletion >= 20000
        if len(i) >= 3 and abs(i[0][-1]) <= 20000 : 
            if i[0][-1] < 0:
                result_list += [[i[0][0],i[-1][2],i[0][1],len(i),i[0][-1]]]
            else:
                result_list += [[i[0][0],i[-1][1],i[0][2],len(i),i[0][-1]]]
    sorted_list = sorted(result_list)
    result_list = [[0,0,0,0]]
    i = 0
    while i < len(sorted_list):
        chr = sorted_list[i][0]
        start = sorted_list[i][1]
        end = sorted_list[i][2]
        length = sorted_list[i][-1]
        if chr == sorted_list[i+1][0] and abs(start-sorted_list[i+1][1]) <=\
                insert_size and abs(end-sorted_list[i+1][2]) <= insert_size\
                and abs(length+sorted_list[i+1][-1]) <= read_length:
            #length of rightsr should be negative
            if length < 0:
                leftsr = sorted_list[i+1][3]
                rightsr = sorted_list[i][3]
            else:
                leftsr = sorted_list[i][3]
                rightsr = sorted_list[i+1][3]
            yield [chr,start,end,end-start,leftsr,rightsr]
            i +=2#skip the paired deletion
        else:
            i +=1

if  __name__=="__main__":
    bamfile = argv[1]
    outfile_name = argv[2]
    insert_size = int(argv[3])
    read_length = int(argv[4])
    pair_reads_distance = insert_size-1.5*read_length     
    reads_file = extract_reads(bamfile,pair_reads_distance)
    raw_list = read_file(reads_file,pair_reads_distance,read_length)
    outfile = open(outfile_name,'w+')
    for i in info(raw_list,pair_reads_distance,read_length):
        #the columns of output should be Chromosome Start End Length
        #Left_SR Right_SR
        outfile.write("%s\t%d\t%d\t%d\t%d\t%d\n" % tuple(i))
    outfile.close()

