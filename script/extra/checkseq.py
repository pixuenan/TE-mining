"""
Xuenan Pi
17-8-2015
Script to check if the sequence is in the fasta file
"""
from sys import argv
import string

def read_fasta(infile_name,inseq):
    seq = {}
    infile = open(infile_name)
    for line in infile:
        if line.startswith('>'):
            header = line.strip()
            seq[header] = ''
        else:
            seq[header] += line
    for key,value in seq.items():
            #if sequen[:5]==sequen[-5:]:
            #    print 'TSD=Y'
            #else:
            #    print sequen[:5],sequen[-5:]
            #  
            #  [:130][-330:-270]
        if inseq in value or inseq in complim(value):
            print key
    return

def complim(seq):
    t = string.maketrans('ATCG','TAGC')
    trans = seq.translate(t)
    compli = list(trans)
    compli.reverse()
    return ''.join(compli)


if __name__ == "__main__":
    inseq1 = 'ggggggtccaggggggcgacgcgccccctggcctgggggtccgggggg'
    inseq2 = 'GGGGGGTCCAGGGGGGCGACGCGCCCCCTGGCCTGGGGGTCCGGGGGG'
    read_fasta(argv[1],inseq2)
