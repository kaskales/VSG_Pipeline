#!/usr/bin/env python
import argparse
import subprocess
from Bio import SeqIO 
from sys import argv
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

parser.add_argument('-r', metavar='reference assembly - FASTA', action="store", dest="r")
parser.add_argument('-n', metavar='assembly number for MULTo, ie, 1293 for \'-a tb1293\' in MULTo run', action ="store", dest = "n") 
parser.add_argument('-p', help='path to MULTo1.0 folder. default is ~/bioinformatics_resources/', action="store", dest='p', default='~/bioinformatics_resources/')

arguments = parser.parse_args()

ref= arguments.r
number=arguments.n
path=arguments.p

def make_multo_dirs(r,n,p):
    samplename = r.split('.')[0]
    tbnumber = str(n)
    path = str(p)

    #make multo file heirarchy from given basename
    subprocess.call(['mkdir -p '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/'], shell=True)

    subprocess.call(['mkdir -p '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)


def make_bed(r,n,p):
    record_dict = SeqIO.index(r,"fasta")

    concatFA = open('chr1.fa', 'w')
    BEDfile = open('chr1.bed', 'w')
    
    name = 'chr1'
    #name of chr

    start = 0
    end = 0

    N = 'N'*100

    concatFA.write('>'+name+'\n')

    for record in record_dict:
            concatFA.write(str(record_dict[record].seq)+str(N))
            seqlen = len(record_dict[record].seq)
            end = start + seqlen
            seqname = record_dict[record].id
            BEDfile.write(name+'\t'+str(start)+'\t'+str(end)+'\t'+str(seqname)+'\t 0 \t + \t'+str(start)+'\t'+str(end)+'\n')
            start = end + 100
            
    concatFA.close()
    BEDfile.close()
    
def move_multo(r,n,p):
    samplename = r.split('.')[0]
    tbnumber = str(n)
    path = str(p)
    #move concat and bedfile into new folders
    subprocess.call(['mv chr1.fa '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)

    subprocess.call(['mv chr1.bed '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/'], shell=True)

    #run MULTo 

    subprocess.call(['python '+path+'MULTo1.0/src/MULTo1.0.py -s tbb -a tb'+tbnumber+' -v 2'], shell=True)

make_multo_dirs(ref,number,path)
make_bed(ref,number,path)
move_multo(ref,number,path)

