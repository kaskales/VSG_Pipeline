import argparse
import subprocess
import vsgfindfunctions as vf
from Bio import SeqIO 
from sys import argv



parser = argparse.ArgumentParser()
# parses user input

parser.add_argument('-i', nargs= '+', metavar='list of files - FASTA', action="store", dest="i")
# nargs = number of arguments to accept, '+' accepts all, but has a requirement of 1 argument

#parser.add_argument('-', metavar='reference assembly - FASTA',  action="store", dest="r")
parser.add_argument('-l', metavar='minimum ORF length for analysis, default is 900',type=int, action ="store", dest = "l", default=900) 
arguments = parser.parse_args()

l=arguments.l

for file in arguments.i:
	print ' *****analyzing '+str(file)+' *****'
	vf.blast_vsg(file,l)
	vf.blast_nonvsg(file,l)
	Vfile = file.split('.')[0]+'.xml' # where are these files created?
	Nfile = file.split('.')[0]+'_nonVSG.xml'
	vf.blast_sort(Vfile,Nfile,file)
	newfile = file.split('.')[0]+'_VSGs.fa'
	vf.cdhitmerge(newfile)