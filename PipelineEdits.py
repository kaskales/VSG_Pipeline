from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
from sys import argv
import VSGFunctions as vsgf
import subprocess
import argparse
import time
import os

# usage is ORF_finder.py infile minimum protein length outfile
# this goes through a FASTA file and finds open reading frames - as soon as 
# an ORF is found, it writes that record to a new FASTA file and stops looking in
# that record.
#if no start is found, beginning -> stop is an ORF
#if no stop is found, start -> end is an ORF

# finds all open reading frames within a specified length

parser = argparse.ArgumentParser()

# arguments about input files and naming
parser.add_argument('-s', nargs= '+',metavar='listed sequence files to be put through the pipeline, fq', action="store", dest="s")
parser.add_argument('-st',metavar='text file with names of sequence files to be put through the pipeline, fq', action="store", dest="st", default='')
parser.add_argument('-d', help='additional descriptive terms to name your run', action="store", dest='d', default='')
# trimming setting
parser.add_argument('-g', metavar='stringency for trim galore', action ="store", dest = "g", default="3") 
# trinity settings
parser.add_argument('-minp', metavar='minimum protein length you are filtering for', action ="store", type=int, dest = "minp", default=300) 
parser.add_argument('-mem', metavar='max memory allocation for trinity', action ="store", dest = "mem", default="10") 
parser.add_argument('-cpu', help='number of processors', action="store", dest='cpu', default='2')
# Blast settings
parser.add_argument('-vsgdb', metavar='name of the vsg database', action ="store", dest = "vsgdb", default="tb427_vsgs")
# cd-hit-est parameters
parser.add_argument('-sit', metavar='sequence identiy threshold - how much the alignment has to match. value is 0.0 through 1.0 ', action ="store", dest = "sit", default=".98")
# MULTO settings
parser.add_argument('-p', help='path to MULTo1.0 folder. default is /Users/mugnierlab/, please dont use "~/", python doesnt like this in the path', action="store", dest='p', default='/Users/mugnierlab/') # default assumes MULTo is in your home dir
parser.add_argument('-v', help='number of mismatches allowed', action="store", dest='v', default='2')
parser.add_argument('-remakeMulto', help='name of the multo files, typically same as header, if default it will make MULTo files. otherwise enter name of multo files to be reused', action="store", dest='rmulto', default='')

#where the pipeline will start
parser.add_argument('-start', help='the step you want the pipeline to start at. 0 = input is raw untrimmed data. 1 = start after trimming, input is already trimmed. 2 = Start after Trinity, input are trinity files. 3 = start after finding ORF, input are ORFs. 4= Start after BLAST/cd-hit-est, input are VSG/ORFs. 5 = continues till the end to MULTo', type=int, action="store", dest='start', default=0)
parser.add_argument('-header', help='input header variable in the format of Y-m-d-H_M , if it doesnot have ', action="store", dest='head', default='')

#where the pipeline will stop
parser.add_argument('-stop', help='the step you want the pipeline to stop at. 1 = stop after trimming. 2 = Stop after Trinity. 3= stop after finding ORF. 4= Stop after BLAST, 5 = continues till the end to MULTo', type=int, action="store", dest='stop', default=5)

# add segment to tkae sequence data, trim with trim_galore and cutadapt

# run trimmed sequences trough trinity

# then take the trimmed trinity sequences and run through this pipeline

arguments = parser.parse_args()


# start from the very beguinning 
if arguments.start == 0:
	header = time.strftime("%Y-%m-%d-%H_%M")  # Year/Month/Day-Hour:Minute , names the folder for output files
	if arguments.d != '':
		for d in arguments.d:
			header = header + "-"+ str(d)
	if not os.path.exists(header):
		os.makedirs(header) # creates the folder
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	vsgf.trimSequences(header, trinityfiles, arguments)
	if arguments.stop >1:
		vsgf.trinity(header, trinityfiles, arguments)
		if arguments.stop > 2:
			vsgf.findORFs(header, trinityfiles, arguments)
			if arguments.stop > 3:
				vsgf.blastCDHIT(header, trinityfiles, arguments)
				if arguments.stop > 4:
					vsgf.makeMulto(header, trinityfiles, arguments)
# start with trimmmed files
elif arguments.start == 1:
	header = arguments.head  # Year/Month/Day-Hour:Minute , names the folder for output files

	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop >1:
		vsgf.trinity(header, trinityfiles, arguments)
		if arguments.stop > 2:
			vsgf.findORFs(header, trinityfiles, arguments)
			if arguments.stop > 3:
				vsgf.blastCDHIT(header, trinityfiles, arguments)
				if arguments.stop > 4:
					vsgf.makeMulto(header, trinityfiles, arguments)
# start with trinity files
elif arguments.start == 2:
	header = arguments.head  # Year/Month/Day-Hour:Minute , names the folder for output files

	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop > 2:
		vsgf.findORFs(header, trinityfiles, arguments)
		if arguments.stop > 3:
			vsgf.blastCDHIT(header, trinityfiles, arguments)
			if arguments.stop > 4:
				vsgf.makeMulto(header, trinityfiles, arguments)
# start with ORFs
elif arguments.start == 3:
	header = arguments.head  # Year/Month/Day-Hour:Minute , names the folder for output files

	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop > 3:
		vsgf.blastCDHIT(header, trinityfiles, arguments)
		if arguments.stop > 4:
			vsgf.makeMulto(header, trinityfiles, arguments)
# start with BLAST-ed VSG ORFs
elif arguments.start == 4:
	header = arguments.head  # Year/Month/Day-Hour:Minute , names the folder for output files
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop > 4:
		vsgf.makeMulto(header, trinityfiles, arguments)
else:
	print("Bad start argument, double check the number you entered")













