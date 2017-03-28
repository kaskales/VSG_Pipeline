from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
from sys import argv
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

def fixSeqRecord(file): # gets rid of not usefull stuff in the file comment ">" section 
	infile = open(file, 'r')
	outfile = open(str(file.split('.')[0])+'_clean.fa', 'w')
	for line in infile:
		if '>' in line:
			names = line.split()
			for name in names:
				if name[0] == '>':
					outfile.write(name+'\n')	
		else: 
			outfile.write(line)
		
	outfile.close()
	infile.close()

def addSeqRecord(recordRD, start, end, count):
	sequence = recordRD.seq[start:end]
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'\n'+str(sequence)+'_'+filename+'\n')
	SeqIO.write(SeqRecord(sequence.translate(), id=str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_'+filename), trans_out_file, "fasta")
def addSeqRecord_RC(recordRD, start, end, count):
	sequence = recordRD.seq[start:end].reverse_complement()
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC'+'\n'+str(sequence)+'_'+filename+'\n')
	SeqIO.write(SeqRecord(sequence.translate(), id=str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC'+'_'+filename), trans_out_file, "fasta")


parser = argparse.ArgumentParser()

parser.add_argument('-t', nargs= '+', metavar='Trinity files to be put through the pipeline', action="store", dest="t") # list of trinity files to be put through the pipline
parser.add_argument('-d', help='additional descriptive terms to name your run', action="store", dest='d', default='')
parser.add_argument('-m', metavar='minimum protein length you are filtering for', action ="store", dest = "m") 
parser.add_argument('-p', help='path to MULTo1.0 folder. default is ~/', action="store", dest='p', default='~/') # default assumes MULTo is in your home dir

arguments = parser.parse_args()

timeRan = time.strftime("%d/%m/%Y-%H:%M_") + str(arguments.d) # Day/Month/Year-Hour:Minute , names the folder for output files
os.makedirs(timeRan) # creates the folder
contig_outfile = open(os.path.join(timeRan, timeRan+"contig.fa"), 'w') # contig output file
global ORF_outfile
ORF_outfile = open(os.path.join(timeRan, timeRan+"orf.fa"), 'w') # orf output file 
global trans_out_file
trans_out_file = open(os.path.join(timeRan, timeRan+'orf_trans.fa'), 'w') # translated orf file
min_pro_len = int(arguments.m)


for file in arguments.t: # loops all the trinity files to find orf, will add everything to the same file
	record_dict = SeqIO.index(file,"fasta") # "parses" a fasta file, creating a dictionary-like object of sequences. not everything is kept in memeory. instead it just records where each record is within the file. parses on demand. 
	# the key is the dictionary is the ">" line in the fasta file
	 # minimum protein length, within typical VSG protein length, in a.a.
	global filename
	filename = str(file).split('.')[0] # takes only the name of the file, this will be added to the ORF description, so we know where it came from
	for record in record_dict: # iterates through the sequences
		#print record_dict[record]
		if len(record_dict[record].seq) > ((min_pro_len)*3): # if the orf length is > 3 times the specified protein length 
			count=1 
			for strand, nuc in [(+1, record_dict[record].seq), (-1,record_dict[record].seq.reverse_complement())]: # get a list of possible ORF translations, loops twice, each tuple
				for frame in [0,1,2]: #loops 3 times, frame will start at 0, 1, and 2

					trans = nuc[frame:].translate() # translates the sequence nuc[frame] into protein sequence for the current frame
					# find a way for it to stop at first stop codon
					trans_len = len(trans) # gets length of the protein sequence
					seq_len = len(record_dict[record].seq) # gets the length of the nucleotide sequence from the fasta file 						
					trans_end = trans.find("*", 0)
					trans_start = trans.find("M", 0)
					if trans_start != -1: #start codon exists
						if trans_end != -1: # end codon exists
							if trans_start < trans_end: # there is a start codon before the 1st stop codon
								if trans_end - trans_start > min_pro_len: # check if this protein is long enought to be considered
									#print "A"
									if strand ==1:
										addSeqRecord(record_dict[record], frame+(trans_start*3), frame+(trans_end*3)+3, count)
									else:
										addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count)
									count += 1
								elif trans_end > min_pro_len : # protein isn't long enough, so just take the entire front end if thats long enough
									#print "B"
									if strand ==1:
										addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count)
									else:
										addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len, count)
									count += 1
								# if neither protein is long enough,or even if it was, we will move on to find where the next start codon is, and find the end codon for that one
								trans_start = trans.find("M", trans_end)
								trans_end = trans.find("*", trans_start)		
							else: # there is a start codon AFTER the first appearing stop codon
								if trans_end > min_pro_len: # is the seq_start to end codon is long enough for the a min length protein
									#print "C"
									if strand ==1:
										addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count)
									else:
										addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len, count)
									count += 1
								trans_end = trans.find("*", trans_start) # will start checking the orf at that first found start codon, sets new trans_end to stop codon after the start	
						else: # start codon exists, but there is no stop codon anywhere! :O
						# this portion of the code will take the entire thing as an ORF, we are being generous
							addSeqRecord_RC(record_dict[record], 0, seq_len, count)
							count += 1
						# this part isn't generous and only takes from the start codon till the end of the sequence
	#						if trans_len-trans_start > min_pro_len: # is the remaining protein segment after the start codon long enought to be a protein?
	#							print "D"
	#							print trans_start
	#							if strand ==1:
	#								addSeqRecord(record_dict[record], frame+(trans_start*3), seq_len, count)
	#							else:
	#								addSeqRecord_RC(record_dict[record], frame, seq_len-frame-(trans_start*3), count)
	#							count += 1
							trans_start = trans_len # won't go into loop, bc the entire rest of the sequence just got taken as an orf, so there is no point	in searching further, no more to search
					elif trans_end != -1: # start codon doesn't exist, but end codon does
						if trans_end > min_pro_len: # take entire sequence from start of seq to the found end codon as ORF if longer than min length
							#print "E"
							if strand ==1:
								addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count)
							else:
								addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count)
							count += 1	
						trans_start = trans_len # won't search for more, since there aren't any start codons anywhere
					else: # there are no start or stop codons anywhere?!
						addSeqRecord_RC(record_dict[record], 0, seq_len, count)
						trans_start = trans_len
					
					trans_max = trans_len - min_pro_len # farthest a codon can be before it's below min protein length

					while trans_start < trans_max:
						if trans_start == -1: # if no more starts are found, doesn't matter if there is a stop or not, game over
							trans_start = trans_len
						elif trans_end == -1: #if start is found but no end is found
							if trans_len-trans_start > min_pro_len: #if ORF from found start to end sequence is long enough
								#print "G"
								if strand ==1:
									addSeqRecord(record_dict[record], frame+(trans_start*3), seq_len, count)
								else:
									addSeqRecord_RC(record_dict[record], frame, seq_len-frame-(trans_start*3), count)
								count += 1
							trans_start = trans_len
						elif trans_end-trans_start > min_pro_len: # if we have both a start and an end and its long enough
							#print "H"
							if strand ==1:
								addSeqRecord(record_dict[record], frame+trans_start*3, frame+(trans_end*3)+3, count)
							else:
								addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count)
							count += 1	
							trans_start = trans.find("M", trans_end) # finds index of start aa sequence
							trans_end = trans.find("*", trans_start)	
						else:
							trans_start = trans.find("M", trans_end) 
							trans_end = trans.find("*", trans_start)			

			if count > 1: # have any new orf been added? if so, add this record to file
				SeqIO.write(record_dict[record], contig_outfile, "fasta") 

ORF_outfile.close()
fixSeqRecord(argv[4])
trans_out_file.close()
fixSeqRecord(argv[4].split('.')[0]+'_trans.fa')
