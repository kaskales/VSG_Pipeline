from Bio import SeqIO 
from sys import argv
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
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

record_dict = SeqIO.index(argv[1],"fasta") # "parses" a fasta file, creating a dictionary-like object of sequences. not everything is kept in memeory. instead it just records where each record is within the file. parses on demand. 
# the key is the dictionary is the ">" line in the fasta file
min_pro_len = int(argv[2]) # minimum protein length, within typical VSG protein length, in a.a.

contig_outfile = open(argv[3], 'w') # initializing contig output file, reference

ORF_outfile = open(argv[4], 'w') # initializing open reading frame output file, what we like

for record in record_dict: # iterates through the sequences
	#print record_dict[record]
	if len(record_dict[record].seq) > ((min_pro_len)*3): # if the orf length is > 3 times the specified protein length 
		count=1 
		for strand, nuc in [(+1, record_dict[record].seq), (-1,record_dict[record].seq.reverse_complement())]: # get a list of possible ORF translations, loops twice, each tuple
			for frame in range(3): #loops 3 times, frame will start at 0, 1, and 2
				trans = nuc[frame:].translate() # translates the sequence nuc[frame] into protein sequence for the current frame
				# find a way for it to stop at first stop codon

				trans_len = len(trans) # gets length of the protein sequence
				seq_len = len(record_dict[record].seq) # gets the length of the nucleotide sequence from the fasta file 
				
				#while aa_start < trans_len: # while the starting amino acid is less than the peptide sequence length
				trans_start = trans.find("M", 0) # finds index of start aa sequence
				trans_end = trans.find("*", trans_start) # finds index of end aa sequence which comes after the start aa seq

				if trans_start == -1: # if no start is found
					trans_end = trans.find("*", 0) # find the end anyway between aa_start and end of seq
					if trans_end != -1: # if end is found
						if trans_end > min_pro_len: #if length from aa_start to the stop codon found is > min length
							# will treat the entire thing as an open reading frame, since there is a stop codon but the start codon isn't to be found... yet in the current frame
							if strand==1: # adds to dictionary for ORFs of 1 strands
								end = min(seq_len, frame+trans_end*3+3) 
								sequence = record_dict[record].seq[0:end]
								ORF_outfile.write('>'+str(record_dict[record].id)+'_'+str(count)+'_'+str(0)+'_'+str(end)+'\n'+str(sequence)+'\n')
							else: # adds to dictionary for ORFs of -1 strands
								start = max(0, seq_len-frame-(trans_end*3)-3)
								sequence = record_dict[record].seq[start:seq_len].reverse_complement()
								ORF_outfile.write('>'+str(record_dict[record].id)+'_'+str(count)+'_'+str(start)+'_'+str(seq_len)+'_RC'+'\n'+str(sequence)+'\n')
							count += 1
				elif trans_end == -1: #if a start is found but no end is found
					if trans_len-trans_start > min_pro_len: #if ORF from found start to end sequence is long enough
						if strand ==1:
							start = frame+trans_start*3
							sequence = record_dict[record].seq[start:seq_len]
							ORF_outfile.write('>'+str(record_dict[record].id)+'_'+str(count)+'_'+str(start)+'_'+str(seq_len)+'\n'+str(sequence)+'\n')
						else:
							end = seq_len-frame-trans_start*3
							sequence = record_dict[record].seq[0:end].reverse_complement()
							ORF_outfile.write('>'+str(record_dict[record].id)+'_'+str(count)+'_'+str(0)+'_'+str(end)+'_RC'+'\n'+str(sequence)+'\n')
						count += 1

				elif trans_end-trans_start > min_pro_len: # if we have both a start and an end and its long enough
					if strand ==1:
						start = frame+trans_start*3
						end = frame+trans_end*3+3
						sequence = record_dict[record].seq[start:end]
						ORF_outfile.write('>'+str(record_dict[record].id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'\n'+str(sequence)+'\n')
					else:
						start = max(0, seq_len-frame-trans_end*3-3)
						end = seq_len-frame-trans_start*3
						sequence = record_dict[record].seq[start:end].reverse_complement()
						ORF_outfile.write('>'+str(record_dict[record].id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC'+'\n'+str(sequence)+'\n')
					count += 1

		if count > 1: # have any new orf been added? if so, add this record to file
			SeqIO.write(record_dict[record], contig_outfile, "fasta") 

ORF_outfile.close()

fixSeqRecord(argv[4])

trans_out_file = open(argv[4].split('.')[0]+'_trans.fa', 'w')

record_dict2 = SeqIO.index(argv[4],"fasta")

for record in record_dict2:
	trans_seq = record_dict2[record].seq.translate()
	trans_rec = SeqRecord(trans_seq, id= record_dict2[record].id)
	SeqIO.write(trans_rec, trans_out_file, "fasta")

trans_out_file.close()

fixSeqRecord(argv[4].split('.')[0]+'_trans.fa')