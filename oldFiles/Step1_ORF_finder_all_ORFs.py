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

noduplicate = [] # initialize array/list

orf_dict = {} # initialize dictionary

for record in record_dict: # iterates through the sequences
	#print record_dict[record]
	orf_dict[record_dict[record].id] = [] # sets the record key/id as the key in the orf dictionary, the ">" line
	if len(record_dict[record].seq) > ((min_pro_len)*3): # if the orf length is > 3 times the specified protein length 
		for strand, nuc in [(+1, record_dict[record].seq), (-1,record_dict[record].seq.reverse_complement())]: # get a list of possible ORF translations, loops twice, each tuple
			for frame in range(3): #loops 3 times, frame will start at 0, 1, and 2
				trans = nuc[frame:].translate() # translates the sequence nuc[frame] into protein sequence for the current frame
				# find a way for it to stop at first stop codon
				new_nuc = nuc[frame:] # gets nucleotide sequence
				trans_len = len(trans) # gets length of the protein sequence
				seq_len = len(record_dict[record].seq) # gets the length of the nucleotide sequence from the fasta file 
				aa_start = 0
				while aa_start < trans_len: # while the starting amino acid is less than the peptide sequence length
					trans_start = trans.find("M", aa_start) # finds index of start aa sequence
					trans_end = trans.find("*", trans_start) # finds index of end aa sequence
					orf_len = trans_end-trans_start # length of current ORF being analyzed
					#print aa_start
					#print "strart: "+str(trans.find("*", aa_start))
					#print trans_start
					#print "-1: "+str(trans.find("*", trans_start))
					if trans_start == -1: # if no start is found
						trans_end = trans.find("*", aa_start) # find the end anyway between aa_start and end of seq
						if trans_end == -1: # if no end is found
							aa_start = trans_len # sets the start to the last index
						else: #if an end is found
							if trans_end - aa_start > min_pro_len: #if length from aa_start to the stop codon found is > min length
								#print trans_end
								#print trans_start
								# will treat the entire thing as an open reading frame, since there is a stop codon but the start codon isn't to be found... yet in the current frame
								if strand==1: # adds to dictionary for ORFs of 1 strands
									#print'1!!!'
									start = aa_start
									end = min(seq_len, frame+trans_end*3+3) # look into this!! whyyyy
									sequence = record_dict[record].seq[start:end]
									orf_dict[record_dict[record].id].append((sequence,start,end,True)) # check datastructure append
									ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_'+str(start)+'_'+str(end))
								if strand==-1: # adds to dictionary for ORFs of -1 strands
									start = max(0, seq_len-frame-(trans_end*3)-3)
									end = seq_len
									sequence = record_dict[record].seq[start:end].reverse_complement()
									ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
									orf_dict[record_dict[record].id].append((sequence,start,end,False))
							aa_start = trans_len
					elif trans_end == -1: #if a start is found but no end is found
						if trans_len-trans_start > min_pro_len: #if ORF from found start to end sequence is long enough
							#print trans_len-trans_start
							# same as above but tweaked for this case
							if strand ==1:
								#print'3!!!'
								start = frame+trans_start*3
								end = seq_len
								sequence = record_dict[record].seq[start:end]
								ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_'+str(start)+'_'+str(end))
								orf_dict[record_dict[record].id].append((sequence,start,end,True))
								
							if strand ==-1:
								#print'4!!!'
								start = 0
								end = seq_len-frame-trans_start*3
								sequence = record_dict[record].seq[start:end].reverse_complement()
								ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
								orf_dict[record_dict[record].id].append((sequence,start,end,False))
							aa_start = trans_len
						else:
							aa_start = trans_len
					elif orf_len > min_pro_len: # if we have both a start and an end and its long enough
						#print str(trans_start)+'\t'+str(trans_end)
						if strand ==1:
							#print'5!!!'
							start = frame+trans_start*3
							end = min(seq_len, frame+trans_end*3+3)
							end = frame+trans_end*3+3
							sequence = record_dict[record].seq[start:end]
							ORF_seq = SeqRecord(sequence, id=str(record_dict[record].id)+'_'+str(start)+'_'+str(end))
							orf_dict[record_dict[record].id].append((sequence,start,end,True))
							
						if strand ==-1:
							start = max(0, seq_len-frame-trans_end*3-3)
							#print start
							end = seq_len-frame-trans_start*3
							sequence = record_dict[record].seq[start:end].reverse_complement()
							#print record_dict[record].seq[start:end].reverse_complement().translate()
							ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
							orf_dict[record_dict[record].id].append((sequence,start,end,False))
						aa_start = trans_len
					else:
						aa_start = trans_len
#print orf_dict
for key in orf_dict.iterkeys(): #keys are contig names '>______'
	count=1 # keep track which ORF is in what contig bcause contigs may contain more than one
	for item in orf_dict[key]: #goes through the items in each key
		if item[3] == True:
			ORF_seq = SeqRecord(item[0], id=key+'_'+str(count)+'_'+str(item[1])+'_'+str(item[2]))
			# key = contig name, + count, + star position within contig, + end position within the contig
			# 
		if item[3] == False:
			ORF_seq = SeqRecord(item[0], id=key+'_'+str(count)+'_'+str(item[1])+'_'+str(item[2])+'_RC') # RC = reverse complement indication
		#print str(ORF_seq.seq)
		ORF_outfile.write('>'+str(ORF_seq.id)+'\n'+str(ORF_seq.seq)+'\n')
		#SeqIO.write(ORF_seq, ORF_outfile, "fasta")
		count = count + 1

ORF_outfile.close()

fixSeqRecord(argv[4])

noduplicate = []
for record in record_dict:
	if not record_dict[record].id in noduplicate:
		if record_dict[record].id in orf_dict.iterkeys():
			noduplicate.append(record_dict[record].id)
			SeqIO.write(record_dict[record], contig_outfile, "fasta")

trans_out_file = open(argv[4].split('.')[0]+'_trans.fa', 'w')

record_dict2 = SeqIO.index(argv[4],"fasta")

for record in record_dict2:
	trans_seq = record_dict2[record].seq.translate()
	trans_rec = SeqRecord(trans_seq, id= record_dict2[record].id)
	SeqIO.write(trans_rec, trans_out_file, "fasta")

trans_out_file.close()

fixSeqRecord(argv[4].split('.')[0]+'_trans.fa')

print len(orf_dict.keys())

	