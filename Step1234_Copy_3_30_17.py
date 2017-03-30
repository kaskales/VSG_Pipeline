from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
from sys import argv
import vsganalysisfunctions as vsg
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
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'\n'+str(sequence)+'\n')
	SeqIO.write(SeqRecord(sequence.translate(), id=str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)), trans_out_file, "fasta")
def addSeqRecord_RC(recordRD, start, end, count):
	sequence = recordRD.seq[start:end].reverse_complement()
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC'+'\n'+str(sequence)+'\n')
	SeqIO.write(SeqRecord(sequence.translate(), id=str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC'), trans_out_file, "fasta")


parser = argparse.ArgumentParser()

parser.add_argument('-t', nargs= '+', metavar='Trinity files to be put through the pipeline', action="store", dest="t") # list of trinity files to be put through the pipline
parser.add_argument('-d', help='additional descriptive terms to name your run', action="store", dest='d', default='')
parser.add_argument('-m', metavar='minimum protein length you are filtering for', action ="store", dest = "m", default=300) 
parser.add_argument('-p', help='path to MULTo1.0 folder. default is ~/', action="store", dest='p', default='~/') # default assumes MULTo is in your home dir
parser.add_argument('-cpu', help='number of processors', action="store", dest='cpu', default='2')
parser.add_argument('-v', help='number of mismatches allowed', action="store", dest='v', default='2')
parser.add_argument('-noqc', help='QC will not be run, default is on', action ='store_false', dest = 'qc' , default = False)
group = parser.add_mutually_exclusive_group()
group.add_argument('-m_only', help='only run rpkmforgenes/multo scripts, input is samfile', action ='store_false', dest = 'm_only', default = False)
group.add_argument('-qc_only', help= 'only run QC, default is off', action= 'store_true', dest='qc_only', default = False)
parser.add_argument('-i', help='list of files - FASTQ, or .sam is bowtie and QC is turned off', action="store", dest="i", required=True)


arguments = parser.parse_args()
global timeRan
timeRan = time.strftime("%d-%m-%Y-%H_%M")  # Day/Month/Year-Hour:Minute , names the folder for output files
if arguments.d != '':
	timeRan = timeRan + "-"+ str(arguments.d)

if not os.path.exists(timeRan):
	os.makedirs(timeRan) # creates the folder

contig_outfile = open(os.path.join(timeRan, timeRan+"_contig.fa"), 'w') # contig output file
global ORF_outfile
ORF_outfile = open(os.path.join(timeRan, timeRan+"_orf.fa"), 'w') # orf output file 
global trans_out_file
trans_out_file = open(os.path.join(timeRan, timeRan+'_orf_trans.fa'), 'w') # translated orf file
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
fixSeqRecord(timeRan + "/"+timeRan+"_orf.fa")
trans_out_file.close()
fixSeqRecord(timeRan + "/"+timeRan+'_orf_trans.fa')

# Step 2 

def blast_sort(v,n,s):
    #v = vsg xml file
    #n = nonvsg xmlfile
    #s = sequence file, contigs from blast searches
    result_handle = open(v)
    nonVSGresult_handle = open(n)
    blast_records = NCBIXML.parse(result_handle) # returns an iterator of the blast results
    #blast_record = blast_records.next()
    record_dict = SeqIO.index(s,"fasta")
    
    outfile = open(s.split('.')[0]+'_VSGs.fa', 'w')
    scorefile = open(s.split('.')[0]+'_VSGs_scores.fa', 'w')

    blast_records_nonVSG = NCBIXML.parse(nonVSGresult_handle)
    blast_record_nonVSG = blast_records_nonVSG.next()
    
    hit_list = [] # list of VSGs we have found! 
    exclude_list = []
        
    print 'Now looking for non-VSG transcripts...'
    
    for blast_record_nonVSG in blast_records_nonVSG:
    	for alignment in blast_record_nonVSG.alignments:
    		for hsp in alignment.hsps: 
    			percent_identity = (100.0 * hsp.identities) / alignment.length # hsp.identities is a tuple(bp matches, total bp in seq) to give percent match of sequence, percent identity is # of bp
    			percent_query_identity = (100.0 * hsp.identities) / blast_record_nonVSG.query_letters
    			#print blast_record_nonVSG.query+'\t'+alignment.title+'\t'+str(percent_identity)+'\t'+str(percent_query_identity)+'\t'
    			if (percent_query_identity > 30 and hsp.identities > 300) or (percent_identity > 90):
    				if not blast_record_nonVSG.query in exclude_list:
    					exclude_list.append(str(blast_record_nonVSG.query))
    					#print 'nonVSG hit!'+'\t'+str(blast_record_nonVSG.query)+' \t '+str(alignment.title)

    print 'VSG hits! - maybe?'
    
    for blast_record in blast_records:
    	for alignment in blast_record.alignments:
    		for hsp in alignment.hsps:
    			if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
    				if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
    					if not blast_record.query in exclude_list: # if the query isn't a fake VSG hit, add it now!
	    					hit_list.append(str(blast_record.query))											# percent query aligned										# percent identity
	    					scorefile.write(str(blast_record.query)+'\t'+str(alignment.title)+'\t'+str((100.0 * hsp.identities) / blast_record.query_letters)+'\t'+str((100.0 * hsp.identities) / alignment.length)+'\t'+str(alignment.length)+'\n')
	    					SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
                    
    outfile.close
    scorefile.close


print ' *****analyzing '+timeRan + "/"+timeRan+"_orf.fa"+' *****'
filename = timeRan + "/"+timeRan+"_orf"
#blast VSG
subprocess.call(['blastn -db tb427_vsgs -query '+str(filename)+'.fa -outfmt 5 -out '+str(filename)+'.xml'], shell=True)
#blast nonVSG
subprocess.call(['blastn -db NOTvsgs -query '+str(filename)+'.fa -outfmt 5 -out '+str(filename)+'_nonVSG.xml'], shell=True)
#get all the blast results which are for ONLY VSGs, get rid of hits which are VSG-similar but not vsgs
blast_sort(str(filename)+'.xml', str(filename)+'_nonVSG.xml',filename+".fa")
# cdhit merge
subprocess.call(['cd-hit-est -i '+filename+'_VSGs.fa '+' -d 0 -o '+filename+'_merged.fa -c 0.9 -n 8 -r 1 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.5'], shell=True)


# Step 3

currDir = os.getcwd()
print currDir	
path = str(arguments.p)
os.path.expanduser('~')
os.chdir(os.path.expanduser('~'))
print os.getcwd()
os.chdir('MULTo1.0')

# make multo dir
tbnumber = timeRan

#make multo file heirarchy from given basename
subprocess.call(['mkdir -p '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/'], shell=True)
subprocess.call(['mkdir -p '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)

# make bed
record_dict = SeqIO.index(currDir +'/' + filename+'_VSGs.fa',"fasta")
concatFA = open('chr1.fa', 'w')
BEDfile = open('chr1.bed', 'w')  
name = 'chr1'
#name of chr
start = 0
end = 0
N = 'N'*100
concatFA.write('>'+name+'\n')
for record in record_dict:
    concatFA.write(str(record_dict[record].seq)+str(N)) # adds lots of N's to the end of the sequence
    seqlen = len(record_dict[record].seq)
    end = start + seqlen
    seqname = record_dict[record].id
    BEDfile.write(name+'\t'+str(start)+'\t'+str(end)+'\t'+str(seqname)+'\t 0 \t + \t'+str(start)+'\t'+str(end)+'\n')
    start = end + 100          
concatFA.close()
BEDfile.close()
   
# move multo 
#move concat and bedfile into new folders
subprocess.call(['mv chr1.fa '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)
subprocess.call(['mv chr1.bed '+path+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/'], shell=True)
#run MULTo 
subprocess.call(['python '+path+'MULTo1.0/src/MULTo1.0.py -s tbb -a tb'+tbnumber+' -v '+arguments.v+' -O -p '+arguments.cpu], shell=True)

os.chdir(currDir) # sets working directory back to VSG____ folder
# Step 4


# file = FASTQ file

n = tbnumber
p = path
qc_only = arguments.qc_only
m_only = arguments.m_only
file = timeRan + "/" + arguments.i

if qc_only==True: # default is false, don't worry about it for now
	#vsg.QC(file)  
	print ""		
elif m_only==True: # defaults is false
	print file
	vsg.count(file,n,p)
else:
	if arguments.qc==True: # default is false
		vsg.QC(file)
		file = file.split('.')[0]+'_trimmed3.fq'
		vsg.bowtie(file,n,p)
		file =file.split('.')[0]+'_align.sam'
	else:
		vsg.bowtie(file,n,p)
		file =file.split('.')[0]+'_align.sam'
	print file
	vsg.count(file,n,p)
    














