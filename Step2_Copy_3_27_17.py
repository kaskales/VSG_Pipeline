import argparse
import subprocess
import vsgfindfunctions as vf
from Bio import SeqIO 
from sys import argv
from sys import argv
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import subprocess

def blast_sort(v,n,s):
    #v = vsg xml file
    #n = nonvsg xmlfile
    #s = sequence file, contigs from blast searches
    result_handle = open(v)
    nonVSGresult_handle = open(n)
    blast_records = NCBIXML.parse(result_handle)
    #blast_record = blast_records.next()
    record_dict = SeqIO.index(s,"fasta")
    
    outfile = open(s.split('.')[0]+'_VSGs.fa', 'w')
    
    blast_records_nonVSG = NCBIXML.parse(nonVSGresult_handle)
    blast_record_nonVSG = blast_records_nonVSG.next()
    
    hit_list = [] # list of VSGs we have found! 
    exclude_list = []
        
    print 'Now looking for non-VSG transcripts...'
    
    for blast_record_nonVSG in blast_records_nonVSG:
    	for alignment in blast_record_nonVSG.alignments:
    		for hsp in alignment.hsps: 
    			percent_identity = (100.0 * hsp.identities) / alignment.length
    			percent_query_identity = (100.0 * hsp.identities) / blast_record_nonVSG.query_letters
    			print blast_record_nonVSG.query+'\t'+alignment.title+'\t'+str(percent_identity)+'\t'+str(percent_query_identity)+'\t'
    			if percent_query_identity > 30 and hsp.identities > 300:
    				if not blast_record_nonVSG.query in exclude_list:
    					exclude_list.append(str(blast_record_nonVSG.query))
    					#print 'nonVSG hit!'+'\t'+str(blast_record_nonVSG.query)+' \t '+str(alignment.title)
    			if percent_identity > 90:
    				if not blast_record_nonVSG.query in exclude_list:
    					exclude_list.append(str(blast_record_nonVSG.query))
                        #print 'nonVSG hit!'+'\t'+str(blast_record_nonVSG.query)+' \t '+str(alignment.title)

    print 'VSG hits! - maybe?'
    
    for blast_record in blast_records:
    	for alignment in blast_record.alignments:
    		for hsp in alignment.hsps:
    			if hsp.expect < 1.0e-10: 
    				if not blast_record.query in hit_list:
    					if not blast_record.query in exclude_list:
	    					hit_list.append(str(blast_record.query))
	    					percent_query_aligned = (100.0 * hsp.identities) / blast_record.query_letters
	    					percent_identity = (100.0 * hsp.identities) / alignment.length
	    					print record_dict[blast_record.query]
	    					SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
                    
    outfile.close


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
#	vf.blast_vsg(file,l)
#	vf.blast_nonvsg(file,l)
#	Vfile = file.split('.')[0]+'.xml' # where are these files created?
#	Nfile = file.split('.')[0]+'_nonVSG.xml'
#	vf.blast_sort(Vfile,Nfile,file)
#	newfile = file.split('.')[0]+'_VSGs.fa'
#	vf.cdhitmerge(newfile)
	filename = str(file).split('.')[0]
	#blast VSG
	subprocess.call(['blastn -db tb427_vsgs -query '+str(filename)+'.fa -outfmt 5 -out '+str(filename)+'.xml'], shell=True)
	#blast nonVSG
	subprocess.call(['blastn -db NOTvsgs -query '+str(filename)+'.fa -outfmt 5 -out '+str(filename)+'_nonVSG.xml'], shell=True)
	#
	blast_sort(file.split('.')[0]+'.xml', file.split('.')[0]+'_nonVSG.xml',file)
	# cdhit merge
	subprocess.call(['cd-hit-est -i '+filename+'_VSGs.fa '+' -d 0 -o '+filename+'_merged.fa -c 0.9 -n 8 -r 1 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.5'], shell=True)










