#!/usr/bin/env python
from Bio import SeqIO 
from sys import argv
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import subprocess


def fix(n):
    infile = open(n, 'r')
    outfile = open(n.split('.')[0]+'_clean.fa', 'w')

    for line in infile:
        if '>' in line:
            names = line.split()
            for name in names:
                if '>' in name:
                    outfile.write(name+'\n')	
        else: 
            outfile.write(line)
                
    outfile.close()
    infile.close()
        

def sizesort(n,l):
    record_dict = SeqIO.index(n,"fasta")	
    outfile = open(n.split('.')[0]+'_'+str(l)+'.fa', 'w')

    for record in record_dict:
        if len(record_dict[record].seq) > l:
            SeqIO.write(record_dict[record], outfile, "fasta")
            
def ORF_find(s,l):
    record_dict = SeqIO.index(s,"fasta")
    min_seq_len = l
    contig_outfile = open(str(s.split('.')[0])+'orf_contigs.fa', 'w')
    ORF_outfile = open(str(s.split('.')[0])+'orfs.fa', 'w')
    noduplicate = []

    for record in record_dict:
        if len(record_dict[record].seq) > (l):
            for strand, nuc in [(+1, record_dict[record].seq), (-1,record_dict[record].seq.reverse_complement())]:
                for frame in range(3):
                    trans = str(nuc[frame:].translate())
                    trans_len = len(trans)
                    seq_len = len(record_dict[record].seq)
                    #print seq_len
                    aa_start = 0
                    while aa_start < trans_len:
                        trans_start = trans.find("M", aa_start)
                        trans_end = trans.find("*", trans_start)
                        #print trans_start
                        #print trans_end
                        #print trans[trans_start:trans_end+1]
                        orf_len = trans_end-trans_start
                        if trans_start == -1:
                            aa_start = trans_len
                        elif trans_end == -1:
                            if trans_len-trans_start > l/3 :
                                #print record_dict[record].id +'\t'+str(trans_len-trans_start) 
                                if not record_dict[record].id in noduplicate:
                                    noduplicate.append(record_dict[record].id)
                                    SeqIO.write(record_dict[record], contig_outfile, "fasta")
                                if strand ==1:
                                    start = frame+trans_start*3
                                    end = seq_len
                                    sequence = record_dict[record].seq[start:end]
                                    ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_'+str(start)+'_'+str(end))
                                    SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                                if strand ==-1:
                                    start = 0
                                    end = seq_len-frame-trans_start*3
                                    sequence = record_dict[record].seq[start:end].reverse_complement()
                                    ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
                                    SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                                aa_start = trans_len
                            else:
                                aa_start = trans_len
                        elif orf_len > l/3 :
                            #print record_dict[record].id +'\t'+str(orf_len)
                            if not record_dict[record].id in noduplicate:
                                noduplicate.append(record_dict[record].id)
                                SeqIO.write(record_dict[record], contig_outfile, "fasta")
                            if strand ==1:
                                start = frame+trans_start*3
                                end = min(seq_len, frame+trans_end*3+3)
                                sequence = record_dict[record].seq[start:end]
                                ORF_seq = SeqRecord(sequence, id=str(record_dict[record].id)+'_'+str(start)+'_'+str(end))
                                SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                            if strand ==-1:
                                start = min(0, seq_len-frame-trans_end*3-3)
                                end = seq_len-frame-trans_start*3
                                sequence = record_dict[record].seq[start:end].reverse_complement()
                                ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
                                SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                            aa_start = trans_end
                        else:
			     aa_start = trans_start + 1


def ORF_find2(s,l):
    noduplicate = []
    record_dict = SeqIO.index(s,"fasta")
    min_seq_len = l
    contig_outfile = open(str(s.split('.')[0])+'orf_contigs.fa', 'w')
    ORF_outfile = open(str(s.split('.')[0])+'orfs.fa', 'w')
    
    for record in record_dict:
        if len(record_dict[record].seq) > (l):
            for strand, nuc in [(+1, record_dict[record].seq), (-1,record_dict[record].seq.reverse_complement())]:
                for frame in range(3):
                    trans = str(nuc[frame:].translate())
                    trans_len = len(trans)
                    seq_len = len(record_dict[record].seq)
                    #print seq_len
                    aa_start = 0
                    while aa_start < trans_len:
                        trans_start = trans.find("M", aa_start)
                        trans_end = trans.find("*", trans_start)
                        #print trans_start
                        #print trans_end
                        #print trans[trans_start:trans_end+1]
                        orf_len = trans_end-trans_start
                        if trans_start == -1:
                            aa_start = trans_len
                        elif trans_end == -1:
                            if trans_len-trans_start > l/3 and trans_len-trans_start < 560:
                                if strand ==1:
                                    start = frame+trans_start*3
                                    end = seq_len
                                    sequence = record_dict[record].seq[start:end]
                                    ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_'+str(start)+'_'+str(end))
                                    print str(record_dict[record].id)+'\t'+str(len(record_dict[record].seq))+'\t'+'\t'+str(end-start)
                                    SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                                    if not record_dict[record].id in noduplicate:
                                            noduplicate.append(record_dict[record].id)
                                            SeqIO.write(record_dict[record], contig_outfile,"fasta")
                                if strand ==-1:
                                    start = 0
                                    end = seq_len-frame-trans_start*3
                                    sequence = record_dict[record].seq[start:end].reverse_complement()
                                    ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
                                    SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                                    print str(record_dict[record].id)+'\t'+str(len(record_dict[record].seq))+'\t'+'\t'+str(end-start)
                                    if not record_dict[record].id in noduplicate:
                                        noduplicate.append(record_dict[record].id)
                                        SeqIO.write(record_dict[record], contig_outfile,"fasta")
                                aa_start = trans_len
                            else:
                                aa_start = trans_len
                        elif orf_len > l/3 and orf_len < 560:
                            if strand ==1:
                                start = frame+trans_start*3
                                end = min(seq_len, frame+trans_end*3+3)
                                sequence = record_dict[record].seq[start:end]
                                ORF_seq = SeqRecord(sequence, id=str(record_dict[record].id)+'_'+str(start)+'_'+str(end-start))
                                print str(record_dict[record].id)+'\t'+str(len(record_dict[record].seq))+'\t'+'\t'+str(end-start)
                                SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                                if not record_dict[record].id in noduplicate:
                                        noduplicate.append(record_dict[record].id)
                                        SeqIO.write(record_dict[record], contig_outfile,"fasta")
                            if strand ==-1:
                                start = min(0, seq_len-frame-trans_end*3-3)
                                end = seq_len-frame-trans_start*3
                                sequence = record_dict[record].seq[start:end].reverse_complement()
                                ORF_seq = SeqRecord(sequence, id=record_dict[record].id+'_RC_'+str(start)+'_'+str(end))
                                print str(record_dict[record].id)+'\t'+str(len(record_dict[record].seq))+'\t'+'\t'+str(end-start)
                                SeqIO.write(ORF_seq, ORF_outfile, "fasta")
                                if not record_dict[record].id in noduplicate:
                                        noduplicate.append(record_dict[record].id)
                                        SeqIO.write(record_dict[record], contig_outfile,"fasta")
                            aa_start = trans_end
                        else:
                             aa_start = trans_start + 1

    

def blast_vsg(s,l):
    filename = s.split('.')[0]
    subprocess.call(['blastn -db tb427_vsgs -query '+str(filename)+'.fa -outfmt 5 -out '+str(filename)+'.xml'], shell=True)
    
def blast_nonvsg(s,l):
     filename = s.split('.')[0]
     subprocess.call(['blastn -db NOTvsgs -query '+str(filename)+'.fa -outfmt 5 -out '+str(filename)+'_nonVSG.xml'], shell=True)
    
def blast_sort(v,n,s):
    #v = vsg xml file
    #n = nonvsg xmlfile
    #s = sequence file, contigs from blast searches
    result_handle = open(v)
    nonVSGresult_handle = open(n)
    blast_records = NCBIXML.parse(result_handle)
    blast_record = blast_records.next()
    record_dict = SeqIO.index(s,"fasta")
    
    outfile = open(s.split('.')[0]+'_VSGs.fa', 'w')
    
    blast_records_nonVSG = NCBIXML.parse(nonVSGresult_handle)
    blast_record_nonVSG = blast_records_nonVSG.next()
    
    hit_list = []
    exclude_list = []
    
    print 'VSG hits! - maybe?'
    
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 1.0e-10: 
                    if not blast_record.query in hit_list:
                        hit_list.append(str(blast_record.query))
                        percent_query_aligned = (100.0 * hsp.identities) / blast_record.query_letters
                        percent_identity = (100.0 * hsp.identities) / alignment.length
                        print str(blast_record.query)+'\t'+str(alignment.title)+'\t'+str(percent_query_aligned)+'\t'+str(percent_identity)+'\t'+str(alignment.length)
    
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
    #print hit_list
    #print exclude_list
    for record in record_dict:
        for VSG in hit_list:
            if VSG not in exclude_list:
                if record_dict[record].id == VSG:
                    SeqIO.write(record_dict[record], outfile, "fasta")
                    
    outfile.close

def cdhitmerge(n):
	sample = n.split('.')[0]
	subprocess.call(['cd-hit-est -i '+n+' -d 0 -o '+sample+'_merged.fa -c 0.9 -n 8 -r 1 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.5'], shell=True)
