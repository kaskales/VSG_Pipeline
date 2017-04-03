#!/usr/bin/env python
import subprocess

def QC(n):
    samplename = n.split('.')[0]
    print 'Trimming galore!...'
    subprocess.call(['trim_galore --length 50 --stringency 3 --dont_gzip '+str(n)], shell=True)

    #cutadapt
    print 'cutadapt...'
    subprocess.call(['cutadapt -b ATTTAGGTGACACTATAG -b CTATAGTGTCACCTAAAT  '+str(samplename)+'_trimmed.fq > '+str(samplename)+'_trimmed2.fq '], shell=True)

    #trim_galore again
    print 'Trimming galore!...'
    subprocess.call(['trim_galore --length 50 --stringency 3 --dont_gzip '+str(samplename)+'_trimmed2.fq '], shell=True)
    subprocess.call(['mv '+str(samplename)+'_trimmed2_trimmed.fq '+str(samplename)+'_trimmed3.fq'], shell=True)

    #remove intermediate trimming files
    print 'removing intermediate trimmed files...'
    subprocess.call(['rm '+str(samplename)+'_trimmed2.fq' ], shell=True)
    subprocess.call(['rm '+str(samplename)+'_trimmed.fq' ], shell=True)

# options to run with trinity
#   --min_contig_length
#   --normalize_reads
# Trinity --seqType fq --max_memory 10G --single 15-19-1_test_trimmed3.fq --min_contig_length 900 --normalize_reads
    
def bowtie(f,n,p):
    #p = path to MULTo1.0 folder, default = ~/bioinformatic_resources/
    #n = multo number
    #f = file for fq input, 
    samplename = f.split('.')[0]
    subprocess.call(['bowtie -v 2 -m 1 -p 6 -S -a --strata --best '+str(p)+'MULTo1.0/files/tbb/tb'+str(n)+'/bowtie_indexes/tb'+str(n)+'_genome/tb'+str(n)+'_no_random '+str(samplename)+'.fq '+str(samplename)+'_align.sam'], shell=True)

def count(f,n,p):
    #runs rpkmforgenes.py
    if f[-3:]!='sam':
        print 'Input must be sam file!'
    else:
        samplename = f.split('.')[0]
        tbnumber = str(n)
        subprocess.call(['python '+str(p)+'MULTo1.0/src/rpkmforgenes.py -i '+samplename+'.sam -samu -bedann -a '+str(p)+'MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/chr1.bed -u '+str(p)+'MULTo1.0/files/tbb/tb'+tbnumber+'/MULfiles/tb'+tbnumber+'_20-255/MULTo_files -o '+samplename+'_MULTo.txt'], shell=True)
    
        #analyze MULTo out	
        infile = open(samplename+'_MULTo.txt','r')
        outfile = open(samplename+'_MULTo_analyzed.txt', 'w')
        
        
        infile2 = infile.read()
        infile_split = infile2.split('\n')
        
        FPKM = float(0)
        values = []
        pairs = {}
     
        for line in infile_split:
             if not line.startswith("#"):
                l = line.split("\t")
                if len(l) > 2:
                    FPKM += float(l[2])
           
        for line in infile_split:
            if not line.startswith("#"):	
                l = line.split("\t")
                if len(l) == 3 :	
                    values.append(float(l[2])/FPKM)
                    pairs[l[0]] = float(l[2])/FPKM
                        
                        
        values.sort(reverse=True)
        
        
        
        for value in values:
            if float(value) > 0:
                for key in pairs.keys():
                    if float(value) == float(pairs[key]):
                        outfile.write(str(key)+"\t"+str(100*float(value))+"\n")
        for key in pairs.keys():
            if float(pairs[key])==0:
                    outfile.write(str(key)+"\t"+str(0.0)+"\n")
  

def make_multo(r,n):
	samplename = r.split('.')[0]
	tbnumber = str(n)
	#make multo file heirarchy from given basename
	subprocess.call(['mkdir -p ~/bioinformatics_resources/MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/'], shell=True)

	subprocess.call(['mkdir -p ~/bioinformatics_resources/MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)


	#concat fq and make bed - BEDmaker.py with input file - 
	subprocess.call(['python ~/bioinformatics_resources/BEDmaker.py '+r+' tb'+tbnumber+'.fa chr1.bed'], shell=True)

	#move concat and bedfile into new folders
	subprocess.call(['mv tb'+tbnumber+'.fa ~/bioinformatics_resources/MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)

	subprocess.call(['mv chr1.bed ~/bioinformatics_resources/MULTo1.0/files/tbb/tb'+tbnumber+'/fastaFiles/annotationFiles/'], shell=True)

	#run MULTo 

	subprocess.call(['python ~/bioinformatics_resources/MULTo1.0/src/MULTo1.0.py -s tbb -a tb'+tbnumber+' -v 2'], shell=True)
