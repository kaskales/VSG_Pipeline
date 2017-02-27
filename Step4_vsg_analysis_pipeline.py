#put in cloud. also put BEDmaker.py! usage is MULTo_pipeline_for_cloud.py inputfile.fa tbnumber (ie, 1293 for '-a tb1293' in MULTo run)
#make two steps - one for trinity file, one for analysis from raw data

import argparse
import subprocess
import vsganalysisfunctions as vsg

parser = argparse.ArgumentParser()

parser.add_argument('-i', nargs= '+', help='list of files - FASTQ, or .sam is bowtie and QC is turned off', action="store", dest="i", required=True)
parser.add_argument('-r', help='reference assembly - FASTA', action="store", dest="r")
parser.add_argument('-n', help='assembly number for MULTo, ie, 1293 for \'-a tb1293\' in MULTo run, default = 4', action ="store", dest = "n", default=4)
parser.add_argument('-noqc', help='QC will not be run, default is on', action ='store_false', dest = 'qc' )
parser.add_argument('-p', help='path to MULTo1.0 folder. default is ~/bioinformatics_resources/', action="store", dest='p', default='~/bioinformatics_resources/')
group = parser.add_mutually_exclusive_group()
group.add_argument('-m_only', help='only run rpkmforgenes/multo scripts, input is samfile', action ='store_true', dest = 'm_only')
group.add_argument('-qc_only', help= 'only run QC, default is off', action= 'store_true', dest='qc_only')


arguments = parser.parse_args()

n = arguments.n
p = arguments.p
qc_only = arguments.qc_only
m_only = arguments.m_only

for file in arguments.i:
	if qc_only==True:
		vsg.QC(file)
		
	elif m_only==True:
		if file.split('.')[1]!='sam':
			print 'Input must be sam file!'
		else:
			vsg.count(file,n,p)
	else:
		if arguments.qc==True:
			vsg.QC(file)
			file = file.split('.')[0]+'_trimmed3.fq'
			vsg.bowtie(file,n,p)
			file =file.split('.')[0]+'_align.sam'
		vsg.count(file,n,p)
		
