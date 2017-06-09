# VSG_Pipeline

Organization Structure

To start:

Current Directory{
	Pipeline.py
	VSGFunctions.py
	[SequencingFile1].fastq
	[SequencingFile2].fastq
}

When Finished:

Current Directory{
	Pipeline.py
	VSGFunctions.py
	[SequencingFile1].fastq
	[SequencingFile2].fastq
	[Y-M-D-H_M]-[OptionalDescriptiveHeaderNames]{
	[SequencingFile1]{
			[SequencingFile1]_trimmed2.fq
		}
	[SequencingFile2]{
			[SequencingFile2]_trimmed2.fq
		}
	}
}