README for bulk_RNAseq directory on Github

Generation of data from bulk RNAseq (known genes, not NTRs).

bulk_RNAseq_fastq_to_combinedcounts_code.txt #code for:
	taking raw fastq files, merging lanes,
	aligning to reference genome,
	QC metrics,
	generation of BAM files,
	and counting reads.
output file is combined counts from GEO repository. 


bulk_RNAseq_combinedcounts_to_analysis.R #code for:
	filtering,
	normalization,
	SAS-based statistical differential expression analysis, and
	FDR correction.
	Input file is combined_counts file from GEO repository.

