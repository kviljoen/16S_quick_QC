# WGS QC
A nextflow-based pipeline for quality assessment and read trimming/filtering. This pipeline accepts raw reads in .fastq format, performs (optional) deduplication, quality filtering, adapter removal and (optional) decontamination

The steps are:
1. FastQC raw reads (where R1 and R2 files are separate files)
2. MultiQC report of 1.
3. For WGS data: optional deduplication step using clumpify.sh from BBtools to remove technical replicates (for PCR-based library prep such as Illumina Nextera, should not be used for PCR-free approaches such as Illumina Truseq where identical reads represent natural duplicates)
4. BBDUK quality and adapter trim + synthetic contaminants (phix) trim


## Basic usage:

    The typical command for running the pipeline is as follows:
    nextflow run kviljoen/6S_quick_QC --rawReads '*_R{1,2}.fastq.gz' --dedup 'yes' --decon 'no' -profile uct_hex

    Mandatory arguments:
      --rawReads			Path to input data (must be surrounded with quotes)
      -profile			  Hardware config to use. uct_hex OR standard
      --dedup      Should read deduplication be performed <'yes', 'no'> (see README for details)
      --decon      Should decontamination be performed <'yes', 'no'> (see README for details)
    BBduk trimming options:
      --qin			Input quality offset: 33 (ASCII+33) or 64 (ASCII+64, default=33
      --kcontaminants		Kmer length used for finding contaminants, default=23	
      --phred			Regions with average quality BELOW this will be trimmed, default=10 
      --minlength		Reads shorter than this after trimming will be discarded, default=60
      --mink			Shorter kmers at read tips to look for, default=11 
      --hdist			Maximum Hamming distance for ref kmers, default=1  
      
    BBwrap parameters for decontamination:	
      --mind			Approximate minimum alignment identity to look for, default=0.95
      --maxindel		Longest indel to look for, default=3
      --bwr			Restrict alignment band to this, default=0.16
	
