# WGS FastQC quality control with optional filtering and trimming options for Illumina paired reads
The steps are:
1. FastQC raw reads (where R1 and R2 files are separate files)
2. MultiQC report of 1.
3. For WGS data: optional deduplication step using clumpify.sh from BBtools to remove technical replicates (for PCR-based library prep such as Illumina Nextera, should not be used for PCR-free approaches such as Illumina Truseq where identical reads represent natural duplicates)
4. BBDUK quality and adapter trim + synthetic contaminants trim

Typical run on UCT hex:
 nextflow run kviljoen/16S_quick_QC --rawReads 'raw_testdata/*{R1,R2}.fastq' --dedup "yes" --decon "no" -profile uct_hex
