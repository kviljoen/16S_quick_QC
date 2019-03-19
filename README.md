# FastQC quality control with optional filtering and trimming options
The steps are:
1. FastQC raw reads
2. MultiQC report of 1.
3. Optional deduplication step using clumpify.sh from BBtools to remove technical replicates (for PCR-based library prep such as Illumina Nextera, should not be used for PCR-free approaches such as Illumina Truseq where identical reads represent natural duplicates)

Typical run on UCT hex:
 nextflow run kviljoen/16S_quick_QC --rawReads 'raw_testdata/*{R1,R2}.fastq' -profile uct_hex
