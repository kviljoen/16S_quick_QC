# 16S_quick_QC
Just a quick QC test running fastqc and multiqc

Typical run on UCT hex:
 nextflow run kviljoen/16S_quick_QC --rawReads 'raw_testdata/*{R1,R2}.fastq' --outdir /path/to/outdir/ -profile hex
