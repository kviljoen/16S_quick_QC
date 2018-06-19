#!/usr/bin/env nextflow

raw_reads = params.rawReads
out_dir = file(params.outDir)

out_dir.mkdir()

read_pair = Channel.fromFilePairs("${raw_reads}/*stripped_again_R[1,2].fastq", type: 'file')

read_pair.into { read_pair_p1; read_pair_p2 }

process runFastQC{
    tag { "${params.projectName}.rFQC.${sample}" }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: false

    input:
        set sample, file(in_fastq) from read_pair_p1

    output:
        file("${sample}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${sample}_fastqc
    fastqc --outdir ${sample}_fastqc \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
    """
}

process runMultiQC{
    tag { "${params.projectName}.rMQC" }
    publishDir "${out_dir}/", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}
