#!/usr/bin/env nextflow

//Validate inputs (dedup has to be specified yes/no)
if (params.dedup == false) { 
	exit 1, "Should Deduplication be performed? Must specify (--dedup) <yes, no>" 
}  
if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
} 

out_dir = file(params.outDir)

out_dir.mkdir()

Channel
    .fromFilePairs( params.rawReads )
    .ifEmpty { error "Cannot find any reads matching: ${params.rawReads}" }
    .set {read_pair_p1; read_pair_p2}


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

Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.
	If the layout is "paired", two FASTQ files are outputted, one for each paired-end.
	If "single", a single FASTQ file will be generated.
	This step is OPTIONAL. De-duplication should be carried on iff you are
    using PCR amplification (in this case identical reads are technical artefacts)
	but not otherwise (identical reads will identify natural duplicates).
    
    process dedup {
	tag { "dedup.${pairId}" }

	input:
	set val(pairId), file(reads) from ReadPairs

	output:
	set val(pairId), file("${pairId}_dedupe_R1.fq"), file("${pairId}_dedupe_R2.fq") into totrim, topublishdedupe
    
    when:
	params.dedup=="yes"
    
	script:
	"""
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	clumpify.sh -Xmx\"\$maxmem\" in1="${reads[0]}" in2="${reads[1]}" out1=${pairId}_dedupe_R1.fq out2=${pairId}_dedupe_R2.fq \
	qin=$params.qin dedupe subs=0 threads=${task.cpus}
	
	"""
}
