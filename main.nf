#!/usr/bin/env nextflow

//Validate inputs (dedup has to be specified yes/no)
if (params.dedup == false) { 
	exit 1, "Should Deduplication be performed? Must specify (--dedup) <yes, no>" 
}  
if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
} 
if (params.decon == false) { 
	exit 1, "Should Decontamination be performed? Must specify (--decon) <yes, no>" 
}  

out_dir = file(params.outDir)

out_dir.mkdir()

Channel
    .fromFilePairs( params.rawReads )
    .ifEmpty { error "Cannot find any reads matching: ${params.rawReads}" }
    .set {ReadPairsToQual; ReadPairs}


process runFastQC{
    tag { "${params.projectName}.rFQC.${pairId}" }
    publishDir "${out_dir}/${pairId}", mode: 'copy', overwrite: false

    input:
        set pairId, file(in_fastq) from ReadPairsToQual

    output:
        file("${pairId}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${sample}_fastqc
    fastqc --outdir ${pairId}_fastqc \
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

/* 
 *	Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.
 *	Two FASTQ files are outputted, one for each paired-end.
 *	If "single", a single FASTQ file will be generated.
 *	This step is OPTIONAL. De-duplication should be carried on iff you are
 *    	using PCR amplification (in this case identical reads are technical artefacts)
 *	but not otherwise (identical reads will identify natural duplicates).
 */  
 
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

/*
 *
 * Step 3: BBDUK: trim + filter (run per sample)
 *
 */

process bbduk {
	tag{ "bbduk.${pairId}" }
	
	input:
	set val(pairId), file("${pairId}_dedupe_R1.fq"), file("${pairId}_dedupe_R2.fq") from totrim
	file adapters from params.adapters
	file artifacts from params.artifacts
	file phix174ill from params.phix174ill

	output:
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq"), file("${pairId}_trimmed_singletons.fq") into todecontaminate
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq") into filteredReadsforQC

	script:
	"""	
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	#Quality and adapter trim:
	bbduk.sh -Xmx\"\$maxmem\" in=${pairId}_dedupe_R1.fq in2=${pairId}_dedupe_R2.fq out=${pairId}_trimmed_R1_tmp.fq \
	out2=${pairId}_trimmed_R2_tmp.fq outs=${pairId}_trimmed_singletons_tmp.fq ktrim=r \
	k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred \
	minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe 
	
	#Synthetic contaminants trim:
	bbduk.sh -Xmx\"\$maxmem\" in=${pairId}_trimmed_R1_tmp.fq in2=${pairId}_trimmed_R2_tmp.fq \
	out=${pairId}_trimmed_R1.fq out2=${pairId}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts \
	qin=$params.qin threads=${task.cpus} 
	#Synthetic contaminants trim for singleton reads:
	bbduk.sh -Xmx\"\$maxmem\" in=${pairId}_trimmed_singletons_tmp.fq out=${pairId}_trimmed_singletons.fq \
	k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=${task.cpus}
	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${pairId}_trimmed*_tmp.fq 
	"""
}


/*
 *
 * Step 4: FastQC post-filter and -trim (run per sample)
 *
 */

process runFastQC_postfilterandtrim {
    tag { "rFQC_post_FT.${pairId}" }
    publishDir "${params.outdir}/FastQC_post_filter_trim", mode: "copy", overwrite: false

    input:
    	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq") from filteredReadsforQC

    output:
        file("${pairId}_fastqc_postfiltertrim/*.zip") into fastqc_files_2

    """
    mkdir ${pairId}_fastqc_postfiltertrim
    fastqc --outdir ${pairId}_fastqc_postfiltertrim \
    ${pairId}_trimmed_R1.fq \
    ${pairId}_trimmed_R2.fq
    """
}

process runMultiQC_postfilterandtrim {
    tag { "rMQC_post_FT" }
    publishDir "${params.outdir}/FastQC_post_filter_trim", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files_2.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

/*
 *
 * Step 5: Decontamination (run per sample)
 *
 */

process decontaminate {
	tag{ "decon.${pairId}" }
	publishDir  "${params.outdir}/decontaminate" , mode: 'copy', pattern: "*_clean.fq.gz", overwrite: false
	cache 'deep'
	
	refForeignGenome_ref = file(params.refForeignGenome, type: 'dir')

	input:
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq"), file("${pairId}_trimmed_singletons.fq") from todecontaminate
	file refForeignGenome from refForeignGenome_ref
	
	output:
	file "*_clean.fq.gz"
	set val(pairId), file("${pairId}_cont.fq") into topublishdecontaminate
	
	script:
	"""
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	#Decontaminate from foreign genomes
	bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=${pairId}_trimmed_R1.fq,${pairId}_trimmed_singletons.fq in2=${pairId}_trimmed_R2.fq,null \
	outu=${pairId}_clean.fq outm=${pairId}_cont.fq minid=$params.mind \
	maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred \
	path=$refForeignGenome qin=$params.qin threads=${task.cpus} untrim quickmatch fast
	
	gzip -c ${pairId}_clean.fq > ${pairId}_clean.fq.gz
	"""
}

