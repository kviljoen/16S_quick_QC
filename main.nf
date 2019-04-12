#!/usr/bin/env nextflow

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info "==================================="
log.info "           FASTQC_BBDUK            "
log.info "==================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.rawReads
summary['OS']		= System.getProperty("os.name")
summary['OS.arch']	= System.getProperty("os.arch") 
summary['OS.version']	= System.getProperty("os.version")
summary['javaversion'] = System.getProperty("java.version") //Java Runtime Environment version
summary['javaVMname'] = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
summary['javaVMVersion'] = System.getProperty("java.vm.version") //Java Virtual Machine implementation version
//Gets starting time		
sysdate = new java.util.Date() 
summary['User']		= System.getProperty("user.name") //User's account name
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Adapter ref']    = params.adapters
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


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

out_dir = file(params.outdir)

out_dir.mkdir()

Channel
    .fromFilePairs( params.rawReads )
    .ifEmpty { error "Cannot find any reads matching: ${params.rawReads}" }
    .into {ReadPairsToQual; ReadPairs}


process runFastQC{
    tag { "${params.projectName}.rFQC.${pairId}" }
    publishDir "${out_dir}/${pairId}", mode: 'copy', overwrite: false

    input:
        set pairId, file(in_fastq) from ReadPairsToQual

    output:
        file("${pairId}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${pairId}_fastqc
    fastqc --outdir ${pairId}_fastqc \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
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
	maxmem=\$(echo ${task.memory} | sed 's/.GB//g')
	maxmem_java=\$((\$maxmem - 4))
	
	
	clumpify.sh -Xmx40G -Xms2G in1="${reads[0]}" in2="${reads[1]}" out1=${pairId}_dedupe_R1.fq out2=${pairId}_dedupe_R2.fq \
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
	publishDir "${params.outdir}/BBDUK", mode: "copy"

	//bbduk reference files
	adapters_ref = file(params.adapters)
	artifacts_ref = file(params.artifacts)
	phix174ill_ref = file(params.phix174ill)
	
	input:
	set val(pairId), file("${pairId}_dedupe_R1.fq"), file("${pairId}_dedupe_R2.fq") from totrim
	file adapters from adapters_ref
	file artifacts from artifacts_ref
	file phix174ill from phix174ill_ref

	output:
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq"), file("${pairId}_trimmed_singletons.fq") into todecontaminate
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq") into filteredReadsforQC
	set val(pairId), file("${pairId}.stats.txt") into multiqc_bbduk_stats
	
	script:
	"""	
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	#Quality and adapter trim:
	bbduk.sh -Xmx40G -Xms2G in=${pairId}_dedupe_R1.fq in2=${pairId}_dedupe_R2.fq out=${pairId}_trimmed_R1_tmp.fq \
	out2=${pairId}_trimmed_R2_tmp.fq outs=${pairId}_trimmed_singletons_tmp.fq \
	stats=${pairId}.stats.txt \
	ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred \
	minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe 
	
	#Synthetic contaminants trim:
	bbduk.sh -Xmx40G -Xms2G in=${pairId}_trimmed_R1_tmp.fq in2=${pairId}_trimmed_R2_tmp.fq \
	out=${pairId}_trimmed_R1.fq out2=${pairId}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts \
	qin=$params.qin threads=${task.cpus} 
	#Synthetic contaminants trim for singleton reads:
	bbduk.sh -Xmx40G -Xms2G in=${pairId}_trimmed_singletons_tmp.fq out=${pairId}_trimmed_singletons.fq \
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
	
	when:
	decon=="yes"
	
	script:
	"""
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	#Decontaminate from foreign genomes
	bbwrap.sh  -Xmx8G mapper=bbmap append=t in1=${pairId}_trimmed_R1.fq,${pairId}_trimmed_singletons.fq in2=${pairId}_trimmed_R2.fq,null \
	outu=${pairId}_clean.fq outm=${pairId}_cont.fq minid=$params.mind \
	maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred \
	path=$refForeignGenome qin=$params.qin threads=${task.cpus} untrim quickmatch fast
	
	gzip -c ${pairId}_clean.fq > ${pairId}_clean.fq.gz
	"""
}

/*
 *
 * Step : MultiQC report
 *
 */


process multiqc {
    tag { "multiqc" }
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    file('*') from fastqc_files.collect().ifEmpty([])
    file('*') from fastqc_files_2.collect().ifEmpty([])
    file('*') from multiqc_bbduk_stats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc -m fastqc -m bbmap .
    """
}
