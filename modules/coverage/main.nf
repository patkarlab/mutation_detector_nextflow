#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process COVERAGE {
	tag "${Sample}"
	//publishDir "${params.outdir}/${Sample}", mode: 'copy', pattern: '*.bam*', includeInputs:true
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
		file (bed_file)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed"), emit: cov
		tuple val (Sample), path ("${Sample}_final.bam")
		tuple val (Sample), path ("${Sample}_final.bam.bai")
	script:
	"""
	bedtools coverage -counts -a ${bed_file} -b ${Sample}_final.bam > ${Sample}.counts.bed
	"""
}

process COVERVIEW {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.coverview_regions.csv'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
		file(bed_file)
	output:
		tuple val (Sample), file("${Sample}.coverview_regions.csv")
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${bed_file} -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}

