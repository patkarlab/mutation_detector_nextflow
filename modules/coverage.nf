#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process COVERAGE {
	tag "${Sample}"
	publishDir "Final_Output/${Sample}", mode: 'copy', pattern: '*.bam*', includeInputs:true
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed"), emit: cov
		tuple val (Sample), path ("${Sample}_final.bam")
		tuple val (Sample), path ("${Sample}_final.bam.bai")
	script:
	"""
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}_final.bam > ${Sample}.counts.bed
	"""
}

process COVERVIEW {
	executor="local"
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
	output:
		tuple val (Sample), ${Sample}.coverview_regions.csv
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}

process COVERVIEW_REPORT {
	errorStrategy 'ignore'
	executor="local"
	input:
		val (Sample)
	output:
		val Sample
	script:
	"""
	python3.6 ${params.coverview_report_path} ${PWD}/Coverview/ ${PWD}/Final_Output/
	"""
}