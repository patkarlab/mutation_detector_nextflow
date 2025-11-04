#!/usr/bin/env nextflow

process GENERATE_FINAL_BAM{
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	
	input:
		tuple val (Sample), file(alignedRecalibratedBam)
	output:
		tuple val(Sample), file ("*.final.bam"),  file ("*.final.bam.bai")
		
	script:
	"""
	${params.samtools} sort ${alignedRecalibratedBam} > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}