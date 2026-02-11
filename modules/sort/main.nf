#!/usr/bin/env nextflow

process GENERATE_BAM{
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	
	input:
		tuple val (Sample), file(alignedRecalibratedBam)
	output:
		tuple val(Sample), file ("*.final.bam"),  file ("*.final.bam.bai")
		
	script:
	"""
	samtools sort ${alignedRecalibratedBam} > ${Sample}.final.bam
	samtools index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}
