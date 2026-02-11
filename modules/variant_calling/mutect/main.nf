#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MUTECT {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
		file (bed_file)
		path (genome_fasta)
		path(ind_files)
	output:
		tuple val (Sample), file ("${Sample}_split_mutect.vcf")
	script:
	"""
	mv ${genome_fasta}.dict ${genome_fasta.simpleName}.dict
	gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 -R ${genome_fasta} -I:tumor ${sortedBam} -O ${Sample}_mutect.vcf -L ${bed_file}
	gatk LeftAlignAndTrimVariants  -R ${genome_fasta} -V ${Sample}_mutect.vcf --split-multi-allelics true -O ${Sample}_split_mutect.vcf
	"""
}
