#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LOFREQ {
	tag "${Sample}"
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
		file (bed_file)
		path (genome_fasta)
	output:
		tuple val (Sample), file ("${Sample}.lofreq.filtered.vcf")
	script:
	"""
	lofreq viterbi -f ${genome_fasta} -o ${Sample}.lofreq.pre.bam ${sortedBam}
	samtools sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	lofreq call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${genome_fasta} -l ${bed_file} -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	lofreq filter -a 0.01 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}
