#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VARDICT {
	tag "${Sample}"
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
		file (bed_file)
		path (genome_fasta)
		path(ind_files)
	output:
		tuple val (Sample), file ("${Sample}_vardict.vcf")
	script:
	"""
	vardict-java -G ${genome_fasta} -f 0.01 -N ${Sample} -b ${sortedBam} -c 1 -S 2 -E 3 -g 4 ${bed_file} | sed '1d' | teststrandbias.R |  var2vcf_valid.pl -N ${Sample} -E -f 0.01 > ${Sample}_vardict.vcf
	"""
}

