#!/usr/bin/env nextflow

process VARSCAN{
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
		path (genome_fasta)
		path (ind_files)
	output:
		tuple val(Sample), file("*.varscan.vcf")
	script:
	"""
	samtools mpileup -f ${genome_fasta} ${finalBam} > ${Sample}.mpileup
	java -Xmx${task.memory.toGiga()}g -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	java -Xmx${task.memory.toGiga()}g -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	bcftools index -t ${Sample}.varscan_snp.vcf.gz
	bcftools index -t ${Sample}.varscan_indel.vcf.gz
	bcftools concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}
