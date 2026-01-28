#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LOFREQ {
	tag "${Sample}"
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
	output:
		tuple val (Sample), file ("${Sample}_lofreq_final.hg19_multianno.csv")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${sortedBam}
	samtools sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.01 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample}.lofreq.filtered.vcf --outfile ${Sample}.avinput --withzyg --includeinfo
	perl ${params.annovar2_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}_lofreq_final --remove --protocol refGene,cosmic84,exac03 --operation g,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout ${params.annovar2_path}/humandb/
	"""
}

process FORMAT_LOFREQ {
	tag "${Sample}"
	input:
		tuple val (Sample), file (multiannoFile)
	output:
		tuple val (Sample), file ("${Sample}_lofreq.csv")
	script:
	"""
	lofreqoutput-format.py ${multiannoFile} ${Sample} ./
	"""
}
