#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MUTECT {
	tag "${Sample}"
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
	output:
		tuple val (Sample), file ("${Sample}_mutect_final.hg19_multianno.csv")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${sortedBam} -o ${Sample}_mutect.vcf -L ${params.bedfile}.bed
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample}_mutect.vcf --outfile ${Sample}.avinput --withzyg --includeinfo
	perl ${params.annovar2_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}_mutect_final --remove --protocol refGene,cosmic84,exac03 --operation g,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout ${params.annovar2_path}/humandb/
	"""
}

process FORMAT_MUTECT {
	tag "${Sample}"
	input:
		tuple val (Sample), file (multiannoFile)
	output:
		tuple val (Sample), file ("${Sample}_mutect.csv")
	script:
	"""
	mutectoutput-format.py ${multiannoFile} ${Sample} ./
	"""
}
