#!/usr/bin/env nextflow

process ANNOVAR{
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.csv'
	input:
		tuple val (Sample), path(Vcf)
		val(variant_caller)
	output:
		 tuple val (Sample), path ("${Sample}_${variant_caller}.out.hg19_multianno.csv")
	script:
	"""
	convert2annovar.pl -format vcf4 ${Vcf}  --outfile ${Sample}_${variant_caller}.avinput --withzyg --includeinfo
	table_annovar.pl ${Sample}_${variant_caller}.avinput --out ${Sample}_${variant_caller}.out --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread ${task.cpus} /databases/humandb/ --xreffile /databases/gene_fullxref.txt

	#somaticseqoutput-format_v2_varscan.py ${Sample}.varscan.hg19_multianno.csv ${Sample}.varscan_.csv
	sleep 5s
	"""
}
