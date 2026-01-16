#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VARDICT {
	tag "${Sample}"
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
		file bed_file
		path (GenFile)
        path (GenDir)
	output:
		tuple val (Sample), file ("${Sample}_vardict_final.hg19_multianno.csv")
	script:
	"""
	VarDict -G ${GenFile} -f 0.01 -N ${Sample} -b ${sortedBam} -c 1 -S 2 -E 3 -g 4 ${bed_file}.bed | sed '1d' | bin/teststrandbias.R | bin/var2vcf_valid.pl -N ${Sample} -E -f 0.01 > ${Sample}_vardict.vcf
	#perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample}_vardict.vcf --outfile ${Sample}.avinput --withzyg --includeinfo
	#perl ${params.annovar2_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}_vardict_final --remove --protocol refGene,cosmic84,exac03 --operation g,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout ${params.annovar2_path}/humandb/
	"""
}

process FORMAT_VARDICT {
	tag "${Sample}"
	input:
		tuple val (Sample), file (multiannoFile)
	output:
		tuple val (Sample), path ("${Sample}_vardict.csv")
	script:
	"""
	python3 ${params.formatVardict_script_path} ${multiannoFile} ${Sample} ./
	"""
}
