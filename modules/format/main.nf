#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FORMAT_VARDICT {
    tag "${Sample}"
    input:
        tuple val (Sample), file (multiannoFile)
    output:
        tuple val (Sample), path ("${Sample}_vardict.csv")
    script:
    """
    vardictoutput-format.py ${multiannoFile} ${Sample} ./
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


process COMBINE_OUTPUT {
	conda '/home/miniconda3/envs/new_base'
	tag "${Sample}"
	input:
		tuple val (Sample), file(vardict_csv), file(mutect_csv), file(lofreq_csv)
	output:
		tuple val (Sample), path ("${Sample}_combined.csv")
	script:
	"""
	combined_csv-format.py ${vardict_csv} ${mutect_csv} ${lofreq_csv} ./  ${Sample}
	search_KDMmutationDB.py ${Sample}_combined.csv ./ ${Sample}
	"""
}

process MERGE_CSVS {
	tag "${Sample}"
	conda '/home/miniconda3/envs/new_base'
	publishDir "Final_Output/${Sample}/", mode: 'copy', pattern: '*_Final.xlsx'        
	input:
		tuple val (Sample), file(combined_csv), file(counts_file) 
	output:
		file("${Sample}_Final.xlsx")
	script:
	"""
	merge_csv.py -com ${combined_csv} -cov ${counts_file} -o ${Sample}_Final.xlsx
	"""
}

process COMBINE_REPLICATES {
	publishDir "Final_Output/", mode: 'copy', pattern: '*xlsx'
	input:
		val xlsx_tuples
	output:
		file("*.xlsx")
	script:
	def xlsx_files = xlsx_tuples.join(" ")
	"""
	mergeA1B1.py ${xlsx_files}
	"""
}
