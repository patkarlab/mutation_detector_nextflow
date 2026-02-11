#!/usr/bin/env nextflow

process REALIGNER_TARGET_CREATOR {
	tag "${Sample}"	
	label 'process_medium'
	input:
		tuple val (Sample), file (bamFile), file(bamBai)
		path (genome_fasta)
		path (ind_files)
		path (known_indels)
	output:
		tuple val(Sample), file ("${Sample}.intervals")
	script:
	"""
	mv ${genome_fasta}.dict ${genome_fasta.simpleName}.dict
	java -Xmx${task.memory.toGiga()}g -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${genome_fasta} -nt ${task.cpus} -I ${bamFile} --known ${known_indels} -o ${Sample}.intervals
	"""
}

process INDEL_REALIGNER{
	tag "${Sample}"
	label 'process_medium'	
	input:
		tuple val(Sample), file (targetIntervals), file(bamFile), file(bamBai)
		path (genome_fasta)
		path (ind_files)
		path (known_indels)
	output:
		tuple val(Sample), file ("${Sample}.realigned.bam")
	script:
	"""
	mv ${genome_fasta}.dict ${genome_fasta.simpleName}.dict
	echo ${Sample} ${targetIntervals} ${bamFile}
	java -Xmx${task.memory.toGiga()}g -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner -R ${genome_fasta} -I ${bamFile} -known ${known_indels} --targetIntervals ${targetIntervals} -o ${Sample}.realigned.bam
	"""
}

process BASE_RECALIBRATOR{
	tag "${Sample}"	
	label 'process_medium'
	input:
		tuple val (Sample), file (realignedBam)
		path (genome_fasta)
		path (ind_files)
		path (known_snps_1)
		path (known_snps_2)
	output:
		tuple val(Sample), file ("${Sample}.recal_data.table")
	script:
	"""
	mv ${genome_fasta}.dict ${genome_fasta.simpleName}.dict
	java -Xmx${task.memory.toGiga()}g -jar /usr/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${genome_fasta} -I ${realignedBam} -knownSites ${known_snps_1} -knownSites ${known_snps_2} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PRINT_READS{
	tag "${Sample}"	
	label 'process_low'
	input:
		tuple val (Sample), file (realignedBam), file (recal_dataTable)
		path (genome_fasta)
		path (ind_files)
	output:
		tuple val (Sample), file ("${Sample}.aligned.recalibrated.bam")
	script:
	"""
	mv ${genome_fasta}.dict ${genome_fasta.simpleName}.dict
	java -Xmx${task.memory.toGiga()}g -jar /usr/GenomeAnalysisTK.jar -T PrintReads -R ${genome_fasta} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample}.aligned.recalibrated.bam
	"""
}
