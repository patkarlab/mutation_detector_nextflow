#!/usr/bin/env nextflow

process REALIGNER_TARGET_CREATOR {
	input:
		tuple val (Sample), file (bamFile), file(bamBai)
        genome_loc
        known_indels
	output:
		tuple val(Sample), file ("${Sample}.intervals")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${genome_loc} -nt 10 -I ${bamFile} --known ${known_indels} -o ${Sample}.intervals
	"""
}

process INDEL_REALIGNER{
	input:
		tuple val(Sample), file (targetIntervals), file(bamFile), file(bamBai)
        genome_loc
        known_indels
	output:
		tuple val(Sample), file ("${Sample}.realigned.bam")
	script:
	"""
	echo ${Sample} ${targetIntervals} ${bamFile}
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${genome_loc} -I ${bamFile} -known ${known_indels} --targetIntervals ${targetIntervals} -o ${Sample}.realigned.bam
	"""
}

process BASE_RECALIBRATOR{
	input:
		tuple val (Sample), file (realignedBam)
        genome_loc
        known_snps_1
        known_snps_2
	output:
		tuple val(Sample), file ("${Sample}.recal_data.table")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${genome_loc} -I ${realignedBam} -knownSites ${known_snps_1} -knownSites ${known_snps_2} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PRINT_READS{
	input:
		tuple val (Sample), file (realignedBam), file (recal_dataTable)
        genome_loc
	output:
		tuple val (Sample), file ("${Sample}.aligned.recalibrated.bam")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${genome_loc} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample}.aligned.recalibrated.bam
	"""
}