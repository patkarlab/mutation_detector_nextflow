#!/usr/bin/env nextflow

process MINIMAP {
	tag "${Sample}"
	label 'process_medium'
	//publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val(Sample), file(read1), file(read2)
		path (minimap_genome)
	output:
		tuple val(Sample), path("${Sample}.sam")
	script:
	"""
	minimap2 -ax sr -t ${task.cpus} ${minimap_genome} ${read1} ${read2} > ${Sample}.sam
	"""
}

process MINIMAP_SORT {
	tag "${Sample}"
	label 'process_medium'
	//publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val(Sample), file(sam)
	output:
		tuple val(Sample), path("${Sample}.chr13.bam")
	script:
	"""
	samtools view -@ ${task.cpus} -b -h ${sam} -o ${Sample}.bam
	samtools sort -@ ${task.cpus} ${Sample}.bam -o ${Sample}.sorted.bam
	samtools index ${Sample}.sorted.bam
	samtools view -@ ${task.cpus} ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	"""
}

process MINIMAP_BAMTOFASTQ {
	tag "${Sample}"
	label 'process_medium'
	//publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val(Sample), file(bam)
	output:
		tuple val(Sample), path("${Sample}_chr13.fastq")
	script:
	"""
	bedtools bamtofastq -i ${bam} -fq ${Sample}_chr13.fastq
	"""
}

process GETITD {
	tag "${Sample}"
	label 'process_high'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val(Sample), path(chr13_fastq)
	output:
		path "*_getitd"
	script:
	"""
	python3 /opt/getitd/getitd.py -reference /opt/getitd/anno/amplicon.txt -anno /opt/getitd/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern ${task.cpus} ${Sample} ${chr13_fastq}
	"""
}