#!/usr/bin/env nextflow

process MINIMAP_GETITD {
	tag "${Sample}"
	label 'process_high'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
			tuple val(Sample), file(read1), file(read2)
			path (minimap_genome)
	output:
			path "*_getitd"
	script:
	"""
	minimap2 -ax sr -t ${task.cpus} ${minimap_genome} ${read1} ${read2} > ${Sample}.sam
	${params.samtools} view -@ ${task.cpus} -b -h ${Sample}.sam -o ${Sample}.bam
	${params.samtools} sort -@ ${task.cpus} ${Sample}.bam -o ${Sample}.sorted.bam
	${params.samtools} index -@ ${task.cpus} ${Sample}.sorted.bam
	${params.samtools} view -@ ${task.cpus} ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}