#!/usr/bin/nextflow
nextflow.enable.dsl=2

// file paths
adapt = file("${params.adaptors}", checkIfExists: true )
seq_locatn = file("${params.sequences}", checkIfExists: true )
genome_loc = file("${params.genome}", checkIfExists: true)
known_SNPs = file("${params.site1}", checkIfExists: true)
known_SNPs_index = file("${params.site1_idx}", checkIfExists: true)
known_INDELS = file("${params.site2}", checkIfExists: true)
known_INDELS_index = file("${params.site2_idx}", checkIfExists: true)

genome_dir = file("${genome_loc.parent}", checkIfExists: true)
genome_fasta = file("${genome_loc.name}")

include { TRIM; MAPBAM; MARK_DUPS; BQSR; APPLY_BQSR; ALIGNMENT_METRICS; INSERT_SIZE_METRICS} from '../modules/bam.nf'

workflow FASTQ_TO_BAM{
	Channel
		.fromPath(params.input)
		.splitCsv(header: false)
		.map { row -> row[0] }
		.map { sample ->
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz")
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz")
			tuple(sample, r1, r2)
		}
		.set { samples_ch }
	TRIM(samples_ch, adapt)
	MAPBAM(TRIM.out, genome_dir, genome_fasta)
	MARK_DUPS(MAPBAM.out)
	BQSR(MARK_DUPS.out, genome_dir, genome_fasta, known_SNPs, known_SNPs_index, known_INDELS, known_INDELS_index)
	APPLY_BQSR(MARK_DUPS.out.join(BQSR.out), genome_dir, genome_fasta)
	ALIGNMENT_METRICS(APPLY_BQSR.out, genome_dir, genome_fasta)
	INSERT_SIZE_METRICS(APPLY_BQSR.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the ${params.output} directory \n" : "Oops .. something went wrong" )
}