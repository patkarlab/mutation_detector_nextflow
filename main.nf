#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

// file paths
bed_file = file("${params.bedfile}.bed", checkIfExists: true )
genome_loc = file("${params.genome}", checkIfExists: true )
known_indels = file("${params.site1}", checkIfExists: true )
known_snps_1 = file("${params.site2}", checkIfExists: true )
known_snps_2 = file("${params.site3}", checkIfExists: true )
minimap_genome = file("${params.genome_minimap_getitd}", checkIfExists: true )

include { FASTQTOBAM } from '/home/pipelines/NextSeq_mutation_detector_leukemia/modules/processes.nf'
include { COVERAGE } from './modules/coverage.nf'
include { VARDICT; FORMAT_VARDICT } from './modules/vardict.nf'
include { LOFREQ; FORMAT_LOFREQ } from './modules/lofreq.nf'
include { MUTECT; FORMAT_MUTECT } from './modules/mutect.nf'
include { COMBINE_OUTPUT; MERGE_CSVS; COMBINE_REPLICATES } from './modules/format.nf'

workflow KDM {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { bam_ch }
	main:
	final_bams_ch = FASTQTOBAM(bam_ch)
	COVERAGE(final_bams_ch)
	VARDICT(final_bams_ch) | FORMAT_VARDICT
	LOFREQ(final_bams_ch) | FORMAT_LOFREQ
	MUTECT(final_bams_ch) | FORMAT_MUTECT
	COMBINE_OUTPUT(FORMAT_VARDICT.out.join(FORMAT_MUTECT.out.join(FORMAT_LOFREQ.out)))
	MERGE_CSVS(COMBINE_OUTPUT.out.join(COVERAGE.out.cov))
	COMBINE_REPLICATES(MERGE_CSVS.out.collect())
}

workflow DICER {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { bam_ch }
	main:
	final_bams_ch = FASTQTOBAM(bam_ch)	
	REALIGNER_TARGET_CREATOR(final_bams_ch, genome_loc, known_indels)
	INDEL_REALIGNER(REALIGNER_TARGET_CREATOR.out.join(final_bams_ch), genome_loc, known_indels)
	BASE_RECALIBRATOR(INDEL_REALIGNER.out, genome_loc, known_snps_1, known_snps_2)
	PRINT_READS(INDEL_REALIGNER.out.join(BASE_RECALIBRATOR.out), genome_loc) | GENERATE_FINAL_BAM
	MINIMAP_GETITD(bam_ch)
	VARSCAN(GENERATE_FINAL_BAM.out)
	COVERVIEW(GENERATE_FINAL_BAM.out)
	COVERVIEW_REPORT(COVERVIEW.out.toList())
	ANNOVAR(VARSCAN.out)
	REMOVE_FILES(ANNOVAR.out)
}


workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	println "Total time taken: ${workflow.duration}"
}