#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_TO_BAM; TRIM_MAPBAM } from './workflows/fastq_bam.nf'
include { COVERAGE; COVERVIEW } from './modules/coverage/main.nf'
include { VARDICT } from './modules/variant_calling/vardict/main.nf'
include { LOFREQ } from './modules/variant_calling/lofreq/main.nf'
include { MUTECT } from './modules/variant_calling/mutect/main.nf'
include { VARSCAN } from './modules/variant_calling/varscan/main.nf'
include { FORMAT_VARDICT; FORMAT_LOFREQ; FORMAT_MUTECT; COMBINE_OUTPUT; MERGE_CSVS; COMBINE_REPLICATES } from './modules/format/main.nf'
include { ANNOVAR as ANNOVAR_VARDICT; ANNOVAR as ANNOVAR_LOFREQ; ANNOVAR as ANNOVAR_MUTECT; ANNOVAR as ANNOVAR_VARSCAN } from './modules/annovar/main.nf'
include { REALIGNER_TARGET_CREATOR; INDEL_REALIGNER; BASE_RECALIBRATOR; PRINT_READS } from './modules/gatk/main.nf'
include { GENERATE_BAM } from './modules/sort/main.nf'
include { MINIMAP; MINIMAP_SORT; MINIMAP_BAMTOFASTQ; GETITD } from './modules/minimap_getitd/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILE PATHS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

bed_file = file("${params.bedfile}", checkIfExists: true )
genome_fasta = file("${params.genome}", checkIfExists: true)
ind_files = file("${params.genome_dir}/${params.ind_files}.*")
dict_file = file("${params.genome_dir}/${params.ind_files}.dict")
vardict = params.vardict
mutect2 = params.mutect2
lofreq  = params.lofreq
varscan = params.varscan
known_snps_1 = file("${params.site1}", checkIfExists: true )
known_snps_2 = file("${params.site3}", checkIfExists: true )
known_indels = file("${params.site2}", checkIfExists: true)
minimap_genome = file("${params.genome_minimap_getitd}", checkIfExists: true )


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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
	kdm_bams_ch = FASTQ_TO_BAM(bam_ch)
	COVERAGE(kdm_bams_ch.final_bams_ch, bed_file )
	COVERVIEW(kdm_bams_ch.final_bams_ch, bed_file )
	VARDICT(kdm_bams_ch.final_bams_ch, bed_file, genome_fasta, ind_files ) 
	ANNOVAR_VARDICT(VARDICT.out, vardict) 
	FORMAT_VARDICT(ANNOVAR_VARDICT.out)
	LOFREQ(kdm_bams_ch.final_bams_ch, bed_file, genome_fasta )  
	ANNOVAR_LOFREQ(LOFREQ.out, lofreq) 
	FORMAT_LOFREQ(ANNOVAR_LOFREQ.out)
	MUTECT(kdm_bams_ch.final_bams_ch, bed_file, genome_fasta, ind_files ) 
	ANNOVAR_MUTECT(MUTECT.out, mutect2) 
	FORMAT_MUTECT(ANNOVAR_MUTECT.out)
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
	dicer_bams_ch = TRIM_MAPBAM(bam_ch)
	REALIGNER_TARGET_CREATOR(dicer_bams_ch.final_bams_ch, genome_fasta, ind_files, known_indels)
	INDEL_REALIGNER(REALIGNER_TARGET_CREATOR.out.join(dicer_bams_ch.final_bams_ch), genome_fasta, ind_files, known_indels)
	BASE_RECALIBRATOR(INDEL_REALIGNER.out, genome_fasta, ind_files, known_snps_1, known_snps_2)
	PRINT_READS(INDEL_REALIGNER.out.join(BASE_RECALIBRATOR.out), genome_fasta, ind_files) | GENERATE_BAM
	VARSCAN(GENERATE_BAM.out, genome_fasta, ind_files)
	ANNOVAR_VARSCAN(VARSCAN.out, varscan)
	MINIMAP(bam_ch, minimap_genome)
	MINIMAP_SORT(MINIMAP.out)
	MINIMAP_BAMTOFASTQ(MINIMAP_SORT.out)
	GETITD(MINIMAP_BAMTOFASTQ.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	println "Total time taken: ${workflow.duration}"
}
