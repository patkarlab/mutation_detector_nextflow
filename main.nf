#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"mkdir Coverview".execute()


log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}

"""


if (params.help) {
    log.info "----------------------------------------------------------------"
    log.info ""
	log.info "						USAGE                                     "
    log.info ""
    log.info "	nextflow run main.nf [options] -entry [workflow options]"
	log.info ""
	log.info "Options"
	log.info "	--input                 path to csv file containing list of sample IDs."
	log.info "	--sequences             path to directory containing fastq files."
	log.info "	--bedfile               path to the BED4 coordinate bedfile (without extension)"
	log.info "	--genome                path to reference genome"
	log.info "	--adaptors              path to FASTA file of adaptors for trimming"
	log.info ""
	log.info "Workflow Options"
	log.info "	MIPS"
	log.info "	ALL"
	log.info "	CLL"
	log.info "	AMPLICON"
	log.info ""
	log.info "Optional Arguments"
	log.info "	--cpus			Max Number of CPUs"
	log.info ""
}


process trimming_fastq_mcf {
	maxForks 15
	publishDir "${PWD}/${Sample}/processed_reads/", mode: 'copy'
	input:
		val (Sample)
	output:
	 tuple val (Sample), file ("*.fastq")
	script:
	"""
	${params.ea_utils_path}/fastq-mcf -o ${Sample}.R1.trimmed.fastq -o ${Sample}.R2.trimmed.fastq -l 53 -k 0 -q 0 ${params.adaptors} ${params.sequences}/${Sample}*_R1_*.fastq.gz ${params.sequences}/${Sample}*_R2_*.fastq.gz
	"""
}



process gzip{
	publishDir "${PWD}/${Sample}/processed_reads/", mode: 'copy'
	input:
		tuple val (Sample), file(trimmedFiles)
	output:
		val Sample
	script:
	"""
	mkdir "$PWD/${Sample}/Annovar_Modified/" 
	gzip -f ${PWD}/${Sample}/processed_reads/${trimmedFiles[0]}
	gzip -f ${PWD}/${Sample}/processed_reads/${trimmedFiles[1]}
	"""
}


process pair_assembly_pear {
	memory '7.0 GB'
	publishDir "${PWD}/${Sample}/assembled_reads/", mode: 'copy'
	input:
		val Sample
	output:
		tuple val (Sample), file("*") 
	script:
	"""
	${params.pear_path} -f ${PWD}/${Sample}/processed_reads/${Sample}.R1.trimmed.fastq.gz -r ${PWD}/${Sample}/processed_reads/${Sample}.R2.trimmed.fastq.gz -o ${Sample} -n 53 -j 25
	#${params.pear_path} -f ${PWD}/sequences/${Sample}*_R1_*.fastq.gz -r ${PWD}/sequences/${Sample}*_R1_*.fastq.gz -o ${Sample} -n 53 -j 25
	"""
}


process mapping_reads{
	maxForks 15
	publishDir "${PWD}/${Sample}/mapped_reads/", mode: 'copy'
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled[0]} > ${Sample}.sam
	"""
}

process MAP_BWA {
	maxForks 15
	publishDir "${PWD}/${Sample}/mapped_reads/", mode: 'copy'
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled[0]} ${pairAssembled[2]} > ${Sample}.sam
	"""
}

process sam_conversion{
	maxForks 15
	publishDir "$PWD/${Sample}/mapped_reads/", mode: 'copy', pattern: '*.fxd_sorted.bam'
	publishDir "$PWD/${Sample}/mapped_reads/", mode: 'copy', pattern: '*.fxd_sorted.bam.bai'

	input:
		tuple val (Sample), file(samFile)
	output:
		tuple val(Sample), file ("*.fxd_sorted.bam"), file ("*.fxd_sorted.bam.bai")
	
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.picard_path} FixMateInformation I= ${samFile} O= ${Sample}.fxd.sam VALIDATION_STRINGENCY=SILENT
	${params.samtools} view -bT ${params.genome} ${Sample}.fxd.sam > ${Sample}.fxd.bam
	${params.samtools} sort ${Sample}.fxd.bam > ${Sample}.fxd_sorted.bam
	${params.samtools} index ${Sample}.fxd_sorted.bam > ${Sample}.fxd_sorted.bam.bai
	"""
}

process RealignerTargetCreator {
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.intervals'
	
	input:
		tuple val (Sample), file (bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.intervals")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}.intervals
	"""
}

process IndelRealigner{
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.realigned.bam'
	input:
		tuple val(Sample), file (targetIntervals), file(bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.realigned.bam")
	script:
	"""
	echo ${Sample} ${targetIntervals} ${bamFile}
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${targetIntervals} -o ${Sample}.realigned.bam
	"""
}

process BaseRecalibrator{
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.recal_data.table'
	input:
		tuple val (Sample), file (realignedBam)
	output:
		tuple val(Sample), file ("*.recal_data.table")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${realignedBam} -knownSites ${params.site2} -knownSites ${params.site3} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PrintReads{
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.aligned.recalibrated.bam'
	input:
		tuple val (Sample), file (realignedBam), file (recal_dataTable)
	output:
		tuple val (Sample), file ("*.aligned.recalibrated.bam")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample}.aligned.recalibrated.bam
	"""
}

process minimap_getitd {
	publishDir "$PWD/${Sample}/", mode: 'copy', pattern: '*_getitd'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
			val Sample
	output:
			path "*_getitd"
	script:
	"""
	minimap2 -ax sr ${params.genome_minimap_getitd} ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz > ${Sample}.sam
	${params.samtools} view -b -h ${Sample}.sam -o ${Sample}.bam
	${params.samtools} sort ${Sample}.bam -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process generatefinalbam{
	publishDir "$PWD/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.final.bam.bai'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	
	input:
		tuple val (Sample), file(alignedRecalibratedBam)
	output:
		tuple val(Sample), file ("*.final.bam"),  file ("*.final.bam.bai")
		
	script:
	"""
	${params.samtools} sort ${alignedRecalibratedBam} > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process mutect2_run{
	maxForks 10
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.mutect2.vcf'
	
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")

	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${finalBam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -contamination 0.02 -mbq 30

	"""
}

process freebayes_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.freebayes.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")

	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} -t ${params.bedfile}.bed > ${Sample}.freebayes.vcf 	

	"""
}


process vardict_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.vardict.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf

	"""
}

process varscan_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf.gz'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf.gz'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val(Sample), file ("*.varscan_snp.vcf"),  file ("*.varscan_indel.vcf"), file("*.varscan.vcf")
		
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf

	"""
}

process Annovar{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.varscan_.csv'
	input:
		tuple val (Sample), file(varscanVcf)
	output:
		 tuple val (Sample), file ("*.hg19_multianno.csv"), file("*.varscan_.csv")
	script:
	"""
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${varscanVcf}  --outfile ${Sample}.varscan.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.varscan.avinput --out ${Sample}.varscan --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	python3 ${PWD}/scripts/somaticseqoutput-format_v2_varscan.py ${Sample}.varscan.hg19_multianno.csv ${Sample}.varscan_.csv
	#cp ${Sample}.varscan_.csv $PWD/Final_Output/${Sample}/
	sleep 5s
	"""
}

process lofreq_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.lofreq.filtered.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${finalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.01 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf

	"""
}

process strelka_run{
	publishDir "$PWD/${Sample}/variants/strelka", mode: 'copy', pattern: '*.strelka.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		val (Sample)
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --callRegions  ${params.bedfile}.bed.gz --targeted --runDir ${PWD}/${Sample}/variants/strelka/
	${PWD}/${Sample}/variants/strelka/runWorkflow.py -m local -j 20
	gunzip -f ${PWD}/${Sample}/variants/strelka/results/variants/variants.vcf.gz
	mv ${PWD}/${Sample}/variants/strelka/results/variants/variants.vcf $PWD/${Sample}/variants/${Sample}.strelka.vcf
	
	${params.strelka_path}/configureStrelkaSomaticWorkflow.py --normalBam ${params.NA12878_bam}  --tumorBam ${finalBam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ${PWD}/${Sample}/variants/strelka-somatic/
	${PWD}/${Sample}/variants/strelka-somatic/runWorkflow.py -m local -j 20
	
	${params.bcftools_path} concat -a ${PWD}/${Sample}/variants/strelka-somatic/results/variants/somatic.indels.vcf.gz ${PWD}/${Sample}/variants/strelka-somatic/results/variants/somatic.snvs.vcf.gz -o ${PWD}/${Sample}/variants/${Sample}.strelka-somatic.vcf

	"""
}

process somaticSeq_run {
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.somaticseq.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	input:
		tuple val (Sample), file(mutectVcf), file(vardictVcf), file(varscanVcf), file(lofreqVcf)
	output:
		tuple val (Sample), file ("*.somaticseq.vcf"), file("*.hg19_multianno.csv")
	script:
	"""
	somaticseq_parallel.py --output-directory ${PWD}/${Sample}/variants/${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${PWD}/${Sample}/gatk38_processing/${Sample}.final.bam --mutect2-vcf ${PWD}/${Sample}/variants/${Sample}.mutect2.vcf --vardict-vcf ${PWD}/${Sample}/variants/${Sample}.vardict.vcf --varscan-vcf ${PWD}/${Sample}/variants/${Sample}.varscan.vcf --lofreq-vcf ${PWD}/${Sample}/variants/${Sample}.lofreq.vcf	--strelka-vcf ${PWD}/${Sample}/variants/${Sample}.strelka.vcf  --sample-name ${Sample}
	
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sSNV.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sSNV.vcf | sort -k1,1V -k2,2g >> ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz
	
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sINDEL.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sINDEL.vcf | sort -k1,1V -k2,2g >> ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz
	
	${params.bcftools_path} concat -a ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf
	
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	"""
}

process platypus_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.platypus.vcf'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		tuple val(Sample), file ("*.platypus.vcf")
	script:
	"""
	python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${params.bedfile}_regions.txt

	"""
}

process coverage {
	publishDir "$PWD/${Sample}/coverage/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	"""
}

process getITD_run{
	publishDir "$PWD/${Sample}/getITD/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	gunzip -c ${params.sequences}/${Sample}*_R1_*.fastq.gz > ${Sample}_R1_.fastq
	gunzip -c ${params.sequences}/${Sample}*_R2_*.fastq.gz > ${Sample}_R2_.fastq
	
	python3 ${params.getitd_path}/getitd.py ${Sample} ${Sample}_R1_.fastq ${Sample}_R2_.fastq -require_indel_free_primers False -anno ${params.getitd_path}/anno/amplicon_kayser.tsv -reference ${params.getitd_path}/anno/amplicon.txt -forward_primer ${params.forward_primer} -reverse_primer ${params.reverse_primer} -forward_adapter  ATACGAGATCCGTAATCGGGAAGCTGAAG -reverse_adapter ACACGCACGATCCGACGGTAGTGT -nkern 25
	"""
}

process getITD_done{
	input:
		tuple val (Sample), file ("*")
	output:
		val Sample
	script:
	"""
	cp -r ${PWD}/${Sample}/getITD/${Sample}_getitd/ ${PWD}/Final_Output/${Sample}/
	"""

}

process pindel{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*pindel_SI.vcf'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*pindel_SI.vcf'
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*.avinput'
        publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*_pindel.hg19_multianno.csv'
        publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_pindel.hg19_multianno.csv'

	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	sh /home/pipelines/mutation_detector_nextflow/scripts/make_config.sh -s ${Sample}
	${params.pindel_path}/pindel -f ${params.genome} -i $PWD/config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel_path}/pindel2vcf -r ${params.genome} -p ${Sample}_pindel_SI -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_pindel_SI.vcf --outfile ${Sample}_pindel.avinput --withzyg --includeinfo

        perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}_pindel.avinput  ${params.annovarLatest_path}/humandb/ -buildver hg19 -out ${Sample}_pindel  --remove -protocol refGene,cytoBand,cosmic84 --operation g,r,f -nastring '.' --otherinfo --csvout --thread 10 --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	"""
}

process coverview_run {
	executor="local"
	publishDir "$PWD/${Sample}/Coverview/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
	output:
		val (Sample)
	script:
	"""
	if [ ! -d "$PWD/${Sample}/Coverview/" ]; then
		mkdir "$PWD/${Sample}/Coverview/"
	fi
	#mkdir $PWD/${Sample}/Coverview/
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${PWD}/${Sample}/Coverview/${Sample}.coverview
	python3 ${params.coverview_script_path} ${PWD}/${Sample}/Coverview/${Sample}.coverview_regions.txt ${PWD}/${Sample}/Coverview/${Sample}.coverview_regions.csv
	cp ${PWD}/${Sample}/Coverview/${Sample}.coverview_regions.csv ${PWD}/Coverview/${Sample}.coverview_regions.csv
	"""
}

process coverview_report {
	errorStrategy 'ignore'
	executor="local"
	input:
		val (Sample)
	output:
		val Sample
	script:
	"""
	python3.6 ${params.coverview_report_path} ${PWD}/Coverview/ ${PWD}/Final_Output/
	"""
}

process combine_variants{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	
	input:
		tuple val (Sample), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val(Sample), file ("*.combined.vcf"),  file ("*.hg19_multianno.csv")
	script:
	"""
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf > ${Sample}.freebayes.sorted.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf | sort -k1,1V -k2,2g >> ${Sample}.freebayes.sorted.vcf
	
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.platypus.vcf > ${Sample}.platypus.sorted.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.platypus.vcf | sort -k1,1V -k2,2g >> ${Sample}.platypus.sorted.vcf
	
	${params.java_path}/java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample}.freebayes.sorted.vcf --variant ${Sample}.platypus.sorted.vcf -o ${Sample}.combined.vcf -genotypeMergeOptions UNIQUIFY
	
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.combined.vcf  --outfile ${Sample}.combined.avinput --withzyg --includeinfo
	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.combined.avinput --out ${Sample}.combined --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	"""
}

process cava {
	publishDir "$PWD/${Sample}/CAVA/", mode: 'copy'
	
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file(combinedVcf)
	
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i $PWD/${Sample}/variants/${Sample}.somaticseq.vcf -o ${Sample}.somaticseq
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i $PWD/${Sample}/variants/${Sample}.combined.vcf -o ${Sample}.combined
	python3.6 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.combined.txt ${Sample}.cava.csv
	"""
}

process format_somaticseq {
	input:
		tuple val (Sample), file(somaticseqVcf), file (multianno)
	output:
		val Sample
	script:
	"""
	python3.6 ${params.format_somaticseq_script} ${PWD}/${Sample}/ANNOVAR/${Sample}.somaticseq.hg19_multianno.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv
	"""
}

process format_combined {
	input:
		tuple val(Sample), file (combinedVcf),  file (hg19_multianno)
	output:
		val Sample
	script:
	"""
	mkdir -p ${PWD}/${Sample}/Annovar_Modified

	if [ -s ${PWD}/${Sample}/ANNOVAR/${hg19_multianno} ]; then 
		python3.6 ${params.format_combined_script} ${PWD}/${Sample}/ANNOVAR/${hg19_multianno} ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv
	else
		touch ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv
	fi
	"""
}

process merge_csv {
	input:
		tuple val (Sample), file (cava_csv)
	output:
		val Sample
	script:
	"""

	python3.6 ${params.merge_csvs_script} ${Sample} ${PWD}/${Sample}/Annovar_Modified/ ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${PWD}/${Sample}/CAVA/ ${PWD}/${Sample}/Coverview/${Sample}.coverview_regions.csv
	"""
}

process Final_Output {
	input:
		tuple val (Sample), file ("*")
	output:
		val Sample
	script:
	"""
	python3.6 ${params.coveragePlot_script} ${Sample} $PWD/${Sample}/coverage/${Sample}.counts.bed $PWD/${Sample}/coverage/
	cp ${PWD}/${Sample}/coverage/${Sample}.Low_Coverage.png ${PWD}/Final_Output/${Sample}/
	"""
}

process remove_files{
	errorStrategy 'ignore'
	input:
		tuple val (Sample)
	script:
	"""
	rm -rf ${PWD}/${Sample}/
	#rm -rf ${PWD}/Coverview/
	"""
}

workflow MIPS {

    Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_fastq_mcf(samples_ch) | gzip | pair_assembly_pear | mapping_reads | sam_conversion
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	platypus_run(generatefinalbam.out)
	coverage(generatefinalbam.out)
	freebayes_run(generatefinalbam.out)
	mutect2_run(generatefinalbam.out)
	vardict_run(generatefinalbam.out)
	varscan_run(generatefinalbam.out)
	lofreq_run(generatefinalbam.out)
	strelka_run(generatefinalbam.out)
	pindel(generatefinalbam.out)
	somaticSeq_run(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out)))))
	getITD_run(generatefinalbam.out)
	getITD_done(getITD_run.out)
	
	coverview_run(generatefinalbam.out)
	coverview_report(coverview_run.out.toList())
	
	combine_variants(freebayes_run.out.join(platypus_run.out))
	cava(somaticSeq_run.out.join(combine_variants.out))
	format_somaticseq(somaticSeq_run.out)
	format_combined(combine_variants.out)
	
	merge_csv(format_somaticseq.out.join(format_combined.out.join(cava.out)))
	Final_Output(coverage.out)
	remove_files(merge_csv.out.join(coverview_run.out.join(Final_Output.out.join(getITD_done.out))))
}

workflow ALL {
	
    Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_fastq_mcf(samples_ch) | gzip | pair_assembly_pear | mapping_reads | sam_conversion
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	platypus_run(generatefinalbam.out)
	coverage(generatefinalbam.out)
	freebayes_run(generatefinalbam.out)
	mutect2_run(generatefinalbam.out)
	vardict_run(generatefinalbam.out)
	varscan_run(generatefinalbam.out)
	lofreq_run(generatefinalbam.out)
	strelka_run(generatefinalbam.out)
	somaticSeq_run(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out)))))
	
	coverview_run(generatefinalbam.out)
	coverview_report(coverview_run.out.toList())
	
	combine_variants(freebayes_run.out.join(platypus_run.out))
	cava(somaticSeq_run.out.join(combine_variants.out))
	format_somaticseq(somaticSeq_run.out)
	format_combined(combine_variants.out)
	
	merge_csv(format_somaticseq.out.join(format_combined.out.join(cava.out)))
	Final_Output(coverage.out)
	remove_files(merge_csv.out.join(coverview_run.out.join(Final_Output.out)))
}

workflow CLL {
	
    Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_fastq_mcf(samples_ch) | gzip | pair_assembly_pear | mapping_reads | sam_conversion
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	platypus_run(generatefinalbam.out)
	coverage(generatefinalbam.out)
	freebayes_run(generatefinalbam.out)
	mutect2_run(generatefinalbam.out)
	vardict_run(generatefinalbam.out)
	varscan_run(generatefinalbam.out)
	lofreq_run(generatefinalbam.out)
	strelka_run(generatefinalbam.out)
	somaticSeq_run(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out)))))
	
	coverview_run(generatefinalbam.out)
	coverview_report(coverview_run.out.toList())
	
	combine_variants(freebayes_run.out.join(platypus_run.out))
	cava(somaticSeq_run.out.join(combine_variants.out))
	format_somaticseq(somaticSeq_run.out)
	format_combined(combine_variants.out)
	
	merge_csv(format_somaticseq.out.join(format_combined.out.join(cava.out)))
	Final_Output(coverage.out)
	remove_files(merge_csv.out.join(coverview_run.out.join(Final_Output.out)))
}




process trimming_trimmomatic {
	maxForks 10
	publishDir "$PWD/${Sample}/trimmed", mode: 'copy'
	input:
		val Sample
	output:
		tuple val (Sample), file("*.fq.gz")
	script:
	"""
	trimmomatic PE \
	${params.sequences}/${Sample}*_R1_*.fastq.gz ${params.sequences}/${Sample}*_R2_*.fastq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	"""
}

process pair_assembly_flash {
	maxForks 20
	publishDir "$PWD/${Sample}/flash", mode: 'copy'
	input:
		tuple val (Sample), file(trimmedFiles)
	output:
		tuple val (Sample), file('*')
	script:
	"""
	${params.flash_path} ${trimmedFiles[0]} ${trimmedFiles[2]} --cap-mismatch-quals -O -M 250 -o ${Sample} 
	"""
}

process samtools_sort{
	maxForks 20
	publishDir "$PWD/${Sample}/Alignment", mode: 'copy', pattern: '*.sorted.bam'
	publishDir "$PWD/${Sample}/Alignment", mode: 'copy', pattern: '*.bam.bai'
	publishDir "$PWD/${Sample}/Alignment", mode: 'copy', pattern: '*.bed'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.sorted.bam'
	publishDir "$PWD/Final_Output/${Sample}", mode: 'copy', pattern: '*.bam.bai'

	input:
		tuple val (Sample), file(samFile)
	output:
		tuple val (Sample), file ("*.sorted.bam"), file ("*.sorted.bam.bai")
	
	script:
	"""
	${params.samtools} view -b ${samFile} > ${Sample}.bam
	${params.samtools} sort ${Sample}.bam > ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai
	${params.bedtools} bamtobed -i ${Sample}.sorted.bam > ${Sample}.bed
	"""
}

process vardict_run_amplicon{
	maxForks 15
	publishDir "$PWD/${Sample}/VariantCalling/Vardict/", mode: 'copy'
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
	output:
		tuple val (Sample), file ("*.hg19_multianno.csv")
	script:
	"""
	VarDict -G ${params.genome} -f 0.01 -N ${Sample} -b ${sortedBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | ${params.vardict_path}/teststrandbias.R | ${params.vardict_path}/var2vcf_valid.pl -N ${Sample} -E -f 0.01 > ${Sample}_vardict.vcf
	#cp ${Sample}_vardict.vcf ${PWD}/${Sample}/VariantCalling/
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample}_vardict.vcf --outfile ${Sample}.avinput --withzyg --includeinfo
	perl ${params.annovar2_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}_final --remove --protocol refGene,cosmic84,exac03 --operation g,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout ${params.annovar2_path}/humandb/
	"""
}

process lofreq_run_amplicon{
	maxForks 15
	publishDir "$PWD/${Sample}/VariantCalling/Lofreq/", mode: 'copy'
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
	output:
		tuple val (Sample), file ("*.hg19_multianno.csv")
	
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${sortedBam}
	samtools sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.01 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	#cp ${Sample}.lofreq.filtered.vcf ${PWD}/${Sample}/VariantCalling/
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample}.lofreq.filtered.vcf --outfile ${Sample}.avinput --withzyg --includeinfo
	perl ${params.annovar2_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}_final --remove --protocol refGene,cosmic84,exac03 --operation g,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout ${params.annovar2_path}/humandb/
	"""
}

process mutect_run_amplicon{
	maxForks 15
	publishDir "$PWD/${Sample}/VariantCalling/Mutect2/", mode: 'copy'
	input:
		tuple val (Sample), file (sortedBam), file (sortedBamBai)
	output:
		tuple val (Sample), file ("*.hg19_multianno.csv")
		
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${PWD}/${Sample}/Alignment/${sortedBam} -o ${Sample}_mutect.vcf -L ${params.bedfile}.bed
	#cp ${Sample}_mutect.vcf ${PWD}/${Sample}/VariantCalling/
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample}_mutect.vcf --outfile ${Sample}.avinput --withzyg --includeinfo
	perl ${params.annovar2_path}/table_annovar.pl ${Sample}.avinput --out ${Sample}_final --remove --protocol refGene,cosmic84,exac03 --operation g,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout ${params.annovar2_path}/humandb/
	"""
}

process format_VardictOutput_amplicon{
	publishDir "$PWD/${Sample}/VariantCalling/", mode: 'copy'
	input:
		tuple val (Sample), file (multiannoFile)
	output:
		val (Sample)
	script:
	"""
	python3 ${params.formatVardict_script_path} ${multiannoFile} ${Sample} ${PWD}/${Sample}/VariantCalling/
	"""
}

process format_MutectOutput_amplicon{
	publishDir "$PWD/${Sample}/VariantCalling/", mode: 'copy'
	input:
		tuple val (Sample), file (multiannoFile)
	output:
		val (Sample)
	script:
	"""
	python3 ${params.formatMutect_script_path} ${multiannoFile} ${Sample} ${PWD}/${Sample}/VariantCalling/
	"""
}

process format_LofreqOutput_amplicon{
	maxForks 20
	publishDir "$PWD/${Sample}/VariantCalling/", mode: 'copy'
	input:
		tuple val (Sample), file (multiannoFile)
	output:
		val (Sample)
	script:
	"""
	python3 ${params.formatLofreq_script_path} ${multiannoFile} ${Sample} ${PWD}/${Sample}/VariantCalling/
	"""
}

process combine_output_amplicon{
	maxForks 20
	input:
		val (Sample)
	output:
		val (Sample)
	script:
	"""
	python3 ${params.combinedCsvFormat_script_path} ${PWD}/${Sample}/VariantCalling/${Sample}_vardict.csv ${PWD}/${Sample}/VariantCalling/${Sample}_mutect.csv ${PWD}/${Sample}/VariantCalling/${Sample}_lofreq.csv ${PWD}/${Sample}/VariantCalling/ ${Sample}
	
	python3 ${params.KDMdb_script_path} ${PWD}/${Sample}/VariantCalling/${Sample}_combined.csv ${PWD}/${Sample}/VariantCalling/ ${Sample}
	"""
}

process merge_csvs_amplicon{
	errorStrategy 'ignore'
	maxForks 20
	input:
		val Sample
	output:
		val Sample
	script:
	"""
	python3 ${params.mergeAmpliconCsv_path} ${Sample} ${PWD}
	"""
}

process combine_replicates {
	input:
		val Sample
	output:
		val Sample
	script:
	"""
	${params.merge_A1B1} ${params.input}
	${params.delete_samples} ${params.input}	# remove_files 
	"""
}

workflow AMPLICON {
    Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { Sample }
	
	main:
	trimming_trimmomatic(Sample) | pair_assembly_flash | mapping_reads | samtools_sort
	// trimming_trimmomatic(Sample) | MAP_BWA | samtools_sort
	coverage(samtools_sort.out)
	vardict_run_amplicon(samtools_sort.out) | format_VardictOutput_amplicon
	lofreq_run_amplicon(samtools_sort.out) | format_LofreqOutput_amplicon
	mutect_run_amplicon(samtools_sort.out) | format_MutectOutput_amplicon
	combine_output_amplicon(format_VardictOutput_amplicon.out.join(format_LofreqOutput_amplicon.out.join(format_MutectOutput_amplicon.out)))
	
	merge_csvs_amplicon(combine_output_amplicon.out)
	combine_replicates (merge_csvs_amplicon.out.collect())
	//remove_files(merge_csvs_amplicon.out)	
}

workflow TRIM {
    Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { Sample }
	
	main:
	trimming_trimmomatic(Sample)
	
}

workflow DICER {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_fastq_mcf(samples_ch) | gzip | pair_assembly_pear | mapping_reads | sam_conversion
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	minimap_getitd(samples_ch)
	varscan_run(generatefinalbam.out)
	coverview_run(generatefinalbam.out)
	coverview_report(coverview_run.out.toList())
	Annovar(varscan_run.out)
	remove_files(Annovar.out)

}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
