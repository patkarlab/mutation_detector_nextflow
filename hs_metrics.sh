#! /usr/bin/bash

for i in `cat sample_list.csv` 
do 
	#java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I=/home/pipelines/mutation_detector_nextflow/Final_Output/$i/$i"_final.bam" O=/home/pipelines/mutation_detector_nextflow/Final_Output/$i".hsmetrics_myo.txt" R=/home/reference_genomes/hg19_broad/hg19_all.fasta BAIT_INTERVALS=/home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd.interval_list TARGET_INTERVALS=/home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd.interval_list VALIDATION_STRINGENCY=LENIENT

	#echo -ne $i'\t'; grep -v '#' ${PWD}/Final_Output/$i/$i"_hsmetrics.txt" | awk 'BEGIN{FS="\t"; OFS="\t"}NR==2{ print $7,$8}'
	echo -ne ${i}'\t'; grep -v '#' ${PWD}/Final_Output/${i}".hsmetrics_myo.txt" | awk 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print $7,$8}'
	#echo -ne $i'\t'; grep -v '#' ${PWD}/Final_Output/$i/${i}"_uncoll_hsmetrics.txt" | awk 'BEGIN{FS="\t"; OFS="\t"}NR==2{ print $7,$8}'
	#echo -ne $i'\t'; grep -v '#' ${PWD}/Final_Output/$i/${i}"_uncoll_hsmetrics.txt" | awk 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print $7,$8}'
done
