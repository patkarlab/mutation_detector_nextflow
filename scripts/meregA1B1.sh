#!/usr/bin/bash

sample_list=$1

awk -F- '{ print $1 }' ${sample_list} | sort | uniq > sample_list.csv.dat
for sampleId in $(cat sample_list.csv.dat)
do
	readarray -t  sample_array <<< $(grep -w "${sampleId}"  ${sample_list})
	sample1=${sample_array[0]} 
	sample2=${sample_array[1]}
	echo ${sample1} ${sample2}
	python3 ${PWD}/scripts/mergeamplicon.py "/home/pipelines/mutation_detector_nextflow/${sample1}/coverage/${sample1}.counts.bed" "/home/pipelines/mutation_detector_nextflow/${sample2}/coverage/${sample2}.counts.bed" "/home/pipelines/mutation_detector_nextflow/${sample1}/VariantCalling/${sample1}_combined.csv"  "/home/pipelines/mutation_detector_nextflow/${sample2}/VariantCalling/${sample2}_combined.csv" "/home/pipelines/mutation_detector_nextflow/Final_Output/${sampleId}.xlsx"
done		
rm sample_list.csv.dat
rm A1B1common.csv
