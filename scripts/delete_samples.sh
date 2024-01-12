#!/usr/bin/bash

sample_list=$1
for i in `cat ${sample_list}`
do 
	rm -r /home/pipelines/mutation_detector_nextflow/${i}
done
