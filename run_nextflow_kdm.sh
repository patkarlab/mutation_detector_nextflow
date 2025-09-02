#!/usr/bin/bash

nextflow -c /home/pipelines/mutation_detector_nextflow/nextflow.config run kdm.nf -entry KDM --bedfile /home/pipelines/mutation_detector_nextflow/bedfile/ablnew --sequences /home/pipelines/mutation_detector_nextflow/sequences/ --input /home/pipelines/mutation_detector_nextflow/sample_list.csv  --adaptors /home/programs/Trimmomatic-0.39/adapters/TruSeq2-PE.fa --output /home/pipelines/mutation_detector_nextflow/Final_Output  -resume -bg --mergeA1B1 /home/pipelines/mutation_detector_nextflow/scripts/mergeA1B1.py 

