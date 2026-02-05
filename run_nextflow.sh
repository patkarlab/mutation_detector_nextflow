#!/usr/bin/bash
##BSUB -J smMIPS_pipeline
##BSUB -n 25
##BSUB -q normal
##-m cn2" to submit jobs to cn2
## or " -m cn3"

##########
#for ENTRY : BEDFILES#
##for AMPLICON [SRSF2/GATA2/U2AF2]: USG
##for AMPLICON [KDM]: ablnew
##for CLL: cll_coordinates_file
##for ALL: ALL_coordinates_byexons
##for MIPS: myeloidrp_coordinates_file
##for DICER: use the second command mentioned below

#nextflow  run main.nf -entry KDM -profile docker -resume -bg

nextflow -c /home/pipelines/mutation_detector_nextflow/nextflow.config run main.nf -entry DICER --sequences /home/pipelines/mutation_detector_nextflow/sequences/ --input /home/pipelines/mutation_detector_nextflow/sample_list.csv -resume -bg
