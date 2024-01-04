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

source activate new_base
nextflow -c /home/pipelines/mutation_detector_nextflow/nextflow.config run main.nf -entry AMPLICON --bedfile /home/pipelines/mutation_detector_nextflow/bedfile/ablnew --sequences /home/pipelines/mutation_detector_nextflow/sequences/ --input /home/pipelines/mutation_detector_nextflow/sample_list.csv -resume -bg

#For AMPLICON: After this pipeline, trigger the following command 
#/home/pipelines/mutation_detector_nextflow/scripts/meregA1B1.sh /home/pipelines/mutation_detector_nextflow/sample_list.csv
