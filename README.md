# KDM pipeline

This is a nextflow pipeline for analysing target DNA sequencing data for detecting BCR-ABL kinase domain mutations.

For running this pipeline, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- flash_path = flash executable path
- java_path = directory containing the java executable
- samtools = samtools executable path
- bedtools = bedtools executable path
- lofreq_path = lofreq executable paths
- GATK38_path = path to the GenomeAnalysisTK-3.8 jar file
- trimmomatic, bwa, VarDict

## Usage:

1. Keep the `fastq` files into the `sequences/` folder.

2. Change the `sample_list.csv`. It should have a list of IDs of the samples. 

3.  Run the following script.

```
./run_nextflow.sh > script.log
```
This script contains the nextflow command used to execute the workflow.

```
source activate new_base

nextflow -c $PWD$/nextflow.config run main.nf -entry AMPLICON --bedfile <bedfile of ABL1 exons> \
--sequences $PWD$/sequences/ --input $PWD$/sample_list.csv -resume -bg

```

For running the DICER workflow, use the `-entry DICER` 
