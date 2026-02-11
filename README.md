# This repo houses 2 workflows
- KDM  
- DICER
## KDM workflow
This is a modular workflow for detection of mutations in the kinase domain of BCR-ABL fusion gene. The pipeline is implemented in Nextflow. It aligns DNA sequencing reads from 2 replicates of each sample to human hg19 reference genome, detects variants and annotates them with functional information. Coverage over target regions is calculated using bedtools and Coverview. The pipeline integrates the annotated outputs of variant detection tools(VarDict, LoFreq and Mutect2 ) as well as the coverage information, for each sample. It then integrates the output from both replicates into a single spreadsheet. Additional outputs are sorted & indexed bam files, insert size metrics and alignment summary metrics.

## Usage

The following parameters need to be modified in the `params` section of the `nextflow.config` -

- *genome* = Complete path to the human genome fasta file(hg19_all.fasta). Please ensure the FASTA index file(hg19_all.fasta.fai) and BWA index files(hg19_all.fasta.amb, hg19_all.fasta.ann, hg19_all.fasta.bwt, hg19_all.fasta.pac, hg19_all.fasta.sa) are also present in the same genome folder

- *bedfile* = bedfile containing the target regions

- *annovar_humandb* = Complete path to the humandb database folder for ANNOVAR (refer https://annovar.openbioinformatics.org/en/latest/user-guide/startup/ )
---

The list of adaptors required by fastp `TruSeq2-PE.fa` should be placed in the `./assets` folder.

## Executing the pipeline

1. Transfer the `fastq.gz` files to the `sequences/` folder.

2. The samplesheet is `sample_list.csv`. The sample_ids, without the file extension, should be mentioned in samplesheet in the following format - <br>
sample1<br>
sample2<br>
sample3<br>
Please check for empty lines in the samplesheet before running the pipeline.

3. The pipeline can be executed with the following command from the nextflow base directory

```bash
./run_nextflow.sh > script.log
```

---

## Output
The outputs are saved in `Final_output/` folder.

## DICER workflow
This workflow is used to identify mutations in the DICER gene. Adapters and low quality bases in the input fastq.gz files are removed using Trimmomatic. Indel realignment and Base recalibration is performed using GATK to obtain a bam file. Variant calling and annotation is done using Varscan and ANNOVAR respectively.  

**Usage and execution of the workflow is similar to the KDM workflow**  
**Modify the `-entry` value from KDM to DICER in the run_nextflow.sh script and enter the following command**  

```bash
./run_nextflow.sh > script.log
```