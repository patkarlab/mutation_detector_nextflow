# Load required packages
library(VariantAnnotation)
library(Espresso)

# Load other packages
library(doParallel)   # processing samples in parallel
library(magrittr)   # piping and readibility
maf_database="MafDb.gnomADex.r2.1.hs37d5" #If you are using a different MafDB annotation package (see "Installation") change it here accordingly
library(maf_database,character.only = TRUE)

# Specify the genome
genome = "hg19"

dir <- "/scratch/hematopath/hemoseq_v2/mutation_detector_nextflow/espresso"
files=list.files(dir, pattern = "*_unique_nochrM.varscan")
files
sample_paths <- paste0(dir, files)
sample_paths
sample_names <- substr(files,1,10)
sample_names

download.file(url = "https://raw.githubusercontent.com/abelson-lab/Espresso_paper/master/heme_COSMIC/COSMIC_heme_freq10.txt", 
 destfile = "COSMIC_heme_freq10.txt")
download.file(url = "https://raw.githubusercontent.com/abelson-lab/Espresso_paper/master/heme_COSMIC/COSMIC_heme_freq3.txt", 
 destfile = "COSMIC_heme_freq3.txt")

hemeCOSMIC_10 <- load_recurrent_mutations("COSMIC_heme_freq10.txt", genome = "hg19")
hemeCOSMIC_3 <- load_recurrent_mutations("COSMIC_heme_freq3.txt", genome = "hg19")
download.file(url = "https://raw.githubusercontent.com/abelson-lab/Espresso_paper/master/heme_COSMIC/hemeCOSMIC_hotspot_range.bed", 
               destfile = "hemeCOSMIC_hotspot_range.bed")
COSMIC_hotspot_range <- load_bed("hemeCOSMIC_hotspot_range.bed", genome)

# Make cluster and start
cl <- makeCluster(4)
registerDoParallel(cl)

# initialize variant calls
variant_calls <- VRangesList()

# Try for all files
for(i in 1:length(sample_paths)){

  # get sample name and path
  samp_path <- sample_paths[i]
  samp_name <- sample_names[i]

  print(samp_name)

  # get sample as VRanges and annotate with sequence context and MAF
  samp <- load_as_VRanges(samp_name, samp_path, genome) %>%
    sequence_context(., genome, context = 3) %>%
    annotate_MAF(., maf_database, genome)
  
# use sample to generate the error models
    samp_models <- samp %>%
        filter_model_input(., MAF_cutoff = 0.001, VAF_cutoff = 0.05, MAPQ_cutoff = 59, recurrent_mutations = hemeCOSMIC_10) %>%    # preprocessing to clean training set (error models)
        generate_all_models()

# call variants using error models, and aggregate together
    variant_calls[[samp_name]] <- samp %>% 
        intersect_VRanges(., COSMIC_hotspot_range) %>% # only keep positions within the hotspot range we want to call in
      filter_MAPQ(MAPQ_cutoff = 59) %>% 
        call_all_variants(., samp_models)

# cleanup
    rm(samp_name, samp, samp_models)
}

# Unregister Cluster
stopCluster(cl)
rm(cl)
registerDoSEQ()

# Output is in vranges list
variant_calls

# unlist and correct pvalue (bonferroni)
variant_calls_unlisted <- variant_calls %>% unlist() %>% correct_pvalues(., method = "bonferroni")
## if we want to correct sample-by-sample
variant_calls_unlisted <- variant_calls %>% correct_pvalues(., method = "bonferroni") %>% unlist()

# show significant calls 
significant_calls <- variant_calls_unlisted[which(variant_calls_unlisted$corrected_pvalue <= 0.05)]
significant_calls
library('tidyverse')
significant_calls_df <- significant_calls %>% 
  tidyr::as_tibble() %>% 
  dplyr::select(-width, -strand, -totalDepth) %>% 
  dplyr::rename("chr" = seqnames) %>% 
  dplyr::mutate(chr = sub(pattern="chr", replacement="",x=chr))
significant_calls_df %>% readr::write_delim("espresso_5_calls.txt", delim = "\t")

significant_calls_df %>% head() %>% print()

significant_calls %>% 
  VariantAnnotation::asVCF() %>% 
  VariantAnnotation::writeVcf("espresso_5_calls.vcf")

library("cellbaseR")
anno <- annotate_variants(significant_calls, genome='hg19')
#anno %>% head()

df1 <- select(anno, chromosome, start, reference, alternate, id, displayConsequenceType)
write.csv(df1, file = "anno_5.txt")
