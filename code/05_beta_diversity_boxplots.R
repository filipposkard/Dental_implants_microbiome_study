rm(list=ls())
library("phyloseq")
library(ggplot2)
library(ggpubr)
# library(knitr)
# library("MicrobiomeStat")
# library(phangorn)
# library(lmerTest)
library(data.table)
library(dplyr)
library(tidyr)
library(RColorBrewer)


project_dir = this.path::here(..=1)
results_dir = file.path(project_dir,"results/beta_diversity")
dir.create(results_dir, showWarnings = F)
in_dir=file.path(project_dir,"data/phyloseq")

#load the phyloseq object and save as "ps"
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
ps <- prune_taxa(taxa_names(ps) != "Unknown", ps)


ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=min(sample_sums(ps)), replace=F)


# Helper function to get sample ID by site and timepoint
get_sample_id <- function(site_val, timepoint_val) {
  samp <- patient_samples %>% 
    filter(site == site_val & timepoint == timepoint_val) %>%
    pull(SampleID)
  if(length(samp) == 0) return(NA) else return(samp[1])
}

measures = c("bray", "jaccard")
names(measures) = c("Bray-Curtis", "Jaccard")
measure="bray"
for(measure in measures){
  measure_full_name = names(measures)[measures == measure]
  # Calculate beta diversity distance matrix (samples x samples)
  dist_bc <- phyloseq::distance(ps.rarefied, method = measure)
  dist_mat <- as.matrix(dist_bc)
  
  # Extract sample metadata for convenience
  sample_meta <- data.frame(sample_data(ps.rel))
  sample_meta$SampleID <- rownames(sample_meta)
  
  # Initialize a list to store results
  results_list <- list()
  
  # Get unique patient IDs
  patients <- unique(sample_meta$patient_id)
  
  # Loop over patients
  for (patient in patients) {
    # Subset metadata to the current patient
    patient_samples <- sample_meta %>% filter(patient_id == patient)
    
    
    
    # Get sample IDs for each condition
    T0 <- get_sample_id("tooth", "baseline")       # tooth, baseline
    T12 <- get_sample_id("tooth", "12m")             # tooth, 12 months
    I0 <- get_sample_id("implant", "baseline")      # implant, baseline
    I12 <- get_sample_id("implant", "12m")           # implant, 12 months
    
    # Helper function to get distance safely (returns NA if samples missing)
    get_dist <- function(s1, s2) {
      if(is.na(s1) | is.na(s2)) return(NA)
      return(dist_mat[s1, s2])
    }
    
    # Calculate requested distances
    T0_T12 <- get_dist(T0, T12)
    T0_I0 <- get_dist(T0, I0)
    I0_I12 <- get_dist(I0, I12)
    T12_I12 <- get_dist(T12, I12)
    
    # Store in list
    results_list[[patient]] <- data.frame(
      patient_id = patient,
      T0_T12 = T0_T12,
      T0_I0 = T0_I0,
      I0_I12 = I0_I12,
      T12_I12 = T12_I12
    )
  }
  
  # Combine into one dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # View results
  print(results_df)
  
  
  # 
  results_long <- results_df %>%
    pivot_longer(cols = c("T0_T12", "T0_I0", "I0_I12", "T12_I12"),
                 names_to = "Comparison",
                 values_to = "Distance")
  
  results_long$Comparison <- factor(
    results_long$Comparison,
    levels = c("T0_I0", "T0_T12","I0_I12", "T12_I12")
  )
  
  # results_long$patient_id = as.character(results_long$patient_id)
  
  results_long$patient_id = factor(
    results_long$patient_id,
    levels = c("P1", "P2", "P3", "P4","P5","P6","P7","P8","P9","P10"))
  
  ggplot(results_long, aes(x = Comparison, y = Distance)) +
    geom_boxplot(fill = "lightgray", outliers = FALSE) +
    geom_line(aes(group = patient_id, color=patient_id), alpha = 0.5) +  # lines grouped by patient
    geom_point(aes(color=patient_id), size = 1, alpha=0.7) +   # points colored by patient
    # stat_summary(color="red", aes(x=Comparison, y= Distance, group = patient_id), fun.y=mean, geom="line", group="patient_id")+
    stat_compare_means(
      paired = TRUE,
      method = "wilcox.test",
      comparisons = combn(levels(results_long$Comparison),2, simplify = FALSE),
      size = 3,                     # smaller font
      bracket.size = 0.5,           # thinner bracket
      step.increase = 0.05,         # smaller step between brackets (default is 0.1)
      tip.length = 0.01             # shorter bracket arms (default is 0.03)
    )+
    theme_minimal() +
    labs(title = paste0(measure_full_name," Distances"), x = "Comparison", y = paste(measure_full_name,"Distance"))+
    # scale_color_manual(values = colorspace::qualitative_hcl(10)) # terrain_hcl(12))
    scale_color_brewer(palette = "BrBG", direction = -1)  # or "Dark2"
    
  ggsave(paste0(results_dir,"/",measure_full_name," boxplot.jpg"), device="jpg")
}
# 

