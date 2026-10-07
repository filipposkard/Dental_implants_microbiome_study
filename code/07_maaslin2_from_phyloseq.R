rm(list = ls())

library("tidyverse")
# packageVersion("tidyverse")
library("dplyr")
# packageVersion("dplyr")
library("Maaslin2")
# packageVersion("Maaslin2")
# BiocManager::install('Maaslin2')
library("knitr")
library("ggplot2")

project_dir <- this.path::here(.. = 1)
in_dir <- file.path(project_dir, "data/phyloseq")
results_dir <- file.path(project_dir, "results/differential_abundance")
dir.create(results_dir, showWarnings = F)

phyloseq_obj <- readRDS(file = file.path(in_dir, "/phyloseq_object_species.rds"))

minab <- 0 # preferred 0.1
prev <- 0 # preferred 0.3
transform <- "LOG" # "AST" LOG NONE
normalization <- "TSS" # CLR Choices: "TSS", "CLR", "CSS", "NONE", "TMM" TMM and CSS only work on counts, and they also return normalized counts
method <- "LM" # NEGBIN LM
level <- "species"
my_timepoint <- "12m"
my_site <- "tooth"

for (minab in c(0, 0.1)) {
  for (prev in c(0.3)) {
    for (method in c("LM", "NEGBIN")) { # All the non-LM models use an intrinsic log link transformation due to their close connection to GLMs and they are recommended to be run with transform = NONE
      if (method == "LM") {
        transformation_options <- c("LOG", "NONE")
      } else {
        transformation_options <- c("NONE")
      }
      for (transform in transformation_options) {
        if (method == "NEGBIN") {
          normalization_options <- c("TMM", "CSS")
        } else {
          normalization_options <- c("TSS", "TMM", "CSS")
        }
        for (normalization in normalization_options) {
          for (level in c("species", "genus", "phylum", "family")) {
            maaslin_run_dir <- paste0(results_dir, "/", level, "_", transform, "_", normalization, "_minab_", minab, "_prev_", prev, "_method_", method)
            dir.create(maaslin_run_dir, showWarnings = F)
            model <- list()
            ps <- speedyseq::tax_glom(phyloseq_obj, taxrank = level, NArm = F)
            taxa_names(ps) <- ps@tax_table[, level]
            metadata <- data.frame(sample_data(ps))
            counts <- data.frame(otu_table(ps))
            for (my_timepoint in c("baseline", "12m")) {
              metadata_filt <- metadata %>% filter(timepoint == my_timepoint)
              if (nrow(metadata_filt) > 0) {
                model[[my_timepoint]] <- Maaslin2(
                  counts,
                  metadata_filt,
                  min_abundance = minab,
                  min_prevalence = prev,
                  min_variance = 0.0,
                  normalization = normalization,
                  transform = transform,
                  analysis_method = method,
                  output = paste0(maaslin_run_dir, "/", my_timepoint),
                  max_significance = 0.25,
                  random_effects = c("patient_id", "run"),
                  fixed_effects = c("site"),
                  correction = "BH",
                  cores = 50,
                  standardize = FALSE,
                  plot_heatmap = TRUE,
                  plot_scatter = TRUE,
                  heatmap_first_n = 100
                )
              }
            }
            for (my_site in c("tooth", "implant")) {
              metadata_filt <- metadata %>% filter(site == my_site)
              input_data <- counts
              if (nrow(metadata_filt) > 0) {
                model[[my_site]] <- Maaslin2(
                  input_data,
                  metadata_filt,
                  min_abundance = minab,
                  min_prevalence = prev,
                  min_variance = 0.0,
                  normalization = normalization,
                  transform = transform,
                  analysis_method = method,
                  max_significance = 0.25,
                  random_effects = c("patient_id"),
                  fixed_effects = c("timepoint"),
                  correction = "BH",
                  cores = 50,
                  standardize = FALSE,
                  plot_heatmap = TRUE,
                  plot_scatter = TRUE,
                  heatmap_first_n = 100,
                  output = paste0(maaslin_run_dir, "/", my_site)
                )
              }
            }
            saveRDS(model, file = paste0(maaslin_run_dir, "/Maaslin2.rds"))
          }
        }
      }
    }
  }
}
