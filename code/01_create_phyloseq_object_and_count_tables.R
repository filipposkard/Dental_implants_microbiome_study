rm(list=ls())
library("phyloseq")
library(data.table)

project_dir = this.path::here(..=1)
out_dir = file.path(project_dir,"data/phyloseq")
dir.create(out_dir, showWarnings = F)

#metadata
metadata = read.csv(file.path(project_dir,"metadata/dental_implants_samples_metadata.tsv"), sep="\t", stringsAsFactors = FALSE)
metadata$timepoint = relevel(as.factor(metadata$timepoint), ref = "baseline")
metadata$site=relevel(as.factor(metadata$site), ref = "tooth")
metadata$run=factor(metadata$run)
metadata$site_timepoint=paste0(metadata$site, metadata$timepoint)
rownames(metadata)= metadata$sample_id
metadata$subgroup <- factor(
  metadata$subgroup,
  levels = c("tooth_baseline", "tooth_12m", "implant_baseline", "implant_12m")
)
metadata$patient_id = factor(
  metadata$patient_id,
  levels = c("P1", "P2", "P3", "P4","P5","P6","P7","P8","P9","P10"))

# taxa counts
# Read the TSV file
level <- "species"
level <- "genus"
for(level in c("species","genus","family","phylum")){
  df <- read.table(paste0("/mnt/raid1b/philip/dental_implants/data/epi2me/kraken2_ehomd/",level,"/abundance_table_",level,".tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  tax_split <- strsplit(df$tax, ";")
  taxonomy <- do.call(rbind, lapply(tax_split, function(x) {
    # length(x) <- 8 # Ensure exactly 8 columns (pad with NA if short)
    return(x)
  }))
  colnames(taxonomy) <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")[1:ncol(taxonomy)]
  rownames(taxonomy) <- taxonomy[,level]
  counts_df <-  df[, !(names(df) %in% c("tax", "total"))]

    #change sample names
  # Create a named vector for mapping old sample IDs to new sample IDs
  mapping <- setNames(metadata$sample_id, metadata$old_sample_id)
  colnames(counts_df) <- mapping[colnames(counts_df)]
  #reorder samples
  counts_df <- counts_df[,metadata$sample_id]
  #add row names to counts_df
  rownames(counts_df) = taxonomy[,level]
  
  # Construct phyloseq object
  #
  otu <- otu_table(counts_df, taxa_are_rows = TRUE)
  # 
  sampledata <- sample_data(metadata)
  #
  taxonomy <- tax_table(as.matrix(taxonomy))
  ps <- phyloseq(otu, sampledata, taxonomy)
  saveRDS(ps,file =paste0(out_dir,"/phyloseq_object_",level,".rds"))

  
  #also save counts table
  # Save raw counts as TSV
  dir.create(paste0(project_dir, "/results/abundance_tables"), showWarnings=F)
  raw_counts_outfile <- paste0(project_dir, "/results/abundance_tables/raw_counts_", level, ".tsv")
  # Add taxon column as first column
  raw_counts_df <- cbind(taxon = rownames(counts_df), counts_df)
  write.table(raw_counts_df, raw_counts_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Calculate relative abundances (each sample's counts divided by sample total)
  rel_abund <- sweep(counts_df, 2, colSums(counts_df), FUN = "/")
  rel_abund <-  round(rel_abund*100,2)
  rel_abund_outfile <- paste0(project_dir, "/results/abundance_tables/relative_percent_abundances_", level, ".tsv")
  # Add taxon column as first column
  rel_abund_df <- cbind(taxon = rownames(rel_abund), rel_abund)
  write.table(rel_abund_df, rel_abund_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
}
