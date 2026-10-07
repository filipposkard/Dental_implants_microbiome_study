rm(list=ls())
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(magick)

project_dir = this.path::here(..=1)
results_dir = file.path(project_dir,"results/differential_abundance")
dir.create(results_dir, showWarnings = F)
in_dir=file.path(project_dir,"data/phyloseq")

add_sig_bracket <- function(p, x1, x2, y, label = "*") {
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y) +
    annotate("segment", x = x1, xend = x1, y = y, yend = y * 0.98) +
    annotate("segment", x = x2, xend = x2, y = y, yend = y * 0.98) +
    annotate("text", x = (x1 + x2) / 2, y = y * 1.03, label = label, size = 8, color = "black")
}

# Construct phyloseq object
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_genus.rds"))
ps@sam_data$subgroup <- factor(
  ps@sam_data$subgroup,
  levels = c("tooth_12m", "tooth_baseline", "implant_baseline", "implant_12m")
)
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})

# Make sure the taxonomic name exactly matches your taxonomy table
taxon_abundance <- otu_table(ps.rel)[tax_table(ps.rel)[, "genus"] == "Gemella", ]

# Transpose if needed to have samples as rows and convert to data frame
taxon_abundance_df <- data.frame(sample_id = colnames(taxon_abundance),
                                   abundance = as.numeric(taxon_abundance))

# Add metadata
taxon_abundance_df <- merge(taxon_abundance_df, data.frame(sample_data(ps.rel)), 
                              by.x = "sample_id", by.y = "sample_id")
p_gemella <- ggplot(taxon_abundance_df, aes(x = subgroup, y = abundance, fill = subgroup)) +
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Relative Abundance of Gemella",
       y = "Relative Abundance",
       x = "Group",
       fill = "Group") +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  expand_limits(y = max(taxon_abundance_df$abundance, na.rm = TRUE) * 1.18)

p_gemella <- add_sig_bracket(
  p_gemella,
  x1 = 3,
  x2 = 4,
  y = max(taxon_abundance_df$abundance, na.rm = TRUE) * 1.08,
  label = "*"
)

ggsave(paste0(results_dir,"/Gemella.jpg"), plot = p_gemella, device = "jpg")



# Construct phyloseq object
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
ps@sam_data$subgroup <- factor(
  ps@sam_data$subgroup,
  levels = c("tooth_12m", "tooth_baseline", "implant_baseline", "implant_12m")
)
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})

# Make sure the taxonomic name exactly matches your taxonomy table
taxon_abundance <- otu_table(ps.rel)[tax_table(ps.rel)[, "species"] == "Streptococcus cristatus", ]

# Transpose if needed to have samples as rows and convert to data frame
taxon_abundance_df <- data.frame(sample_id = colnames(taxon_abundance),
                                 abundance = as.numeric(taxon_abundance))

# Add metadata
taxon_abundance_df <- merge(taxon_abundance_df, data.frame(sample_data(ps.rel)), 
                            by.x = "sample_id", by.y = "sample_id")
p_strep <- ggplot(taxon_abundance_df, aes(x = subgroup, y = abundance, fill = subgroup)) +
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Relative Abundance of Streptococcus cristatus",
       y = "Relative Abundance",
       x = "Group",
       fill = "Group") +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  expand_limits(y = max(taxon_abundance_df$abundance, na.rm = TRUE) * 1.18)

p_strep <- add_sig_bracket(
  p_strep,
  x1 = 2,
  x2 = 3,
  y = max(taxon_abundance_df$abundance, na.rm = TRUE) * 1.08,
  label = "*"
)

ggsave(paste0(results_dir,"/Streptococcus cristatus.jpg"), plot = p_strep, device = "jpg")


img_gemella <- image_read(file.path(results_dir, "Gemella.jpg"))
img_strep <- image_read(file.path(results_dir, "Streptococcus cristatus.jpg"))

panel_label_size <- 60
img_gemella <- image_annotate(img_gemella, "A", size = panel_label_size, color = "black", gravity = "northwest", location = "+20+20")
img_strep <- image_annotate(img_strep, "B", size = panel_label_size, color = "black", gravity = "northwest", location = "+20+20")

info_gemella <- image_info(img_gemella)
info_strep <- image_info(img_strep)
max_height <- max(info_gemella$height, info_strep$height)

img_gemella <- image_extent(img_gemella, geometry = paste0(info_gemella$width, "x", max_height), gravity = "northwest", color = "white")
img_strep <- image_extent(img_strep, geometry = paste0(info_strep$width, "x", max_height), gravity = "northwest", color = "white")

combined <- image_append(c(img_gemella, img_strep), stack = FALSE)

image_write(combined, path = file.path(results_dir, "Gemella_Streptococcus_combined.jpg"), format = "jpg", quality = 95)




