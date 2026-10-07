rm(list=ls())
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(tibble)

project_dir = this.path::here(..=1)
in_dir=file.path(project_dir,"data/phyloseq")
results_dir = file.path(project_dir,"results/beta_diversity")
dir.create(results_dir, showWarnings = F)

#load the phyloseq object and save as "ps"
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
#remove unknown taxon
ps <- prune_taxa(taxa_names(ps) != "Unknown", ps)
ps@sam_data$subgroup <- factor(
  ps@sam_data$subgroup,
  levels = c("tooth_12m","tooth_baseline","implant_baseline", "implant_12m")
)
# ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=min(sample_sums(ps)), replace=F)

#plot beta
# Calculate Bray-Curtis distance matrix
bray_dist <- distance(ps.rarefied, method = "bray")

# Perform PCoA ordination using the Bray-Curtis distance
pcoa_bray <- ordinate(ps.rarefied, method = "PCoA", distance = bray_dist)

# Plot the ordination with ggplot2 customization
p1 <- plot_ordination(ps.rarefied, pcoa_bray, color = "subgroup", shape = "site") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = subgroup, color = subgroup), type = "norm", linetype = 2, size = 0.5, level=0.95) +
  geom_text_repel(aes(label = paste0(patient_id)), size = 2, show.legend = FALSE, max.overlaps = Inf,color="black") +
  scale_color_brewer(palette = "Spectral", direction = -1, labels=function(x) gsub("_", " ", gsub("12m", "12 months",x))) +  # or "Dark2"
  theme_minimal() +
  labs(title = "Bray-Curtis distances", color = "Group", shape="Site")

# Print the plot
print(p1)
ggsave(paste0(results_dir,"/Bray_curtis_distance.jpg"), device="jpg",height = 6, width=7, dpi=600)


# Plot the ordination with ggplot2 customization
p3 <- plot_ordination(ps.rarefied, pcoa_bray, color = "patient_id", shape = "subgroup") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = patient_id, color = patient_id), type = "norm", linetype = 2, size = 0.3, level=0.75) +
  geom_text_repel(aes(label = paste0(patient_id)), size = 2, show.legend = FALSE, max.overlaps = Inf,color="black") +
  scale_color_brewer(palette = "BrBG", direction = -1, labels=function(x) gsub("_", " ", gsub("12m", "12 months",x))) +  # or "Dark2"
  theme_minimal() +
  labs(title = "Bray-Curtis distances", color = "Subject", shape="Group")

# Print the plot
print(p3)
ggsave(paste0(results_dir,"/Bray_curtis_distance_subject.jpg"), device="jpg",height = 6, width=7, dpi=600)

#Jaccard
jac_dist <- distance(ps.rarefied, method = "jaccard")
pcoa_jac <- ordinate(ps.rarefied, method = "PCoA", distance = jac_dist)
p2 <- plot_ordination(ps.rarefied, pcoa_jac, color = "subgroup", shape = "site") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = subgroup, color = subgroup), type = "norm", linetype = 2, size = 0.5) +
  geom_text_repel(aes(label = paste0(patient_id)), size = 2, show.legend = FALSE, max.overlaps = Inf,color="black") +
  scale_color_brewer(palette = "Spectral", direction = -1, labels=function(x) gsub("_", " ", gsub("12m", "12 months",x))) +  # or "Dark2"
  theme_minimal() +
  labs(title = "Jaccard distances", color = "Group", shape="Site")  # or use your variable name


# Print the plot
print(p2)
ggsave(paste0(results_dir,"/Jaccard_distances.jpg"), device="jpg",height = 6, width=7, dpi=600)


p4 <- plot_ordination(ps.rarefied, pcoa_jac, color = "patient_id", shape = "subgroup") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = patient_id, color = patient_id), type = "norm", linetype = 2, size = 0.3, level=0.75) +
  geom_text_repel(aes(label = paste0(patient_id)), size = 2, show.legend = FALSE, max.overlaps = Inf,color="black") +
  scale_color_brewer(palette = "BrBG", direction = -1, labels=function(x) gsub("_", " ", gsub("12m", "12 months",x))) +  # or "Dark2"
  theme_minimal() +
  labs(title = "Jaccard distances", color = "Subject", shape="Group")

# Print the plot
print(p4)
ggsave(paste0(results_dir,"/Jaccard_distances_subject.jpg"), device="jpg",height = 6, width=7, dpi=600)

# library(patchwork)
# Combine all plots with patchwork and add labels A, B, C
# combined_plot <- p1 / p2 + plot_annotation(tag_levels = "A")
# Display
# print(combined_plot)

#Combine
p1 <- p1 + guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
p2 <- p2 + guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
# Then combine with ggarrange:
pAB <- ggpubr::ggarrange(
  p1, p2,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)
pAB
# Save combined plot to file
ggsave(plot=pAB, filename = paste0(results_dir,"/combined_beta_diversity_plots_groups.jpg"), width = 5, height = 8, dpi = 600)

#Combine

# Then combine with ggarrange:
pCD <- ggpubr::ggarrange(
  p3, p4,
  labels = c("C", "D"),
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "right"
)
pCD
ggsave(plot=pCD, filename = paste0(results_dir,"/combined_beta_diversity_plots_subject.jpg"), width = 7, height = 8, dpi = 600)
# Then combine with ggarrange:
P <- ggpubr::ggarrange(
  pAB, pCD,
  ncol = 2, nrow = 1,
  common.legend = T)
P
# Save combined plot to file
ggsave(plot=P, filename = paste0(results_dir,"/combined_beta_diversity_plots.jpg"), width = 10, height = 8, dpi = 600)



#

# Calculate Bray-Curtis and Jaccard dissimilarity matrices
bray_dist <- phyloseq::distance(ps.rarefied, method = "bray")
jaccard_dist <- phyloseq::distance(ps.rarefied, method = "jaccard", binary = TRUE)

# Get metadata (make sure row order matches distance matrix)
meta_df <- data.frame(sample_data(ps.rarefied))
# Ensure order matches
meta_df <- meta_df[match(rownames(as.matrix(bray_dist)), rownames(meta_df)), ]

# PERMANOVA for Bray-Curtis
adonis_bray <- vegan::adonis2(as.matrix(bray_dist) ~ subgroup, data = meta_df)
results1 <- as_tibble(as.data.frame(adonis_bray))
results1$Distance <- "Bray-Curtis"
results1$Factor <- "Group"
adonis_bray <- vegan::adonis2(as.matrix(bray_dist) ~ patient_id, data = meta_df)
results2 <- as_tibble(as.data.frame(adonis_bray))
results2$Distance <- "Bray-Curtis"
results2$Factor <- "Subject"
adonis_bray <- vegan::adonis2(as.matrix(bray_dist) ~ subgroup+patient_id, data = meta_df)
results3 <- as_tibble(as.data.frame(adonis_bray))
results3$Distance <- "Bray-Curtis"
results3$Factor <- "Group+Subject"

# PERMANOVA for Jaccard
adonis_jaccard <- vegan::adonis2(as.matrix(jaccard_dist) ~ subgroup, data = meta_df)
results4 <- as_tibble(as.data.frame(adonis_jaccard))
results4$Distance <- "Jaccard"
results4$Factor <- "Group"
adonis_jaccard <- vegan::adonis2(as.matrix(jaccard_dist) ~ patient_id, data = meta_df)
results5 <- as_tibble(as.data.frame(adonis_jaccard))
results5$Distance <- "Jaccard"
results5$Factor <- "Subject"
adonis_jaccard <- vegan::adonis2(as.matrix(jaccard_dist) ~ subgroup+patient_id, data = meta_df)
results6 <- as_tibble(as.data.frame(adonis_jaccard))
results6$Distance <- "Jaccard"
results6$Factor <- "Group+Subject"


# Combine results
permanova_results <- dplyr::bind_rows(results1[1,],results2[1,],results3[1,],results4[1,],results5[1,],results6[1,])

# Write to TSV file
write.table(permanova_results, paste0(results_dir,"/PERMANOVA_results.tsv"), row.names = F, sep="\t")
