rm(list=ls())
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggh4x)
library("dplyr")
# install.packages("ggplot2")
# install.packages("tidyverse")

# install.packages("ggrepel")
# install.packages("ggh4x")

# Remove
# remove.packages(c("ggplot2","Seurat","patchwork"))

# Reinstall
# install.packages(c("ggplot2","Seurat","patchwork"))


project_dir = this.path::here(..=1)
in_dir=file.path(project_dir,"data/phyloseq")
results_dir = file.path(project_dir,"results/composition_plots")
dir.create(results_dir, showWarnings = F)
# 

ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})



my_colors <- c("#BCD6D3",
                        "orchid1", "maroon","yellow3",
                        "#FDBF6F", "#E31A1C",       
                        "palegreen2","#FB9A99",       # lt pink
                        "skyblue2",  "gold1","darkorange4",
                        "#FF7F00","#6A3D9A","green4",
                        "khaki2","dodgerblue2", "yellow4", "#CAB2D6","brown","green1", "darkturquoise",
                        "steelblue4","blue1","deeppink1"
)
phyloseq_obj = speedyseq::tax_glom(ps.rel, taxrank="species", NArm = F)
toptaxa = names(sort(taxa_sums(phyloseq_obj),TRUE)[1:15])
tax_table(phyloseq_obj) = cbind(tax_table(phyloseq_obj), top_taxon = ' Other')
tax_table(phyloseq_obj)[toptaxa,"top_taxon"] <- as(tax_table(phyloseq_obj)[toptaxa, "species"], "character")
tax_table(phyloseq_obj)[toptaxa,"top_taxon"]  = gsub("Unknown", "  Unknown", tax_table(phyloseq_obj)[toptaxa,"top_taxon"] )
phyloseq_df <- psmelt(phyloseq_obj)
# Summarize Abundance by Species for each sample
phyloseq_df_summarized <- phyloseq_df %>%
  group_by(patient_id, site, timepoint, top_taxon) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")  # Sum abundances of same species


# 1. Calculate mean abundance per species across all samples
species_means <- phyloseq_df_summarized %>%
  group_by(top_taxon) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

# 2. Reorder factor levels of top_taxon by mean abundance
ordered_levels <- species_means$top_taxon
first_taxon <- " Other"
rest_taxa <- rev(ordered_levels[-1])
# 3. Combine: first taxon stays, rest reversed
new_levels <- c(first_taxon, rest_taxa)

# 4. Set the new factor order
phyloseq_df_summarized$top_taxon <- factor(
  phyloseq_df_summarized$top_taxon,
  levels = new_levels
)



phyloseq_df_summarized$Abundance <- as.numeric(phyloseq_df_summarized$Abundance)
plot = ggplot(phyloseq_df_summarized, aes(x = patient_id, y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(site ~ timepoint, scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 15 Species") +
  xlab("Patient") + ylab("Relative Abundance %") +
  # , axis.text.x = element_text(angle = -90, hjust = 1), axis.text.y = element_text(angle = -90, hjust = 1)) +
  guides(
    fill = guide_legend(
      ncol = 3,
      label.theme = element_text(size = 11),
      title.theme = element_text(size = 11),
      reverse = TRUE
    )
  )
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )

plot
ggsave(filename = "Species_barplot.jpg", path = results_dir, width = 11, height = 9)
plot = ggplot(phyloseq_df_summarized, aes(x = paste(patient_id,site,timepoint), y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(rows="patient_id", scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 15 Species") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11)))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 8),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 12)
  )
plot
ggsave(filename = "Species_barplot_long.jpg", path = results_dir, width = 14, height = 15)



#genus level
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_genus.rds"))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})
phyloseq_obj = speedyseq::tax_glom(ps.rel, taxrank="genus")
toptaxa = names(sort(taxa_sums(phyloseq_obj),TRUE)[1:15])
tax_table(phyloseq_obj) = cbind(tax_table(phyloseq_obj), top_taxon = ' Other')
tax_table(phyloseq_obj)[toptaxa,"top_taxon"] <- as(tax_table(phyloseq_obj)[toptaxa, "genus"], "character")
tax_table(phyloseq_obj)[toptaxa,"top_taxon"]  = gsub("Unknown", "  Unknown", tax_table(phyloseq_obj)[toptaxa,"top_taxon"] )
phyloseq_df <- psmelt(phyloseq_obj)
# Summarize Abundance by genus for each sample
phyloseq_df_summarized <- phyloseq_df %>%
  group_by(patient_id, site, timepoint, top_taxon) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")  # Sum abundances of same species

# 1. Calculate mean abundance per species
species_means <- phyloseq_df_summarized %>%
  group_by(top_taxon) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

# 2. Reverse the ordered species except " Other"
ordered_levels <- rev(species_means$top_taxon)
ordered_levels_no_other <- ordered_levels[ordered_levels != " Other"]

# 3. Put " Other" at the front
final_levels <- c(" Other", ordered_levels_no_other)

# 4. Reorder factor levels in your data frame
phyloseq_df_summarized$top_taxon <- factor(
  phyloseq_df_summarized$top_taxon,
  levels = final_levels
)

plot = ggplot(phyloseq_df_summarized, aes(x = patient_id, y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(site ~ timepoint, scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 15 Genera") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11),
                             reverse = TRUE))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
plot
ggsave(filename = "Genus_barplot.jpg", path = results_dir, width = 11, height = 9)

plot = ggplot(phyloseq_df_summarized, aes(x = paste(patient_id,site,timepoint), y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(rows="patient_id", scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 15 Genera") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11)))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 12)
  )
plot
ggsave(filename = "Genus_barplot_long.jpg", path = results_dir, width = 11, height = 15)

##########
#family level
my_colors <- c("#BCD6D3",
                        "yellow3",
                        "#FDBF6F", "#E31A1C",       
                        "palegreen2","#FB9A99",       # lt pink
                        "skyblue2",  "gold1","darkorange4",
                        "#FF7F00","#6A3D9A","green4",
                        "khaki2","dodgerblue2", "yellow4", "#CAB2D6","brown","green1", "darkturquoise",
                        "steelblue4","blue1","deeppink1"
)
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_family.rds"))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})
phyloseq_obj = speedyseq::tax_glom(ps.rel, taxrank="family")
toptaxa = names(sort(taxa_sums(phyloseq_obj),TRUE)[1:13])
tax_table(phyloseq_obj) = cbind(tax_table(phyloseq_obj), top_taxon = ' Other')
tax_table(phyloseq_obj)[toptaxa,"top_taxon"] <- as(tax_table(phyloseq_obj)[toptaxa, "family"], "character")
tax_table(phyloseq_obj)[toptaxa,"top_taxon"]  = gsub("Unknown", "  Unknown", tax_table(phyloseq_obj)[toptaxa,"top_taxon"] )
phyloseq_df <- psmelt(phyloseq_obj)
# Summarize Abundance by fammily for each sample
phyloseq_df_summarized <- phyloseq_df %>%
  group_by(patient_id, site, timepoint, top_taxon) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")  # Sum abundances of same species


# 1. Calculate mean abundance per species
species_means <- phyloseq_df_summarized %>%
  group_by(top_taxon) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

# 2. Reverse the ordered species except " Other"
ordered_levels <- rev(species_means$top_taxon)
ordered_levels_no_other <- ordered_levels[ordered_levels != " Other"]

# 3. Put " Other" at the front
final_levels <- c(" Other", ordered_levels_no_other)

# 4. Reorder factor levels in your data frame
phyloseq_df_summarized$top_taxon <- factor(
  phyloseq_df_summarized$top_taxon,
  levels = final_levels
)


plot = ggplot(phyloseq_df_summarized, aes(x = patient_id, y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(site ~ timepoint, scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 12 Families") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11)))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
plot
ggsave(filename = "Family_barplot.jpg", path = results_dir, width = 11, height = 9)

plot = ggplot(phyloseq_df_summarized, aes(x = paste(patient_id,site,timepoint), y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(rows="patient_id", scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 12 Families") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11)))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 12)
  )
plot
ggsave(filename = "Family_barplot_long.jpg", path = results_dir, width = 11, height = 15)





##########
#Phylum level
my_colors <- c("#BCD6D3",
                        "yellow3",
                        "#FDBF6F", "#E31A1C",       
                        "palegreen2","#FB9A99",       # lt pink
                        "skyblue2",  "gold1","darkorange4",
                        "#FF7F00","#6A3D9A","green4",
                        "khaki2","dodgerblue2", "yellow4", "#CAB2D6","brown","green1", "darkturquoise",
                        "steelblue4","blue1","deeppink1"
)
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_phylum.rds"))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})
phyloseq_obj = speedyseq::tax_glom(ps.rel, taxrank="phylum")
toptaxa = names(sort(taxa_sums(phyloseq_obj),TRUE)[1:11])
tax_table(phyloseq_obj) = cbind(tax_table(phyloseq_obj), top_taxon = ' Other')
tax_table(phyloseq_obj)[toptaxa,"top_taxon"] <- as(tax_table(phyloseq_obj)[toptaxa, "phylum"], "character")
tax_table(phyloseq_obj)[toptaxa,"top_taxon"]  = gsub("Unknown", "  Unknown", tax_table(phyloseq_obj)[toptaxa,"top_taxon"] )
phyloseq_df <- psmelt(phyloseq_obj)
# Summarize Abundance by fammily for each sample
phyloseq_df_summarized <- phyloseq_df %>%
  group_by(patient_id, site, timepoint, top_taxon) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")  # Sum abundances of same species


# 1. Calculate mean abundance per species
species_means <- phyloseq_df_summarized %>%
  group_by(top_taxon) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

# 2. Reverse the ordered species except " Other"
ordered_levels <- rev(species_means$top_taxon)
ordered_levels_no_other <- ordered_levels[ordered_levels != " Other"]

# 3. Put " Other" at the front
final_levels <- c(" Other", ordered_levels_no_other)

# 4. Reorder factor levels in your data frame
phyloseq_df_summarized$top_taxon <- factor(
  phyloseq_df_summarized$top_taxon,
  levels = final_levels
)


plot = ggplot(phyloseq_df_summarized, aes(x = patient_id, y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(site ~ timepoint, scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 12 Families") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11)))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
plot
ggsave(filename = "Phylum_barplot.jpg", path = results_dir, width = 11, height = 9)

plot = ggplot(phyloseq_df_summarized, aes(x = paste(patient_id,site,timepoint), y = Abundance, fill = top_taxon)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(rows="patient_id", scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=my_colors, name = "Top 12 Families") +
  xlab("Patient") + ylab("Relative Abundance %") +
  guides(fill = guide_legend(ncol = 4, 
                             label.theme = element_text(size=11),
                             title.theme = element_text(size=11)))
plot = plot +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.text.x = element_text(size = 15),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.spacing.y = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 12)
  )
plot
ggsave(filename = "Phylum_barplot_long.jpg", path = results_dir, width = 11, height = 15)

# Combine Genus and Species barplots with titles
library(magick)

# Read the two images
img_genus <- image_read(file.path(results_dir, "Genus_barplot.jpg"))
img_species <- image_read(file.path(results_dir, "Species_barplot.jpg"))

# Get dimensions
info_genus <- image_info(img_genus)
info_species <- image_info(img_species)

# Determine max width and match both images
max_width <- max(info_genus$width, info_species$width)
img_genus <- image_resize(img_genus, paste0(max_width, "x"))
img_species <- image_resize(img_species, paste0(max_width, "x"))

# Get updated dimensions
info_genus <- image_info(img_genus)
info_species <- image_info(img_species)

# Create title bars above each image
title_height <- 160
white_bg <- image_blank(info_genus$width, title_height, color = "white")
title_genus <- image_annotate(white_bg, "Genus-level composition", size = 72, color = "black", gravity = "center", location = "+0+0")
white_bg2 <- image_blank(info_species$width, title_height, color = "white")
title_species <- image_annotate(white_bg2, "Species-level composition", size = 72, color = "black", gravity = "center", location = "+0+0")

# Append titles to the images
img_genus_with_title <- image_append(c(title_genus, img_genus), stack = TRUE)
img_species_with_title <- image_append(c(title_species, img_species), stack = TRUE)

# Add panel labels A and B to each
label_size <- 90
img_genus_with_title <- image_annotate(img_genus_with_title, "A", size = label_size, color = "black", gravity = "northwest", location = "+20+20")
img_species_with_title <- image_annotate(img_species_with_title, "B", size = label_size, color = "black", gravity = "northwest", location = "+20+20")

# Stack vertically
combined <- image_append(c(img_genus_with_title, img_species_with_title), stack = TRUE)

# Save the combined image
image_write(combined, path = file.path(results_dir, "Genus_Species_composition_combined.jpg"), format = "jpg", quality = 95)


