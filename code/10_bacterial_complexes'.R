rm(list=ls())
library("lme4")
library("emmeans")
# library(mia)
library("dplyr")
library("phyloseq")
# remotes::install_github("kstagaman/phyloseqCompanion")
library("phyloseqCompanion")
library("ggplot2")
library("ggh4x")
library("lmerTest")


project_dir = this.path::here(..=1)
in_dir=file.path(project_dir,"data/phyloseq")
results_dir = file.path(project_dir,"results/bacterial_complexes_of_interest")
dir.create(results_dir, showWarnings = F)

####color complexes

# Define colors per complex
complex_colors <- c(
  "Red Complex" = "red",
  "Orange Complex" = "orange",
  "Opportunistic" = "mediumpurple",
  " Other" = "lightgrey"
)

bacteria_complex <- c(
  "Porphyromonas gingivalis" = "Red Complex",
  "Tannerella forsythia" = "Red Complex",
  "Treponema denticola" = "Red Complex",
  
  "Fusobacterium nucleatum" = "Orange Complex",
  "Fusobacterium nucleatum_sensu_stricto" =  "Orange Complex", #new nomeclature
  "Fusobacterium periodonticum" = "Orange Complex",
  "Prevotella intermedia" = "Orange Complex",
  "Prevotella nigrescens" = "Orange Complex",
  "Campylobacter rectus" = "Orange Complex",
  "Campylobacter showae" = "Orange Complex",
  "Campylobacter gracilis" = "Orange Complex",
  "Parvimonas micra" = "Orange Complex",
  "Peptostreptococcus micros" = "Orange Complex", #old name, not in eHOMD
  "Eubacterium nodatum" = "Orange Complex",
  "Hornefia nodatum" = "Orange Complex", # eHOMD taxon name instead of "Eubacterium nodatum"
  
  "Streptococcus constellatus" = "Orange Complex",
  
  "Filifactor alocis"="Opportunistic",
  "Parvimonas micra"="Opportunistic",
  "Streptococcus mitis"="Opportunistic",
  "Veillonella dispar"="Opportunistic", #not found in samples
  "Actinomyces israelii"="Opportunistic",
  "Rothia"="Opportunistic",
  "Gemella"="Opportunistic",
  "Neisseria"="Opportunistic",
  "Corynebacterium"="Opportunistic",
  "Staphylococcus epidermidis"="Opportunistic"
)

ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)*100})
phyloseq_df <- psmelt(ps.rel)

sum(phyloseq_df$species %in% names(bacteria_complex))

s <- phyloseq_df$species[phyloseq_df$species %in% names(bacteria_complex) ]
phyloseq_df$Bacterial_complex[phyloseq_df$species %in% names(bacteria_complex) ]  <- bacteria_complex[s]
phyloseq_df$taxon[phyloseq_df$species %in% names(bacteria_complex) ] <- names(bacteria_complex[s])

g <- phyloseq_df$genus[phyloseq_df$genus %in% names(bacteria_complex) ]
phyloseq_df$Bacterial_complex[phyloseq_df$genus %in% names(bacteria_complex) ]  <- bacteria_complex[g]
phyloseq_df$taxon[phyloseq_df$genus %in% names(bacteria_complex) ] <- phyloseq_df$species[phyloseq_df$genus %in% names(bacteria_complex) ]

# Summarize Abundance by Species for each sample
phyloseq_df_summarized <- phyloseq_df %>%
  group_by(sample_id,patient_id,subgroup, site, timepoint,taxon,Bacterial_complex) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")  # Sum abundances of same species
# If not matched (NA), set as "Other"
phyloseq_df_summarized$Bacterial_complex[is.na(phyloseq_df_summarized$Bacterial_complex)] <- " Other"

plot = ggplot(phyloseq_df_summarized, aes(x = patient_id, y = Abundance, fill = Bacterial_complex)) + 
  geom_bar(stat = "identity", color = "black") +  # Removes black outline   , color = NA
  coord_flip() +
  facet_grid2(site ~ timepoint, scales="free_y", space = "free", strip = strip_nested()) +
  scale_fill_manual(values=complex_colors, name = "Bacterial Complex: ") +
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
ggsave(filename = "Bacterial_complexes_barplot.jpg", path = results_dir, width = 11, height = 9)

taxa_plotted <- unique(phyloseq_df_summarized[,c("taxon","Bacterial_complex")])
write.table(taxa_plotted, file = paste0(results_dir, "/Taxon-Complex_pairs_in_this_study.tsv"), sep="\t", row.names = F)

library(tidyr)
df_wide <- phyloseq_df_summarized %>%
  select(-patient_id, -subgroup, -site, -timepoint) %>%  # Drop columns
  pivot_wider(
    names_from = sample_id, 
    values_from = Abundance,
    names_sort = F  # Optional: sorts sample names alphabetically
  )
write.table(df_wide, file = paste0(results_dir, "/Bacterial_complexes_rel_abundance.tsv"), sep="\t", row.names = F)


# # summarise to subgroup level (mean across patients, per complex)
# 1) Sum abundances per patient, subgroup, and complex
df_patient_complex <- phyloseq_df_summarized %>%
  group_by(patient_id, subgroup, Bacterial_complex) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 2) Mean across patients within each subgroup and complex
df_subgroup_mean <- df_patient_complex %>%
  group_by(subgroup, Bacterial_complex) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")


# order subgroups as you used before
df_subgroup_mean$subgroup <- factor(
  df_subgroup_mean$subgroup,
  levels = c("tooth_12m", "tooth_baseline", "implant_baseline", "implant_12m")
)

# barplot of mean relative abundance per subgroup
plot_mean <- ggplot(df_subgroup_mean,
                    aes(x = subgroup, y = Abundance, fill = Bacterial_complex)) +
  geom_bar(stat = "identity", color = "black") +
  xlab("Group") +
  ylab("Mean relative abundance (%)") +
  scale_fill_manual(values = complex_colors, name = "Bacterial complex: ",
                    labels = function(x) gsub(" ", "\n", gsub("-","-\n",x))) +
  guides(fill = guide_legend(
    ncol = 2,
    label.theme = element_text(size = 10),
    title.theme = element_text(size = 11)
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )+
  scale_x_discrete(labels = function(x) gsub("_", "\n", gsub("12m", "12 months",x)))

plot_mean
ggsave(filename = "Bacterial_complexes_mean_barplot.jpg", path = results_dir, width = 5, height = 7)










orange_species <- c("Fusobacterium nucleatum", "Fusobacterium nucleatum_sensu_stricto", "Fusobacterium periodonticum", 
                    "Prevotella intermedia", "Prevotella nigrescens", "Campylobacter rectus", "Campylobacter showae", 
                    "Campylobacter gracilis", "Parvimonas micra", "Peptostreptococcus micros", "Eubacterium nodatum", 
                    "Hornefia nodatum", "Streptococcus constellatus") 

peri_species <- c("Filifactor alocis", "Parvimonas micra", "Streptococcus mitis", "Veillonella dispar", 
                  "Actinomyces israelii", "Rothia", "Gemella", "Neisseria", "Corynebacterium", "Staphylococcus epidermidis") 




ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
ps.rel = ps %>% transform_sample_counts(function(x) {x/sum(x)})
ps.clr <- microbiome::transform(ps, transform = "clr")





# 1. Prepare Data
otu <- as.data.frame(otu_table(ps))
otu <- otu + 0.5 # Simple pseudocount replacement
tax <- as.data.frame(tax_table(ps))

# 2. Identify Rows for the 2 Groups
# Helper to get OTU IDs
get_ids <- function(names) rownames(tax)[tax$species %in% names]
#or
get_ids <- function(names) {
  # 1. Create a regex pattern by joining all names with the OR operator "|"
  # This tells R to look for "Name1" OR "Name2" OR "Name3"
  pattern <- paste(names, collapse = "|")
  
  # 2. Use grepl to find rows where the species column contains the pattern
  # distinct from %in%, this finds substrings (e.g. "Fusobacterium" matches "Fusobacterium nucleatum")
  matches <- grepl(pattern, tax$species)
  
  return(rownames(tax)[matches])
}
ids_orange <- get_ids(orange_species)
ids_peri   <- get_ids(peri_species)


# Define "The Rest" (All taxa NOT in the target group)
# For the Orange analysis, Group A = Orange, Group B = Everything else
ids_rest_orange <- rownames(otu)[!rownames(otu) %in% ids_orange]
ids_rest_peri   <- rownames(otu)[!rownames(otu) %in% ids_peri]


otu_table <- otu 
ids_pos <- ids_orange
ids_neg <- ids_rest_orange


# 3. Function to Calculate Rivera-Pinto Balance --------------------------------
calculate_balance <- function(otu_table, ids_pos, ids_neg) {
  
  # Get sub-tables
  X_pos <- otu_table[ids_pos, , drop=FALSE]
  X_neg <- otu_table[ids_neg, , drop=FALSE]
  
  # Number of parts (k)
  k_pos <- nrow(X_pos)
  k_neg <- nrow(X_neg)
  
  # Calculate Geometric Means (per sample)
  # formula: exp(mean(log(values)))
  # We do this column-wise (per sample)
  g_pos <- apply(X_pos, 2, function(x) exp(mean(log(x))))
  g_neg <- apply(X_neg, 2, function(x) exp(mean(log(x))))
  
  # The Scaling Coefficient (The "Isometric" Part):
  coeff <- sqrt((k_pos * k_neg) / (k_pos + k_neg))
  
  # Calculate Balance
  # B = coeff * log( g_pos / g_neg )
  balance <- coeff * log(g_pos / g_neg)
  
  return(balance)
}


# 4. Calculate Balances --------------------------------------------------------
# Balance 1: Orange Complex vs The Rest
bal_orange <- calculate_balance(otu, ids_orange, ids_rest_orange)

# Balance 2: Peri-implantitis vs The Rest
bal_peri   <- calculate_balance(otu, ids_peri, ids_rest_peri)

# Create Dataframe for Statistics
df_stats <- data.frame(
  Sample = colnames(otu),
  Orange_Balance = bal_orange,
  Peri_Balance = bal_peri
)

# Merge with Metadata
meta <- as.data.frame(sample_data(ps))
df_final <- merge(df_stats, meta, by.x="Sample", by.y="row.names")
df_final$Condition <- paste(df_final$site, df_final$timepoint, sep="_")

# 5. Run Statistics (LMER) -----------------------------------------------------
# Now you are running LMER on the exact Isometric Log-Ratio coordinate
model <- lmer(Orange_Balance ~ Condition + (1|patient_id), data=df_final)

print(anova(model))

# Post-hoc
emm <- emmeans(model, ~ Condition)
pairs(emm, adjust="tukey")


model2 <- lmer(Peri_Balance ~ Condition + (1|patient_id), data=df_final)

print(anova(model2))

# Post-hoc
emm <- emmeans(model2, ~ Condition)
pairs(emm, adjust="tukey")



