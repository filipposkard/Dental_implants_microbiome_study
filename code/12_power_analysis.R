## Install once if needed
# install.packages("pwr")

library(pwr)
library(dplyr)
library(phyloseq)

project_dir = this.path::here(..=1)
results_dir = file.path(project_dir,"results/alpha_diversity")
dir.create(results_dir, showWarnings = F)

in_dir=file.path(project_dir,"data/phyloseq")
ps = readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
#remove unknown taxon
ps <- prune_taxa(taxa_names(ps) != "Unknown", ps)
ps@sam_data$subgroup <- factor(
  ps@sam_data$subgroup,
  levels = c("tooth_12m", "tooth_baseline", "implant_baseline", "implant_12m")
)

#rarefy:
set.seed(1)
ps.rarefied <- rarefy_even_depth(ps, rngseed = 1, sample.size = min(sample_sums(ps)), replace = FALSE)
alpha <- estimate_richness(ps.rarefied, split = TRUE,
                           measures = c("Observed", "Shannon", "Simpson"))
alpha$sample_id <- rownames(alpha)
alpha <- merge(alpha,
               data.frame(sample_data(ps.rarefied)),
               by.x = "sample_id", by.y = "row.names")

alpha$subgroup <- factor(alpha$subgroup)  # 4 levels

# Helper: compute Cohen's f from one-way ANOVA
effect_size_f <- function(y, group){
  df <- data.frame(y = y, group = group)
  df <- df[complete.cases(df), ]
  aov_fit <- aov(y ~ group, data = df)
  ss <- summary(aov_fit)[[1]]
  SSb <- ss["group", "Sum Sq"]              # between-groups SS
  SSt <- sum(ss[,"Sum Sq"])                 # total SS
  eta2 <- SSb / SSt                         # eta-squared
  f <- sqrt(eta2 / (1 - eta2))              # Cohen's f
  return(f)
}

# Function to compute power-based N for one alpha metric
power_for_metric <- function(metric_name){
  y <- alpha[[metric_name]]
  f <- effect_size_f(y, alpha$subgroup)
  
  cat("\nMetric:", metric_name, "\n")
  cat("Estimated effect size (Cohen's f):", round(f, 3), "\n")
  
  # If effect is ~0 (no signal), power calc is not meaningful
  if(is.na(f) || f <= 0){
    cat("Effect size is zero or undefined; cannot compute meaningful sample size.\n")
    return(invisible(NULL))
  }
  
  # Power analysis: 4 groups, 80% power, alpha = 0.05
  pwr_res <- pwr.anova.test(k = nlevels(alpha$subgroup),
                            f = f,
                            sig.level = 0.05,
                            power = 0.80)
  
  cat("Required total sample size for 80% power:", ceiling(pwr_res$n * pwr_res$k), "\n")
  cat("Required sample size per group (assuming equal n):", ceiling(pwr_res$n), "\n")
  
  invisible(pwr_res)
}

# Run for your three alpha-diversity metrics
res_observed <- power_for_metric("Observed")
res_shannon  <- power_for_metric("Shannon")
res_simpson  <- power_for_metric("Simpson")

################
#beta diversity measures, power analysis

bray_dist <- phyloseq::distance(ps.rarefied, method = "bray")
jaccard_dist <- phyloseq::distance(ps.rarefied, method = "jaccard", binary = TRUE)
meta_df <- data.frame(sample_data(ps.rarefied))
meta_df <- meta_df[match(rownames(as.matrix(bray_dist)), rownames(meta_df)), ]

# PERMANOVA for Bray-Curtis
adonis_bray <- vegan::adonis2(as.matrix(bray_dist) ~ subgroup, data = meta_df)
# PERMANOVA for Jaccard
adonis_jaccard <- vegan::adonis2(as.matrix(jaccard_dist) ~ subgroup, data = meta_df)


R2_bray <- adonis_bray[1, "R2"]    # Bray-Curtis R²
R2_jaccard <- adonis_jaccard[1, "R2"]  # Jaccard R²

cat("Bray-Curtis: R² =", round(R2_bray, 3), "\n")
cat("Jaccard:    R² =", round(R2_jaccard,3), "\n")

# Convert R² to Cohen's f
f_bray    <- sqrt(R2_bray    / (1 - R2_bray))
f_jaccard <- sqrt(R2_jaccard / (1 - R2_jaccard))

cat("Bray-Curtis:  R² =", round(R2_bray, 3),    " → f =", round(f_bray, 3), "\n")
cat("Jaccard:      R² =", round(R2_jaccard, 3), " → f =", round(f_jaccard, 3), "\n")


# Function to compute required N from f
beta_power_from_f <- function(f, k = 4, power = 0.80, alpha = 0.05){
  if(is.na(f) || f <= 0){
    cat("Effect size f is zero/NA; cannot compute meaningful sample size.\n")
    return(invisible(NULL))
  }
  res <- pwr.anova.test(k = k, f = f, sig.level = alpha, power = power)
  total_n <- ceiling(res$n * k)
  per_group <- ceiling(res$n)
  list(result = res, total_n = total_n, per_group = per_group)
}

# Bray-Curtis
pow_bray <- beta_power_from_f(f_bray, k = nlevels(meta_df$subgroup))
cat("\nBray-Curtis PERMANOVA power:\n")
cat("Required total N for 80% power:", pow_bray$total_n, "\n")
cat("Required N per group:", pow_bray$per_group, "\n")

# Jaccard
pow_jac <- beta_power_from_f(f_jaccard, k = nlevels(meta_df$subgroup))
cat("\nJaccard PERMANOVA power:\n")
cat("Required total N for 80% power:", pow_jac$total_n, "\n")
cat("Required N per group:", pow_jac$per_group, "\n")
