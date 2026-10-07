rm(list=ls())
library("phyloseq")
library(ggplot2)
# library(ggpubr)
require('plyr') # ldply
require('reshape2') # melt
project_dir = this.path::here(..=1)
in_dir=file.path(project_dir,"data/phyloseq")
results_dir = file.path(project_dir,"results/alpha_homd_species")
dir.create(results_dir,  showWarnings = F)

#load the phyloseq object and save as "ps"
ps = readRDS(file =file.path(in_dir,"/phyloseq_object.rds"))

set.seed(1)
depth = 100
# measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
calculate_rarefaction_curves <- function(ps, measures, depths) {
  estimate_rarified_richness <- function(ps, measures, depth) {
    if(max(sample_sums(ps)) < depth) return()
    ps <- prune_samples(sample_sums(ps) >= depth, ps)
    rarified_ps <- rarefy_even_depth(ps, depth, verbose = FALSE, replace=F)
    alpha_diversity <- estimate_richness(rarified_ps, measures = measures)
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- reshape2::melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, ps = ps, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}

#run function
rarefaction_curve_data <- calculate_rarefaction_curves(ps, c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), rep(c(100, 500, 1000, 2000, 3500, 5000, 10000, 20000, 50000, 80000, 110000), each = 5))
# summary(rarefaction_curve_data)
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(ps)), by.x = 'Sample', by.y = 'row.names')
#save

measure = "Shannon"
rarefaction_curve_data= data.frame(rarefaction_curve_data_summary_verbose)
measures=c("Shannon","Observed","Simpson")
names(measures) = c("Shannon", "Observed", "Simpson")
for(measure in measures){
  measure_full_name = names(measures)[measures == measure]
  data_to_plot = rarefaction_curve_data %>% dplyr::filter(Measure == measure)
  if (nrow(data_to_plot)){
    plot = ggplot(
      data = data_to_plot,
      mapping = aes(
        x = Depth/1000,
        y = Alpha_diversity_mean,
        ymin = Alpha_diversity_mean - Alpha_diversity_sd,
        ymax = Alpha_diversity_mean + Alpha_diversity_sd,
        colour = patient_id,
        group = sample_id)
      ) + 
      geom_line() + 
      geom_pointrange(fatten = .5, size = 1) + 
      facet_wrap(facets = ~ subgroup) +
      xlim(0,10) + 
      # scale_x_continuous(trans = "log") +
      ggtitle(paste0(measure_full_name," rarefaction curves")) + 
      xlab("Library depth in Thousands") +
      ylab(paste0(measure, " (mean)")) +
      guides(colour=guide_legend(ncol = 2, title = "Patient ID", override.aes = list(size=0.2))) +
      theme(legend.key.size = unit(1, "lines"))
    print(plot)
    ggsave(paste0(results_dir,"/",measure_full_name ,"_rarefaction_curves.jpg"), height = 6, width=7)
  }
}

