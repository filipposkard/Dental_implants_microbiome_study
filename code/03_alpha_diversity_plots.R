rm(list=ls())
library("phyloseq")
library(ggplot2)
library("ggpubr")
library("ggrepel")
library("knitr")
library("ggsignif")

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
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=min(sample_sums(ps)), replace=F)
alpha <- estimate_richness(ps.rarefied, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
alpha$sample_id <- rownames(alpha)
alpha = merge(alpha, data.frame(sample_data(ps.rarefied)), by.x = "sample_id", by.y = "row.names")

measures=c("Observed", "Shannon", "Simpson")
names(measures) = c("Species Richness", "Shannon Index", "Simpson Index")
x_feature="subgroup"
measure="Shannon"

#Run kruskal walis test
# Create a data frame to store test results
kruskal_results <- data.frame(
  measure = character(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
for (m in measures) {
  test <- kruskal.test(alpha[[m]] ~ alpha$subgroup)
  kruskal_results <- rbind(
    kruskal_results,
    data.frame(
      measure = names(measures[measures==m]),
      statistic = test$statistic,
      p_value = test$p.value
    )
  )
}
kruskal_results
#  save as TSV or TSV file
write.table(kruskal_results, file = file.path(results_dir, "kruskal_wallis_alpha_diversity_results.tsv"), sep="\t",row.names = FALSE)

plot_list<- list()
for(x_feature in c("subgroup")){
  for(measure in measures){
    measure_full_name = names(measures)[measures == measure]
    max_y = max(alpha[,measure])
    min_y = min(alpha[,measure])
    
    #p values
    # test = pairwise.wilcox.test(alpha[,measure], alpha$subgroup, p.adjust.method="bonferroni", paired=T)
    stat_pvalue <- alpha %>% 
      rstatix::pairwise_wilcox_test(formula = as.formula(paste(measure, "~ subgroup")), 
                                                           p.adjust.method = "bonferroni", paired=T) 
    write.table(stat_pvalue, file=paste0(results_dir,"/",measure_full_name,"_wilcox_paired_test.tsv"), sep="\t", row.names = F)
    stat_pvalue <- stat_pvalue %>%
      dplyr::filter(p<=0.05) %>%
      mutate(displayed_p = "") %>%
      mutate(y.position = "")
    if(nrow(stat_pvalue)){
      stat_pvalue <- stat_pvalue %>% 
        rstatix::add_y_position() %>% 
        # mutate(y.position = seq(max_y+(max_y-min_y)/15, max_y+(max_y-min_y)/15+(max_y-min_y)/6,length.out = length(p))) %>%
        mutate(y.position = seq(max_y+(max_y-min_y)/15, max_y+(max_y-min_y)/15+(max_y-min_y)/36*length(p),length.out = length(p))) %>%
        mutate(displayed_p= paste0("p=",round(p,2),", p.adj=",p.adj.signif))
    }
    

    p=ggplot(alpha, aes(x = .data[[x_feature]], y = .data[[measure]], fill =.data[[x_feature]] , group=.data[[x_feature]])) +
      geom_boxplot(alpha=0.4, outliers = F)+
      geom_line(aes(group = patient_id), alpha = 0.5, color = "gray") +  # lines grouped by patient
      # geom_point(aes(color=.data[[x_feature]]))+
      geom_point(aes(color=patient_id), alpha=0.8)+
      # geom_jitter(aes(color = patient_id), width = 0.15, alpha = 0.8, size = 1.5) +  # Jitter points horizontally with width
      geom_text_repel(
        aes(label = patient_id), 
        size = 2, 
        # nudge_y = 0.05,          # small vertical nudge (adjust as needed)
        segment.size = 0.2,      # thinner connector lines
        force = 0.1,             # reduce repulsion force to keep labels nearer
        max.overlaps = Inf
      )+
      stat_summary(color="red", aes(x=.data[[x_feature]], y= .data[[measure]], group = patient_id), fun.y=mean, geom="line", group="patient_id")+
      ggpubr::stat_pvalue_manual(stat_pvalue, label = "displayed_p", 
                         size = 3,                     # smaller font
                         bracket.size = 0.3,           # thinner bracket
                         step.increase = 0.01,         # smaller step between brackets (default is 0.1)
                         tip.length = 0.01)+       # shorter bracket arms (default is 0.03))+
      scale_fill_brewer(palette = "Spectral", direction = -1) +  # or "Dark2"
      scale_color_brewer(palette = "BrBG", direction = -1) +  # or "Dark2"
      theme(legend.position="none", axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=12), plot.title = element_text(hjust = 0.5))+
      labs(title = paste0("Alpha diversity: ",measure_full_name), y = measure_full_name, x = x_feature)
    p = p +
      theme(
        strip.text = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 11),
        panel.spacing.y = unit(0.3, "cm"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 11)
      ) +
      scale_x_discrete(labels = function(x) gsub("_", "\n", gsub("12m", "12 months",x)))
    p
    ggsave(paste0(results_dir,"/",measure_full_name ,"_per_",x_feature,".jpg"), height = 6, width=7, dpi=600)
    plot_list[[measure]] <- p
  }
}
library(patchwork)
# Combine all plots with patchwork and add labels A, B, C
combined_plot <- (plot_list[["Observed"]] + labs(title = NULL, x=NULL)) +
  (plot_list[["Shannon"]] + labs(title = NULL, x=NULL)) +
  (plot_list[["Simpson"]] + labs(title = NULL, x=NULL)) +
  plot_annotation(tag_levels = "A")
# Display
# print(combined_plot)
# Save combined plot to file
ggsave(filename = paste0(results_dir,"/combined_alpha_diversity_plots.jpg"), plot = combined_plot, width = 12, height = 6, dpi = 600)


