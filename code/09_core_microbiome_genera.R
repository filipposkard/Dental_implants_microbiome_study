#core
rm(list=ls())
library("phyloseq")
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(microbiome)

project_dir = this.path::here(..=1)
results_dir = file.path(project_dir,"results/core_microbiome")
dir.create(results_dir, showWarnings = F)
in_dir=file.path(project_dir,"data/phyloseq")

#load the phyloseq object and save as "ps"
pseq <-  readRDS(file =file.path(in_dir,"/phyloseq_object_genus.rds"))
# ps <- prune_taxa(taxa_names(ps) != "Unknown", ps)
pseq.rel <-  microbiome::transform(pseq, "compositional")


#########
ps.rel.t12 <- subset_samples(pseq.rel, subgroup == "tooth_12m")
ps.rel.t0 <- subset_samples(pseq.rel, subgroup == "tooth_baseline")
ps.rel.i0 <- subset_samples(pseq.rel, subgroup == "implant_baseline")
ps.rel.i12 <- subset_samples(pseq.rel, subgroup == "implant_12m")

get_core_subgroup <- function(ps, detection=0.01, prevalence=0.5) {
  core_members(ps, detection = detection, prevalence = prevalence, include.lowest = T)
}

taxa_tooth12   <- get_core_subgroup(ps.rel.t12)
taxa_tooth0    <- get_core_subgroup(ps.rel.t0)
taxa_implant0  <- get_core_subgroup(ps.rel.i0)
taxa_implant12 <- get_core_subgroup(ps.rel.i12)

# union of all taxa that meet the criteria in at least one subgroup
taxa_union <- union(
  union(taxa_tooth12, taxa_tooth0),
  union(taxa_implant0, taxa_implant12)
)

taxa_union
pseq.rel.pruned <- prune_taxa(taxa_union, pseq.rel)

#######
# ensure same taxa order across all subsets
taxa_order <- names(sort(prevalence(pseq.rel.pruned)))
pseq.rel <- prune_taxa(taxa_order, pseq.rel.pruned)

make_core_heatmap_noy <- function(ps_sub, title, show_y = FALSE) {
  detections  <- c(1, 2, 5, 10, 15, 20, 30) / 100
  colours     <- rev(RColorBrewer::brewer.pal(10, "Spectral"))
  g <- plot_core(
    ps_sub,
    plot.type   = "heatmap",
    colours     = colours,
    detections  = detections,
    min.prevalence = 0, 
    taxa.order = taxa_order
  ) +
    labs(
      title = title,
      x = "Detection threshold\n(Relative abundance (%))",
      y = if (show_y) "Taxa" else NULL
    ) +
    scale_x_discrete(labels = function(x) paste0(100 * as.numeric(x),"%")) +
    scale_fill_gradientn("Prevalence", breaks = seq(from = 0, 
                                                    to = 1, by = 0.1), labels = scales::percent, 
                         colours = colours, limits = c(0, 1))
  if (!show_y) {
    g <- g +
      theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.key.width  = unit(1.8, "cm"),   # wider color bar
        legend.key.height = unit(0.6, "cm"),   # taller bar
        legend.title      = element_text(size = 12),
        legend.text       = element_text(size = 9),
        legend.spacing.x  = unit(0.4, "cm")    # space between bar and labels
        # more space between labels and bar
      )
  }
  g
}

p_tooth12   <- make_core_heatmap_noy(subset_samples(pseq.rel.pruned, subgroup == "tooth_12m"),
                                     "Tooth 12 months", show_y = T)
p_tooth0    <- make_core_heatmap_noy(subset_samples(pseq.rel.pruned, subgroup == "tooth_baseline"),
                                     "Tooth baseline", show_y = F)
p_implant0  <- make_core_heatmap_noy(subset_samples(pseq.rel.pruned, subgroup == "implant_baseline"),
                                     "Implant baseline", show_y = F)
p_implant12 <- make_core_heatmap_noy(subset_samples(pseq.rel.pruned, subgroup == "implant_12m"),
                                     "Implant 12 months", show_y = F)

##### fix legend
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(p_tooth0)   

#arranging the legend and plots in a grid:
p <- grid.arrange(arrangeGrob(p_tooth12 + theme(legend.position="none"), 
                              p_tooth0 + theme(legend.position="none"), 
                              p_implant0 + theme(legend.position="none"), 
                              p_implant12 + theme(legend.position="none"), 
                              nrow=1, widths=c(1.65,1,1,1)), 
                  p2_legend, 
                  nrow=2,heights=c(10, 1))


ggsave(p, file=paste0(results_dir,"/core_genera_detection.jpg"), width=10.47, height=5.17, dpi=600)
