#core
rm(list=ls())
library("phyloseq")
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(RColorBrewer)
# BiocManager::install("microbiome")
library(microbiome)
library(gridExtra)

project_dir = this.path::here(..=1)
results_dir = file.path(project_dir,"results/core_microbiome")
dir.create(results_dir, showWarnings = F)
in_dir=file.path(project_dir,"data/phyloseq")

#load the phyloseq object and save as "ps"
pseq <-  readRDS(file =file.path(in_dir,"/phyloseq_object_species.rds"))
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
# pseq.rel <- prune_taxa(taxa_order, pseq.rel.pruned)
pseq.rel.pruned <- prune_taxa(taxa_order, pseq.rel.pruned)


ps_sub <- subset_samples(pseq.rel.pruned, subgroup == "tooth_12m")
title=""
show_y=T
make_core_heatmap_noy <- function(ps_sub, title, show_y = FALSE) {
  detections  <- c(1, 2, 5, 10, 15, 20, 30) / 100
  colours     <- rev(RColorBrewer::brewer.pal(10, "Spectral"))
  g <- microbiome::plot_core(
    ps_sub,
    plot.type   = "heatmap",
    colours     = colours,
    detections  = detections,
    min.prevalence = -1, 
    taxa.order = taxa_order
  ) 
  g <- g +
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
# pseq.rel.pruned@otu_table = pseq.rel.pruned@otu_table+0.000000000001
p_tooth12   <- make_core_heatmap_noy(subset_samples(pseq.rel.pruned, subgroup == "tooth_12m"),
                                     title="Tooth 12 months", show_y = T)
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
                         nrow=1, widths=c(2,1,1,1)), 
             p2_legend, 
             nrow=2,heights=c(10, 1))


ggsave(p, file=paste0(results_dir,"/core_species_detection_v2.jpg"), width=10.47, height=5.17, dpi=600)



# -------------------------------------------------------------------------
# Euler Diagram Visualization
# -------------------------------------------------------------------------

# Install the package if you don't have it yet:
# install.packages("eulerr")

library(eulerr)

# 1. Create a named list of your core taxa
# (These variables: taxa_tooth0, taxa_tooth12, etc. were calculated earlier in your script)
core_input_list <- list(
  "Tooth Baseline"   = taxa_tooth0,
  "Tooth 12m"        = taxa_tooth12,
  "Implant Baseline" = taxa_implant0,
  "Implant 12m"      = taxa_implant12
)

# 2. Fit the Euler diagram model
# We use shape = "ellipse" because 4-way overlaps are mathematically impossible 
# to represent perfectly with circles. Ellipses provide a better fit.
fit <- euler(core_input_list, shape = "ellipse")

# 3. Create the plot object
# quantities = TRUE adds the numbers inside the segments
c <- RColorBrewer::brewer.pal(4,name="Spectral")  # or "Dark2"
p_euler <- plot(fit,
                quantities = list(type = c("counts", "percent")), # Show count and %
                # fills = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), # Distinct colors
                fills = rev(RColorBrewer::brewer.pal(4,name="Spectral")), 
                alpha = 0.6,
                legend = list(side = "right"),
                main = "Core Microbiome Overlaps (Prev > 50%, Detect > 1%)")


# 4. Save the plot
# Note: Euler plots are base/grid graphics, so we use distinct save commands
# Saving as PDF
pdf(file = file.path(results_dir, "core_euler_diagram.pdf"), width = 8, height = 6)
print(p_euler)
dev.off()

# Saving as high-res JPG
jpeg(file = file.path(results_dir, "core_euler_diagram.jpg"), width = 2400, height = 1800, res = 300)
print(p_euler)
dev.off()

