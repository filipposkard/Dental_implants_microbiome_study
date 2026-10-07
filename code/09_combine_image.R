#!/usr/bin/env Rscript
## Combine two core-microbiome JPEGs into a single vertical figure
rm(list=ls())
options(stringsAsFactors = FALSE)

## dependencies: magick, this.path
if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Please install the 'magick' package: install.packages('magick')")
}

project_dir = this.path::here(..=1)
results_dir = file.path(project_dir, "results/core_microbiome")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

img_species <- file.path(results_dir, "core_species_detection_v2.jpg")
img_genera  <- file.path(results_dir, "core_genera_detection.jpg")
out_file    <- file.path(results_dir, "core_combined_panels.jpg")

if (!file.exists(img_species)) stop(paste("Missing file:", img_species))
if (!file.exists(img_genera))  stop(paste("Missing file:", img_genera))

library(magick)

# read images
i1 <- image_read(img_species)
i2 <- image_read(img_genera)

# resize both to the same width (keep aspect ratio)
info1 <- image_info(i1)
info2 <- image_info(i2)
maxw <- max(info1$width, info2$width)
i1 <- image_resize(i1, paste0(maxw, "x"))
i2 <- image_resize(i2, paste0(maxw, "x"))

# annotate panels A (top) and B (bottom) in white on semi-opaque black box for contrast
label_size <- floor(maxw / 20)
# remove background box and use black labels
i1a <- image_annotate(i1, "A", size = label_size, color = "black", gravity = "northwest", location = "+20+20")
i2a <- image_annotate(i2, "B", size = label_size, color = "black", gravity = "northwest", location = "+20+20")

# stack vertically
combined <- image_append(c(i1a, i2a), stack = TRUE)

# write output (high quality)
image_write(combined, path = out_file, format = "jpg", quality = 95)

cat("Wrote combined image to:", out_file, "\n")
