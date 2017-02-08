# Figure 4 - Genetic structuring
#
# NOTE: specific to analysis detailed in github.com/pimbongaerts/bermuda_rad
#
# Author: Pim Bongaerts
# Copyright (C) 2016 Pim Bongaerts
# License: GPL

## Dependencies ==============================================================

suppressMessages(library("ggplot2"))
suppressMessages(library("reshape2"))
suppressMessages(library("gridExtra"))
suppressMessages(library("cowplot"))

## Functions =================================================================

structure_plot <- function(struct_melt, cluster_order, colours, title){
  # STRUCTURE plot function
  plot_struct <- ggplot(struct_melt, aes(x = sample, 
                                         y = value, 
                                         fill = variable)) +
                    geom_bar(stat = "identity", position = "fill", width = 1, 
                             size = 0.25, colour = "black") +
                    scale_fill_manual(name="Cluster",values = colours, 
                                      breaks = cluster_order) +
                    scale_y_continuous(expand = c(0,0)) +
                    facet_grid(~ pop, space = "free_x", 
                               scales = "free_x", switch = "x") +
                    xlab("") +
                    ylab("") +
                    guides(fill = FALSE) +
                    theme(plot.title = element_text(size = 8, 
                                                    vjust = 0, 
                                                    face = "bold"),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 6),
                          axis.ticks = element_blank(),
                          legend.title = element_blank(),
                          panel.margin.x = unit(-0.075, "lines"),
                          panel.margin.y = unit(-2, "lines"),
                          panel.border = element_rect(fill = NA, 
                                                      colour = "black", 
                                                      size = 1, 
                                                      linetype="solid"),
                          strip.background = element_blank(),
                          strip.text.x = element_text(size = 0, 
                                                margin = margin(0, 0, 10, 0)))
  return(plot_struct)
}

pca_plot <- function(scatter_data, colours, shapes, manual_scale, title){
  # Scatter plot function
  plot_scatter <- ggplot(scatter_data, aes(x = pca1,
                                           y = pca2, 
                                           shape = pop, 
                                           color = pop, 
                                           fill = pop)) +
      geom_point(size = 2, color = "black") +
      scale_shape_manual(values = shapes) +
      scale_colour_manual(values = colours) +
      scale_fill_manual(values = colours) +
      scale_x_continuous(limits = c(manual_scale[1], manual_scale[2])) +
      scale_y_continuous(limits = c(manual_scale[3], manual_scale[4])) +
      theme(plot.title = element_text(size = 8, vjust = 0, face="bold"),
            legend.title = element_blank(),
            legend.position = "none",
            axis.text = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
  return(plot_scatter)
}

## Main code =================================================================

### Initial set-up ===========================================================

# Read in Agaricia datasets
afra_struct <- read.csv("afra_data/afra_clumpp_K3.out.csv", header = FALSE)
afra_structneut <- read.csv("afra_data/afra_neut_clumpp_K3.out.csv", 
                            header = FALSE)
afra_pca <- read.csv("afra_data/afra_2f_pca.txt", header = FALSE)

# Read in Stephanocoenia datasets
sint_struct <- read.csv("sint_data/sint_clumpp_K2.out.csv", header = FALSE)
sint_dapcassign <- read.csv("sint_data/sint_bi_2f_posterior.txt", 
                            header = FALSE)
sint_pca <- read.table("sint_data/sint_bi_2f_pca.txt", header = FALSE)

# Add headers
colnames(afra_struct)  <- c("sample", "pop", "clu1", "clu2", "clu3")
colnames(afra_structneut)  <- c("sample", "pop", "clu1", "clu2", "clu3")
colnames(afra_pca) <- c("sample", "pop", "pca1", "pca2", "col", "pch")
colnames(sint_struct)  <- c("sample", "pop", "clu1", "clu2")
colnames(sint_dapcassign)  <- c("sample", "clu1", "clu2", "clu3", "clu4", 
                                "clu5", "clu6", "clu7")
colnames(sint_pca) <- c("sample", "pop", "pca1", "pca2")

# Order of populations
pop_order <- c("GS", "GD", "JS", "JD", "PS", "PD", "WD")

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(-1/24, -1/12, 8/12, 4/12), # Sint - DAPC ASSIGN
                       c(-1/24,  2/12, 8/12, 4/12),  # Sint - STRUCTURE
                       c(16/24,  0/12, 8/24, 6/12),  # Sint - PCA
                       c(-1/24,  5/12, 8/12, 4/12),  # Afra - STRUCTURE
                       c(-1/24,  8/12, 8/12, 4/12),  # Afra - STRUCTURE NEUTR
                       c(16/24,  6/12, 8/24, 6/12))  # Afra - PCA
colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("sint_dapcassign", "sint_struct", "sint_pca",
                             "afra_structneut", "afra_struct", "afra_pca")

### Produce individual plots =================================================

#### Sint - STRUCTURE
title <- "STRUCTURE - Stephanocoenia"
cluster_order <- c("clu1", "clu2")
sint_struct_colours <- c("#9ecae1", "#3182bd")
sint_struct$pop <- factor(sint_struct$pop, levels = pop_order)
sint_struct_melt <- melt(sint_struct, id.vars = c("sample", "pop"), 
                        measure.vars = cluster_order)
plot_sint_struct <- structure_plot(sint_struct_melt, cluster_order, 
                                   sint_struct_colours, title)

#### Sint - DAPC-ASSIGN
title <- "DAPC - Stephanocoenia"
cluster_order <- c("clu1", "clu2", "clu3", "clu4", "clu5", "clu6", "clu7")
sint_dapcassign_colours <- c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", 
                             "#4292c6", "#2171b5", "#084594")
sint_dapcassign$pop <- substr(sint_dapcassign$sample, start = 4, stop = 5)
sint_dapcassign$pop <- factor(sint_dapcassign$pop, levels = pop_order)
sint_dapcassign_melt <- melt(sint_dapcassign, id.vars = c("sample", "pop"), 
                             measure.vars = cluster_order)
plot_sint_dapcassign <- structure_plot(sint_dapcassign_melt, cluster_order, 
                                       sint_dapcassign_colours, title)

#### Sint - PCA
title <- "PCA - Stephanocoenia"
sint_pca_colours <- c("#3182bd", "#9ecae1", "#3182bd", "#9ecae1", 
                      "#3182bd", "#9ecae1", "#3182bd")
sint_shapes <- c(24, 24, 23, 23, 25, 25, 22)
sint_scale <- c(-2.5, 2, -4, 3) # NOTE: some outliers outside of this range
plot_sint_pca <- pca_plot(sint_pca, sint_pca_colours, sint_shapes, 
                          sint_scale, title)

#### Afra - STRUCTURE
title <- "STRUCTURE - Agaricia"
cluster_order <- c("clu1", "clu2", "clu3")
afra_struct_colours <- c("#74c476", "#bae4b3", "#238b45")
afra_struct$pop <- factor(afra_struct$pop, levels = pop_order)
afra_struct_melt <- melt(afra_struct, id.vars = c("sample", "pop"), 
                         measure.vars = cluster_order)
plot_afra_struct <- structure_plot(afra_struct_melt, cluster_order, 
                                   afra_struct_colours, title)

#### Afra - STRUCTURE Neutral
title <- "STRUCTURE NEUTRAL - Agaricia"
cluster_order <- c("clu2", "clu1", "clu3")
afra_structneut_colours <- c("#bae4b3", "#238b45", "#74c476")
afra_structneut$pop <- factor(afra_structneut$pop, levels = pop_order)
afra_structneut_melt <- melt(afra_structneut, id.vars = c("sample", "pop"), 
                              measure.vars = cluster_order)
plot_afra_structneut <- structure_plot(afra_structneut_melt, cluster_order, 
                                       afra_structneut_colours, title)

#### Afra - PCA
title <- "PCA - Agaricia"
afra_pca_colours <- c("#238b45", "#bae4b3", "#238b45", "#bae4b3", 
                      "#238b45", "#bae4b3", "#74c476")
afra_shapes <- c(24, 24, 23, 23, 25, 25, 22)
afra_scale <- c(-2.5, 1, -2.5, 2) # NOTE: some outliers outside of this range
plot_afra_pca <- pca_plot(afra_pca, afra_pca_colours, afra_shapes, 
                          afra_scale, title)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_sint_dapcassign, 
            screen_coords["sint_dapcassign", "x"], 
            screen_coords["sint_dapcassign", "y"], 
            screen_coords["sint_dapcassign", "width"], 
            screen_coords["sint_dapcassign", "height"]) +
  draw_plot(plot_sint_struct, 
            screen_coords["sint_struct", "x"],
            screen_coords["sint_struct", "y"], 
            screen_coords["sint_struct", "width"],
            screen_coords["sint_struct", "height"]) +
  draw_plot(plot_sint_pca, 
            screen_coords["sint_pca", "x"], 
            screen_coords["sint_pca", "y"], 
            screen_coords["sint_pca", "width"], 
            screen_coords["sint_pca", "height"]) +
  draw_plot(plot_afra_structneut, 
            screen_coords["afra_structneut", "x"],
            screen_coords["afra_structneut", "y"], 
            screen_coords["afra_structneut", "width"],
            screen_coords["afra_structneut", "height"]) +
  draw_plot(plot_afra_struct, 
            screen_coords["afra_struct", "x"],
            screen_coords["afra_struct", "y"], 
            screen_coords["afra_struct", "width"],
            screen_coords["afra_struct", "height"]) +
  draw_plot(plot_afra_pca, 
            screen_coords["afra_pca", "x"], 
            screen_coords["afra_pca", "y"], 
            screen_coords["afra_pca", "width"], 
            screen_coords["afra_pca", "height"])

ggsave("figure4.pdf", width = 18, height = 12, units = "cm")
ggsave("figure4.png", width = 18, height = 12, units = "cm")