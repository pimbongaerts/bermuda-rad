# Figure 5 - Admixture plot
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


heatmap_plot <- function(matrix_melt, colours, title){
  # Heatmap plot function
  plot_heatmap <- ggplot(matrix_melt, aes(x = variable, 
                                          y = chrom_pos, 
                                          fill = value)) + 
                    geom_raster() +
                    coord_flip() +
                    ggtitle(title) +
                    labs(x="Individuals") +
                    scale_fill_manual(values = colours) +
                    guides(fill = FALSE) +
                    theme(plot.title = element_text(size = 8, vjust = 0),
                          axis.text.x = element_text(angle = 90, 
                                                     hjust = 1, 
                                                     size = 6),
                          axis.text.y = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.ticks = element_blank(),
                          panel.border = element_rect(colour = "black", 
                                                      fill = NA, 
                                                      size = 1, 
                                                      linetype = "solid"))
  return(plot_heatmap)
}

structure_plot <- function(struct_melt, cluster_order, colours, title){
  # STRUCTURE plot function
  plot_struct <- ggplot(struct_melt, aes(x = sample, 
                                         y = value, 
                                         fill = variable)) +
                  geom_bar(stat = "identity", 
                           position = "fill", 
                           width = 1.2, 
                           size = 1) +
                  coord_flip() +
                  scale_fill_manual(name="Cluster",values = colours, 
                                    breaks = cluster_order) +
                  scale_y_continuous(expand = c(0,0)) +
                  ggtitle(title) +
                  xlab("") +
                  ylab("") +
                  guides(fill = FALSE) +
                  theme(plot.title = element_text(size = 8, vjust = 0),
                        axis.text.x = element_text(size = 6),
                        axis.text.y =  element_blank(),
                        axis.ticks = element_blank(),
                        legend.title = element_blank(),
                        plot.margin = unit(c(2.6, 6, 3.5, 0), "mm"),
                        panel.border = element_rect(colour = "black", 
                                                    fill = NA, 
                                                    size = 1, 
                                                    linetype = "solid"))
  return(plot_struct)
}


## Main code =================================================================

### Initial set-up ===========================================================

# Read in Agaricia datasets
afra_gts <- read.table("afra_data/afra_matrix_5b.txt", header = TRUE)
afra_assign <- read.csv("afra_data/afra_assign_5b.txt", header = FALSE)

# Add headers
colnames(afra_assign)  <- c("sample", "pop", "clu1", "clu2")

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(19/24, 0/12,  5/24, 12/12),  # Afra - assign
                       c( 0/12, 0/12,  2/12, 12/12),  # Afra - gts outlier
                       c( 7/48, 0/12, 16/24, 12/12))  # Afra - gts deltap

colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("afra_assign", "afra_gts", "afra_gts_deltap")

### Produce individual plots =================================================

#### Afra - Assign
title <- "Overall ancestry"
cluster_order <- c("clu1", "clu2")
afra_assign_colours <- c("#238b45", "#bae4b3")
afra_assign$sample <- factor(afra_assign$sample, 
                             levels = rev(afra_assign$sample))
afra_assign_melt <- melt(afra_assign, id.vars = c("sample"), 
                         measure.vars = cluster_order)
plot_afra_assign <- structure_plot(afra_assign_melt, 
                                   cluster_order, 
                                   afra_assign_colours, 
                                   title)

#### Afra - Outlier genotypes
title <- "Outliers"
afra_gts_outlier_colours = c("#ffffff", "#de2d26", "#fc9272", "#fee5d9")
afra_gts_melt <- melt(afra_gts, id.vars = c("type", "chrom_pos"))
afra_gts_melt$value <- as.character(afra_gts_melt$value)
afra_gts_melt$variable <- factor(afra_gts_melt$variable, 
                                 levels = rev(afra_gts_melt$variable))
afra_gts_melt_outliers <- subset(afra_gts_melt, type == "outlier")
plot_afra_gts <- heatmap_plot(afra_gts_melt_outliers, 
                              afra_gts_outlier_colours, 
                              title)

#### Afra - deltap genotypes
title <- "Delta p"
afra_gts_melt_deltap <- subset(afra_gts_melt, type == "deltap")
afra_gts_deltap_colours = c("#ffffff", "#222222", "#888888", "#dddddd")
plot_afra_gts_deltap <- heatmap_plot(afra_gts_melt_deltap, 
                                     afra_gts_deltap_colours, 
                                     title)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_afra_assign, 
            screen_coords["afra_assign", "x"], 
            screen_coords["afra_assign", "y"],
            screen_coords["afra_assign", "width"], 
            screen_coords["afra_assign", "height"]) +
  draw_plot(plot_afra_gts, 
            screen_coords["afra_gts", "x"], 
            screen_coords["afra_gts", "y"],
            screen_coords["afra_gts", "width"], 
            screen_coords["afra_gts", "height"]) +
  draw_plot(plot_afra_gts_deltap, 
            screen_coords["afra_gts_deltap", "x"], 
            screen_coords["afra_gts_deltap", "y"],
            screen_coords["afra_gts_deltap", "width"], 
            screen_coords["afra_gts_deltap", "height"])
  
ggsave("figure5.pdf", width = 18, height = 8, units = "cm")
ggsave("figure5.png", width = 18, height = 8, units = "cm")