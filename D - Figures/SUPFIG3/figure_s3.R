# Figure S3 - Genetic structuring plots
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
                    facet_grid(~ pop, space = "free_x", scales = "free_x", 
                               switch = "x") +
                    ggtitle(title) +
                    xlab("") +
                    ylab("") +
                    guides(fill = FALSE) +
                    theme(plot.title = element_text(size = 8, 
                                                    vjust = 0, 
                                                    face = "bold"),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 8),
                          axis.ticks = element_blank(),
                          legend.title = element_blank(),
                          panel.margin.x = unit(-0.075, "lines"),
                          panel.margin.y = unit(-2, "lines"),
                          panel.border = element_rect(fill = NA, 
                                                      colour = "black", 
                                                      size = 1, 
                                                      linetype = "solid"),
                          strip.background = element_blank(),
                          strip.text.x = element_text(size = 8, 
                                               margin = margin( 0, 0, 10, 0)))
  return(plot_struct)
}

## Main code =================================================================

### Initial set-up ===========================================================

# Read in Stephanocoenia datasets
sint_struct <- read.csv("sint_data/clumpp_K7.out.csv", header = FALSE)
sint_structshuffle <- read.csv("sint_data/clumpp_shuffle_K7.out.csv", 
                               header = FALSE)

# Add headers
colnames(sint_struct)  <- c("sample", "pop", "clu1", "clu2", "clu3", "clu4", 
                            "clu5", "clu6", "clu7")
colnames(sint_structshuffle)  <- c("sample","pop", "clu1", "clu2", "clu3", 
                                   "clu4", "clu5", "clu6", "clu7")

# Order of populations
pop_order <- c("GS", "GD", "JS", "JD", "PS", "PD", "WD")

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(0/12, 4/12, 12/12, 8/12), # Sint - STRUCT shuffle
                       c(0/12, -2/12, 12/12, 8/12)) # Sint - STRUCT
                       
colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("sint_struct", "sint_structshuffle")

### Produce individual plots =================================================

cluster_order <- c("clu1", "clu2", "clu3", "clu4", "clu5", "clu6", "clu7")
struct_colours <- c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", 
                    "#2171b5", "#084594")

#### Sint - STRUCTURE
title <- "Stephanocoenia (90-100th FST perc)"
sint_struct$pop <- factor(sint_struct$pop, levels = pop_order)
sint_struct_melt <- melt(sint_struct, id.vars = c("sample", "pop"), 
                        measure.vars = cluster_order)
plot_sint_struct <- structure_plot(sint_struct_melt, 
                                   cluster_order, 
                                   struct_colours, 
                                   title)

#### Sint - STRUCTURE
title <- "Stephanocoenia (90-100th FST perc / shuffled)"
sint_structshuffle$pop <- factor(sint_structshuffle$pop, levels = pop_order)
sint_structshuffle_melt <- melt(sint_structshuffle, 
                                id.vars = c("sample", "pop"), 
                                measure.vars = cluster_order)
plot_sint_structshuffle <- structure_plot(sint_structshuffle_melt, 
                                          cluster_order, 
                                          struct_colours, 
                                          title)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_sint_struct, 
            screen_coords["sint_struct", "x"], 
            screen_coords["sint_struct", "y"], 
            screen_coords["sint_struct", "width"], 
            screen_coords["sint_struct", "height"]) +
  draw_plot(plot_sint_structshuffle, 
            screen_coords["sint_structshuffle", "x"], 
            screen_coords["sint_structshuffle", "y"], 
            screen_coords["sint_structshuffle", "width"], 
            screen_coords["sint_structshuffle", "height"]) +
  
ggsave("figure_s3.pdf", width = 18, height = 6, units = "cm")
ggsave("figure_s3.png", width = 18, height = 6, units = "cm")