# Figure S2 - Genetic differentiation plot
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
fst_boxplot <- function(fst_data, colours, title){
  # Genotyple plot function
  plot_fst <- ggplot(fst_data, aes(x = variable, 
                     	           y = value, 
                     	           fill = factor(variable))) +  
    geom_boxplot(lwd = 0.3, outlier.size = 0.1) +
    ggtitle(title) +
    xlab("") +
    ylab(expression(F[ST])) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(-0.1, 0.9)) +
    theme(plot.title = element_text(size = 10, vjust = 0),
          legend.position="none",
          axis.text.x = element_text(size = 6, angle = 90),
          axis.title.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          plot.margin = unit(c(4, 0, 0, 0), "mm"))
  return(plot_fst)
}

## Main code =================================================================

### Initial set-up ===========================================================

# Read in Agaricia datasets
afra_fst <- read.csv("afra_data/afra_5b_fst_compare.csv", header = TRUE)
sint_fst <- read.csv("sint_data/sint_5a_fst_compare.csv", header = TRUE)

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(0/12, 0/12, 6/12, 12/12),  # Afra - Fst Boxplot
                       c(6/12, 0/12, 6/12, 12/12))  # Sint - Fst Boxplot

colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("afra_fst", "sint_fst")

### Produce individual plots =================================================

#### Afra - Fst boxplot
title <- "Agaricia"
afra_fst_colours <- c("lightgray", "lightgray", "lightgray",
                      "#bae4b3", "#bae4b3", "#bae4b3", 
                      "#238b45", "#238b45", "#238b45",
                      "#74c476", "#74c476", "#74c476")
melt_afra_fst <- melt(afra_fst, id = "X")
plot_afra_fst <- fst_boxplot(melt_afra_fst, afra_fst_colours, title)

#### Sint - Fst boxplot
title <- "Stephanocoenia"
sint_genotype_colours <- c("#000000", "#999999", "#CCCCCC")
sint_fst_colours <- c("lightgray", "lightgray", "lightgray",
                      "#9ecae1", "#9ecae1", "#9ecae1",
                      "#3182bd", "#3182bd", "#3182bd",
                      "#3182bd", "#3182bd", "#3182bd")
melt_sint_fst <- melt(sint_fst, id = "X")
plot_sint_fst <- fst_boxplot(melt_sint_fst, sint_fst_colours, title)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_afra_fst, 
            screen_coords["afra_fst", "x"], 
            screen_coords["afra_fst", "y"],
            screen_coords["afra_fst", "width"], 
            screen_coords["afra_fst", "height"]) +
  draw_plot(plot_sint_fst, 
            screen_coords["sint_fst", "x"], 
            screen_coords["sint_fst", "y"],
            screen_coords["sint_fst", "width"], 
            screen_coords["sint_fst", "height"]) +
  
ggsave("figure_s2.pdf", width = 18, height = 8, units = "cm")
ggsave("figure_s2.png", width = 18, height = 8, units = "cm")