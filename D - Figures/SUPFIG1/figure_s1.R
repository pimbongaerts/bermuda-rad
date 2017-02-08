# Figure S1 - Outlier frequency plot
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
genotype_plot <- function(genotype_data, colours){
  # Genotyple plot function
  genotype_data$loc <- substr(genotype_data$pop, start = 1, stop = 1)
  genotype_data$freq <- with(genotype_data, 
                             ifelse(group == "shallow", -freq, freq))
  plot_genotype <- ggplot(genotype_data, aes(x = loc, y = freq, fill = gt)) +  
    geom_bar(stat="identity") +
    facet_grid(group ~ chrom_pos, space = "free_x", scales = "free") +
    xlab("") +
    ylab("Frequency") +
    scale_fill_manual(values = colours) +
    theme(plot.title = element_text(size = 10, vjust = 0, face="bold"),
          legend.position="none",
          axis.text.x = element_text(size = 8, angle= 90),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          panel.margin.x = unit(0.2, "lines"),
          panel.margin.y = unit(0, "lines"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 8)
    ) +
    geom_blank(aes(y = ifelse(group == "shallow", -16, 16)))
  return(plot_genotype)
}

## Main code =================================================================

### Initial set-up ===========================================================

# Read in Agaricia datasets
afra_gtfreqs <- read.table("afra_data/afra_b3c_other_freqs.txt", 
                           header = TRUE)

### Produce plot =================================================

afra_genotype_colours <- c("#000000", "#999999", "#CCCCCC")
genotype_plot(afra_gtfreqs, afra_genotype_colours)

ggsave("figure_s1.pdf", width = 18, height = 6, units = "cm")
ggsave("figure_s1.png", width = 18, height = 6, units = "cm")