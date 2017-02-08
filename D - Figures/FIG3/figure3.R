# Figure 3 - Outlier identification and frequencies
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

hist_plot <- function(fst_data, title, colour){
  # Histogram function
  plot_hist <- ggplot(fst_data, aes(x = fst)) +
                geom_histogram(fill= colour, colour = "white", bins = 25) +
                ggtitle(title) + 
                ylab("Frequency") +
                coord_flip() +
                scale_x_continuous(limits = c(-0.05, 0.61)) +
                theme(plot.title = element_text(size = 10, 
                                                vjust = 0, 
                                                face = "bold"),
                      axis.title.x = element_text(size = 8),
                      axis.title.y = element_blank(),
                      axis.text.x = element_text(size = 8),
                      axis.text.y = element_text(size = 6),
                      plot.margin = unit(c(2.5, 0, 3, 3), "mm"))
  return(plot_hist)
}

outlier_plot <- function(outlier_data, colours, title){
  # Outlier scatter plot function
  plot_outlier <- ggplot(outlier_data, aes(x = het, 
                                           y = fst, 
                                           color = outlier_type)) +
                    geom_point(size = 0.6) +
                    ggtitle(title) +
                    xlab(expression(Heterozygosity ~ "(" * H[E] *")")) +
                    ylab(expression(F[ST])) +
                    scale_color_manual(values=colours) +
                    scale_x_continuous(limits = c(0.0, 0.55)) +
                    scale_y_continuous(limits = c(-0.05, 0.61)) +
                    theme(plot.title = element_text(size = 10, 
                                                    vjust = 0, 
                                                    face = "bold"),
                          axis.title.x = element_text(size = 8),
                          axis.title.y = element_text(size = 8),
                          axis.text.x = element_text(size = 8),
                          axis.text.y = element_text(size = 8),
                          panel.border = element_rect(colour = "black", 
                                                      fill = NA, 
                                                      size = 0.5, 
                                                      linetype = "solid"),
                          legend.position="none")
  return(plot_outlier)
}

genotype_plot <- function(genotype_data, colours, title){
  # Genotype plot function
  genotype_data$loc <- substr(genotype_data$pop, start = 1, stop = 1)
  genotype_data$freq <- with(genotype_data, 
                             ifelse(group == "shallow",  -freq, freq))
  plot_genotype <- ggplot(genotype_data, aes(x = loc, y = freq, fill = gt)) +  
                    geom_bar(stat="identity") +
                    facet_grid(group ~ chrom_pos, 
                               space = "free_x", 
                               scales = "free") +
                    ggtitle(title) +
                    xlab("") +
                    ylab("Frequency") +
                    scale_fill_manual(values = colours) +
                    theme(plot.title = element_text(size = 10, 
                                                    vjust = 0, 
                                                    face = "bold"),
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
afra_otl_types <- read.table("afra_data/outlier_summary_types.txt", 
                                header = TRUE)
afra_gtfreqs <- read.table("afra_data/afra_b3c_freqs.txt", header = TRUE)

# Read in Stephanocoenia datasets
sint_otl_types <- read.table("sint_data/outlier_summary_types.txt", 
                                header = TRUE)
sint_gtfreqs <- read.table("sint_data/sint_3c_freqs.txt", header = TRUE)

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(17/48, 0/12,  5/24, 25/48),  # Sint - Histogram
                       c(0/12,  0/12,  5/12, 25/48),  # Sint - Outliers
                       c(7/12,  0/12,  5/12, 25/48),  # Sint - Allele freqs
                       c(17/48, 12/24, 5/24, 25/48),  # Afra - Histogram
                       c(0/12,  12/24, 5/12, 25/48),  # Afra - Outliers
                       c(7/12,  12/24, 5/12, 25/48))  # Afra - Allele freqs

colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("sint_hist", "sint_outlier", "sint_genotype",
                             "afra_hist", "afra_outlier", "afra_genotype")

### Produce individual plots =================================================

#### Afra - histogram
title <- ""
afra_hist_colour <- "#74c476"
plot_afra_hist <- hist_plot(afra_otl_types, title, afra_hist_colour)

#### Afra - outliers
title <- ""
afra_outl_colours <- c("#de2d26", "#000000", "#888888", "#CCCCCC")
plot_afra_outlier <- outlier_plot(afra_otl_types, afra_outl_colours, title)

#### Afra - genotype freqs
title <- ""
afra_gt_colours <- c("#de2d26", "#fc9272", "#fee5d9")
plot_afra_genotype <- genotype_plot(afra_gtfreqs, afra_gt_colours, title)

#### Sint - histogram
title <- ""
sint_hist_colour <- "#9ecae1"
plot_sint_hist <- hist_plot(sint_otl_types, title, sint_hist_colour)

#### Sint - outliers
title <- ""
sint_outl_colours <- c("#000000", "#888888", "#CCCCCC")
plot_sint_outlier <- outlier_plot(sint_otl_types, sint_outl_colours, title)

#### Sint - genotype freqs
title <- ""
sint_gt_colours <- c("#000000", "#999999", "#CCCCCC")
plot_sint_genotype <- genotype_plot(sint_gtfreqs, sint_gt_colours, title)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_afra_hist, 
            screen_coords["afra_hist", "x"], 
            screen_coords["afra_hist", "y"],
            screen_coords["afra_hist", "width"], 
            screen_coords["afra_hist", "height"]) +
  draw_plot(plot_afra_outlier, 
            screen_coords["afra_outlier", "x"], 
            screen_coords["afra_outlier", "y"],
            screen_coords["afra_outlier", "width"], 
            screen_coords["afra_outlier", "height"]) +
  draw_plot(plot_afra_genotype, 
            screen_coords["afra_genotype", "x"], 
            screen_coords["afra_genotype", "y"],
            screen_coords["afra_genotype", "width"], 
            screen_coords["afra_genotype", "height"]) +
  draw_plot(plot_sint_hist, 
            screen_coords["sint_hist", "x"], 
            screen_coords["sint_hist", "y"],
            screen_coords["sint_hist", "width"], 
            screen_coords["sint_hist", "height"]) +
  draw_plot(plot_sint_outlier, 
            screen_coords["sint_outlier", "x"], 
            screen_coords["sint_outlier", "y"],
            screen_coords["sint_outlier", "width"], 
            screen_coords["sint_outlier", "height"]) +
  draw_plot(plot_sint_genotype, 
            screen_coords["sint_genotype", "x"], 
            screen_coords["sint_genotype", "y"],
            screen_coords["sint_genotype", "width"], 
            screen_coords["sint_genotype", "height"])

ggsave("figure3.pdf", width = 18, height = 12, units = "cm")
ggsave("figure3.png", width = 18, height = 12, units = "cm")