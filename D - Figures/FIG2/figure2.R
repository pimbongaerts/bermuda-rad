# Figure 2 - Genetic distance heatmaps
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

heatmap_plot <- function(matrix, colours, custom_scale, title, pop_order){
  # Heatmap plot function
  plot_heatmap <- ggplot(matrix, aes(x = sample, 
                                     y = variable, 
                                     fill = value)) +
    facet_grid(pop2 ~ pop1, scales = "free", space = "free", switch = "y") +
    geom_tile() +
    ggtitle(title) +
    ylab("Individuals (by population)") +
    scale_fill_gradient2(low = colours[[1]], 
                         mid = colours[[2]], 
                         high = colours[[3]],
                         limit = c(custom_scale[1], custom_scale[3]),
                         midpoint = custom_scale[2],
                         space="Lab") +
    theme(plot.title = element_text(size = 8, vjust = 0, face="italic"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.ticks = element_blank(),
          axis.ticks.margin = unit(0, "lines"),
          panel.margin.x = unit(-0.075, "lines"),
          panel.margin.y = unit(-0.075, "lines"),
          panel.border = element_rect(fill = NA, colour = "black",
                                      size = 1, linetype="solid"),
          strip.background = element_blank(),
          strip.background = element_rect(fill = "grey", colour = "black",
                                          size = 1, linetype="solid"),
          strip.text = element_text(size = 8, margin = margin( 0, 0, 0, 0)),
          strip.switch.pad.grid = unit(-0.23, "cm"))
  return(plot_heatmap)
}


bar_plot <- function(perform_data, colour){
  # Barplot function
  perform_data$samplenr <- substr(perform_data$INDIVIDUAL, 
                                  start = 5, 
                                  stop = 9)
  plot_bar <- ggplot(perform_data, aes(x = samplenr, y = GENO)) +
    geom_bar(stat = "identity", width = 1, size = 0.1, 
             fill = colour, colour = "white") +
    ylab("Number of SNPs") +
    scale_y_reverse() +
    scale_x_discrete(breaks = NULL) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6, margin = margin(0, 0, 0, 0)),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8))
  return(plot_bar)
}

## Main code =================================================================

### Initial set-up ===========================================================

# Read in Agaricia datasets
afra_gd <- read.table("afra_data/afra_wclon_depth_gd_4c.txt", header = TRUE)
afra_perf <- read.table("afra_data/afra_wclon_2f_stats.txt", header = TRUE)

# Read in Stephanocoenia datasets
sint_gd <- read.table("sint_data/sint_depth_withCU_gd_4c.txt", header = TRUE)
sint_perf <- read.table("sint_data/sint_2f_stats.txt", header = TRUE)

# Order of populations
pop_order <- c("GS", "GD", "JS", "JD", "PS", "PD", "WD", "B1")

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(0/12, 0/12, 6/12, 10/48),  # Afra - Missing data
                       c(0/12, 2/12, 6/12, 10/12),  # Afra - GD
                       c(6/12, 0/12, 6/12, 10/48),  # Sint - Missing data
                       c(6/12, 2/12, 6/12, 10/12))  # Sint - GD

colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("afra_perf", "afra_gd",
                             "sint_perf", "sint_gd")

### Produce individual plots =================================================

#### Afra - GD
title <- "Agaricia"
afra_gd_colours <- c("#ff0000", "#ffffff", "#222222")
afra_gd_melt <- melt(afra_gd, id.vars = "sample")
afra_gd_scale <- c(min(afra_gd_melt$value), 0.14, max(afra_gd_melt$value))
afra_gd_melt$value[afra_gd_melt$value == 0.0] <- max(afra_gd_melt$value)
afra_gd_melt$variable = with(afra_gd_melt, 
                             factor(variable, levels = rev(levels(variable))))
afra_gd_melt$pop1 <- substr(afra_gd_melt$sample, start = 4, stop = 5)
afra_gd_melt$pop1[afra_gd_melt$pop1 == "PX"] <- "PD"
afra_gd_melt$pop1 <- factor(afra_gd_melt$pop1, levels = pop_order)
afra_gd_melt$pop2 <- substr(afra_gd_melt$variable, start = 4, stop = 5)
afra_gd_melt$pop2[afra_gd_melt$pop2 == "PX"] <- "PD"
afra_gd_melt$pop2 <- factor(afra_gd_melt$pop2, levels = pop_order)
plot_afra_gd <- heatmap_plot(afra_gd_melt, 
                             afra_gd_colours, 
                             afra_gd_scale, 
                             title, 
                             pop_order)

#### Afra - PERFORMANCE
afra_perf_colour <- "#74c476"
plot_afra_perf <- bar_plot(afra_perf, afra_perf_colour)

#### Sint - GD
title <- "Stephanocoenia"
sint_gd_colours <- c("#ff0000", "#ffffff", "#222222")
sint_gd_scale <- afra_gd_scale
sint_gd_melt <- melt(sint_gd, id.vars = "sample")
sint_gd_melt$value[sint_gd_melt$value == 0.0] <- max(afra_gd_melt$value)
sint_gd_melt$variable = with(sint_gd_melt, 
                             factor(variable, levels = rev(levels(variable))))
sint_gd_melt$pop1 <- substr(sint_gd_melt$sample, start = 4, stop = 5)
sint_gd_melt$pop1[sint_gd_melt$pop1 == "PX"] <- "PD"
sint_gd_melt$pop1 <- factor(sint_gd_melt$pop1, levels = pop_order)
sint_gd_melt$pop2 <- substr(sint_gd_melt$variable, start = 4, stop = 5)
sint_gd_melt$pop2[sint_gd_melt$pop2 == "PX"] <- "PD"
sint_gd_melt$pop2 <- factor(sint_gd_melt$pop2, levels = pop_order)
plot_sint_gd <- heatmap_plot(sint_gd_melt, 
                             sint_gd_colours, 
                             sint_gd_scale, 
                             title, 
                             pop_order)

#### Sint - PERFORMANCE
sint_perf_colour <- "#9ecae1"
plot_sint_perf <- bar_plot(sint_perf, sint_perf_colour)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_afra_gd, 
            screen_coords["afra_gd", "x"], 
            screen_coords["afra_gd", "y"],
            screen_coords["afra_gd", "width"], 
            screen_coords["afra_gd", "height"]) +
  draw_plot(plot_afra_perf, 
            screen_coords["afra_perf", "x"], 
            screen_coords["afra_perf", "y"],
            screen_coords["afra_perf", "width"], 
            screen_coords["afra_perf", "height"]) +
  
  draw_plot(plot_sint_gd, 
            screen_coords["sint_gd", "x"], 
            screen_coords["sint_gd", "y"],
            screen_coords["sint_gd", "width"], 
            screen_coords["sint_gd", "height"]) +
  draw_plot(plot_sint_perf, 
            screen_coords["sint_perf", "x"], 
            screen_coords["sint_perf", "y"],
            screen_coords["sint_perf", "width"], 
            screen_coords["sint_perf", "height"])

ggsave("figure2.pdf", width = 18, height = 12, units = "cm")
ggsave("figure2.png", width = 18, height = 12, units = "cm")