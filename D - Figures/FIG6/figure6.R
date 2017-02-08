# Figure 6 - Species abundances and community structure
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
suppressMessages(library("plyr"))

## Functions =================================================================
species_plot <- function(data, colours){
  # Species abundances plot
  plot_species <-  ggplot(data, aes(x = location, 
                                y = mean, 
                                fill = speciesdepth)) + 
                          geom_bar(stat = "identity", 
                                   position = "dodge", 
                                   color = "black") +
                          scale_fill_manual(values = colours) +
                          ggtitle("Shallow / deep") + 
                          xlab("") +
                          ylab(expression(coral ~ colonies ~ m^{-2})) +
                          theme(plot.title = element_text(size = 10, 
                                vjust = 0, face="plain"),
                                legend.position = "none",
                                axis.text.x = element_text(size = 8),
                                axis.text.y = element_text(size = ifelse(
                                              data$depth == "shallow", 8, 0)),
                                axis.title.x = element_text(size = 10),
                                axis.title.y = element_text(size = 10),
                                strip.background = element_blank(),
                                strip.text.x = element_text(size = 10),
                                strip.text.y = element_blank(),
                                panel.margin.y = unit(0.2, "lines")) +
                          geom_blank(aes(y = ifelse(
                                                depth == "shallow", -8, 8))) +
                          geom_errorbar(aes(ymin = mean - sterr, 
                                            ymax = mean + sterr), 
                                            width = 0.3) +
                          facet_grid(variable ~ ., scale = "free") + 
                          coord_flip()
  return(plot_species)
}
community_plot <- function(data, colours){
  # Community structure plot
  plot_community <-  ggplot(data, aes(x = 1, y = mean, fill = variable)) +
                      geom_bar(position = "fill", 
                               stat = "identity", 
                               color = "black") + 
                      facet_grid(depth ~ location) + 
                      coord_polar("y", start = 0) +
                      xlab("") +
                      ylab("") +
                      scale_fill_manual(values = colours) +
                      theme(legend.position = "none",
                            axis.title = element_text(size = 10),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            strip.background = element_blank(),
                            strip.text.x = element_text(size = 10, 
                                               margin = margin( 0, 0, 10, 0)),
                            strip.text.y = element_text(size = 10, 
                                               margin = margin( 0, 0, 0, 10)))
  return(plot_community)
}

## Main code =================================================================

### Initial set-up ===========================================================
# Read in datasets
transect <- read.csv("data/transect_data.txt", header = TRUE, sep = "\t")

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c(0/12, 0/12, 4/12, 12/12),  # Species abundances (SH)
                       c(13/48, 0/12, 7/24, 12/12),  # Species abundances (DE)
                       c(6/12, 0/12, 6/12, 12/12))  # Community composition

colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("species_shallow", "species_deep", "community")

mcav <- "M. cavernosa"
afra <- "A. fragilis"
sint <- "S. intersepta"
ofav <- "O. faveolata"

### Produce individual plots =================================================
## Transect matrix modifications
tra_melt <- melt(transect, id = c("location", "depth", "transect"))
levels(tra_melt$location)[levels(tra_melt$location) == "G"] <- "Gurnet"
levels(tra_melt$location)[levels(tra_melt$location) == "J"] <- "John Smith"
levels(tra_melt$location)[levels(tra_melt$location) == "P"] <- "Princess"
tra_melt$depth <- factor(tra_melt$depth, levels=c("shallow", "deep"))
tra_melt_summary <- ddply(tra_melt, c("location", "depth", "variable"), 
                          summarise, mean = mean(value),
                          sterr = sd(value) / sqrt(length(value)))


## Species plots
spp_melt <- subset(tra_melt_summary, 
                   variable %in% c("mcav_sa", "ofav_sa", 
                                   "sint_sa", "afra_sa"))
spp_melt$variable <- factor(spp_melt$variable)
spp_melt$speciesdepth <- paste(spp_melt$variable, spp_melt$depth, sep = "")
spp_melt$mean <- with(spp_melt, ifelse(depth == "shallow", -mean, mean))
levels(spp_melt$variable)[levels(spp_melt$variable) == "mcav_sa"] <- mcav
levels(spp_melt$variable)[levels(spp_melt$variable) == "afra_sa"] <- afra
levels(spp_melt$variable)[levels(spp_melt$variable) == "sint_sa"] <- sint
levels(spp_melt$variable)[levels(spp_melt$variable) == "ofav_sa"] <- ofav
spp_melt$variable <- factor(spp_melt$variable, levels = c(mcav, ofav, 
                                                          afra, sint))
spp_melt_shallow <- subset(spp_melt, depth == "shallow")
spp_colours <- c("#74c476", "#999999", "#cccccc", "#3182bd")
plot_species_shallow <- species_plot(spp_melt_shallow, spp_colours)
spp_melt_deep <- subset(spp_melt, depth == "deep")
plot_species_deep <- species_plot(spp_melt_deep, spp_colours)

#### Community
comm_melt <- subset(tra_melt_summary, 
                         variable %in% c("mcav_prop", "afra_prop", 
                                         "sint_prop", "ofav_prop", 
                                         "other_prop"))
comm_melt$location <- factor(comm_melt$location, 
                             levels = c("Princess", "John Smith", "Gurnet"))
comm_melt$variable <- factor(comm_melt$variable, 
                             levels = c("mcav_prop", "ofav_prop", "sint_prop", 
                                        "afra_prop", "other_prop"))
comm_colours <- c("#999999", "#cccccc", "#3182bd", "#74c476", "#eeeeee")
plot_community <- community_plot(comm_melt, comm_colours)
 
### Composite of individual plots ============================================

ggdraw() +
  draw_plot(plot_species_shallow, 
            screen_coords["species_shallow", "x"], 
            screen_coords["species_shallow", "y"], 
            screen_coords["species_shallow", "width"], 
            screen_coords["species_shallow", "height"]) +
  draw_plot(plot_species_deep, 
            screen_coords["species_deep", "x"], 
            screen_coords["species_deep", "y"], 
            screen_coords["species_deep", "width"], 
            screen_coords["species_deep", "height"]) +
  draw_plot(plot_community, 
            screen_coords["community", "x"], 
            screen_coords["community", "y"], 
            screen_coords["community", "width"], 
            screen_coords["community", "height"])

ggsave("figure6.pdf", width = 18, height = 6, units = "cm")
ggsave("figure6.png", width = 18, height = 6, units = "cm")