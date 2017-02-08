# Figure 1 - Morphometrics (remainder of figure created in Adobe Illustrator)
#
# NOTE: specific to analysis detailed in github.com/pimbongaerts/bermuda_rad
#
# Author: Pim Bongaerts
# Copyright (C) 2016 Pim Bongaerts
# License: GPL

## Dependencies ==============================================================

suppressMessages(library("ggplot2"))
suppressMessages(library("plyr"))

## Functions =================================================================
morpho_plot <- function(data, colours, shapes){
  # Scatter plot function
  plot_morpho <-  ggplot(morpho_summary, aes(x = cordiam_mean, 
                                             y = cordens_mean, 
                                             fill = speciespop,
                                             shape = speciespop)) +
                  geom_errorbar(aes(ymin = cordens_errmin, 
                                    ymax = cordens_errmax), 
                                    colour = "#AAAAAA", 
                                    width = 0) +
                  geom_errorbarh(aes(xmin = cordiam_errmin, 
                                     xmax = cordiam_errmax), 
                                     colour = "#AAAAAA", 
                                     height = 0) +
                  geom_point(size = 3) +
                  scale_shape_manual(values = morpho_shapes) +
                  scale_fill_manual(values = morpho_colours) +
                  xlab("corallite diameter (mm)") +
                  ylab("corallite density (# / cm2)") +
                  theme(legend.position = "none",
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8),
                        axis.title.x = element_text(size = 8),
                        axis.title.y = element_text(size = 8))
  return(plot_morpho)
}

## Main code =================================================================

### Initial set-up ===========================================================

# Read in datasets
morphometrics <- read.csv("data/morpho_data.csv", header = TRUE)

### Produce plot =================================================

#### Morphometrics
morpho_summary <- ddply(morphometrics,  c("species", "pop"), summarise,
                        N = length(avg_cor_dens),
                        cordens_mean = mean(avg_cor_dens),
                        cordens_se   = sd(avg_cor_dens) / sqrt(N),
                        cordens_errmin = cordens_mean - cordens_se,
                        cordens_errmax = cordens_mean + cordens_se,
                        cordiam_mean = mean(avg_cor_diam),
                        cordiam_se   = sd(avg_cor_diam) / sqrt(N),
                        cordiam_errmin = cordiam_mean - cordiam_se,
                        cordiam_errmax = cordiam_mean + cordiam_se)
morpho_summary <- subset(morpho_summary, pop != "PX")
morpho_shapes <- c(24, 24, 23, 25, 25, 22, 24, 24, 23, 25, 25, 22)
morpho_colours <- c("#238b45", "#bae4b3", "#238b45", 
                    "#238b45", "#bae4b3", "#74c476",
                    "#3182bd", "#9ecae1", "#3182bd", 
                    "#3182bd", "#9ecae1", "#3182bd")
morpho_summary$speciespop <- paste(morpho_summary$species, 
                                   morpho_summary$pop, 
                                   sep = "_")
morpho_plot(morpho_summary, morpho_colours, morpho_shapes)

### Output plots ============================================
ggsave("figure1.pdf", width = 18, height = 12, units = "cm")
ggsave("figure1.png", width = 18, height = 12, units = "cm")