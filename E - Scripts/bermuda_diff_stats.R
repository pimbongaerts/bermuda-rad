# Differentiation and diversity stats
#
# NOTE: specific to analysis detailed in github.com/pimbongaerts/bermuda_rad
#
# Authors: Cynthia Riginos, Pim Bongaerts
# Copyright (C) 2015 Cynthia Riginos & Pim Bongaerts
# License: GPL

## Dependencies ==============================================================

suppressMessages(library("adegenet"))
suppressMessages(library("vcfR"))
suppressMessages(library("hierfstat"))
suppressMessages(library("reshape"))

## Constants =================================================================

COL_OTHER <- "lightgray"

## Functions =================================================================

per_loc_fsts <- function(genind_pop1, genind_pop2){
  # Create matrix of per-locus F-statistics for combination of populations
  return(as.data.frame(wc(repool(genind_pop1, genind_pop2))$per.loc))
}

per_loc_pairwise_fsts <- function(coral.genind.bypop){
  # Create matrix with different pairwise comparisons

  # Between habitats (within location)
  Fstats.Habt.G <- per_loc_fsts(coral.genind.bypop$GD, coral.genind.bypop$GS)
  Fstats.Habt.J <- per_loc_fsts(coral.genind.bypop$JD, coral.genind.bypop$JS)
  Fstats.Habt.P <- per_loc_fsts(coral.genind.bypop$PD, coral.genind.bypop$PS)
  
  #Between locations (within shallow)
  Fstats.Shal.GJ <- per_loc_fsts(coral.genind.bypop$GS, coral.genind.bypop$JS)
  Fstats.Shal.JP <- per_loc_fsts(coral.genind.bypop$JS, coral.genind.bypop$PS)
  Fstats.Shal.GP <- per_loc_fsts(coral.genind.bypop$GS, coral.genind.bypop$PS)
  
  #Between locations (within deep)
  Fstats.Deep.GJ <- per_loc_fsts(coral.genind.bypop$GD, coral.genind.bypop$JD)
  Fstats.Deep.JP <- per_loc_fsts(coral.genind.bypop$JD, coral.genind.bypop$PD)
  Fstats.Deep.GP <- per_loc_fsts(coral.genind.bypop$GD, coral.genind.bypop$PD)
  
  #Between WD and other locations (within deep)
  Fstats.West.GW <- per_loc_fsts(coral.genind.bypop$GD, coral.genind.bypop$WD)
  Fstats.West.JW <- per_loc_fsts(coral.genind.bypop$JD, coral.genind.bypop$WD)
  Fstats.West.PW <- per_loc_fsts(coral.genind.bypop$PD, coral.genind.bypop$WD)
  
  # Create data frame to hold pairwise results
  Fst_compare <- cbind(Habt.G <- Fstats.Habt.G$FST,
                       Habt.J <- Fstats.Habt.J$FST,
                       Habt.P <- Fstats.Habt.P$FST,
                       Shal.GJ <- Fstats.Shal.GJ$FST,
                       Shal.JP <- Fstats.Shal.JP$FST,
                       Shal.GP <- Fstats.Shal.GP$FST,
                       Deep.GJ <- Fstats.Deep.GJ$FST,
                       Deep.JP <- Fstats.Deep.JP$FST,
                       Deep.GP <- Fstats.Deep.GP$FST,
                       West.GW <- Fstats.West.GW$FST,
                       West.JW <- Fstats.West.JW$FST,
                       West.PW <- Fstats.West.PW$FST)
  return(Fst_compare)
}
## Main code =================================================================

### Initial set-up ===========================================================

# Obtain filenames and colours through command-line arguments
args <- commandArgs(TRUE)
vcf_filename <- args[1]
pop_filename <- args[2]
col_shal <- args[3]
col_deep <- args[4]
col_west <- args[5]

# Import with vcfR
coral.genind <- vcfR2genind(read.vcfR(vcf_filename))
popfile <- read.csv(pop_filename, header = TRUE)
pop(coral.genind) <- popfile[, 2]

# Create data subsets: single pop and eastern pops (minus WD)
coral.genind.bypop <- seppop(coral.genind)
coral.genind.eastpops <- repool(coral.genind.bypop[1:6])

# Create separate dataset where missing values are replaced with mean
coral.genind.noNA <- coral.genind
coral.genind.noNA$tab <- tab(coral.genind.noNA, NA.method = "mean")
coral.genind.noNA$tab <- apply(coral.genind.noNA$tab, 1:2, as.integer)
coral.genind.bypop.noNA <- seppop(coral.genind.noNA)
coral.genind.eastpops.noNA <- repool(coral.genind.bypop.noNA[1:6])

### Calculation of different statistics ======================================

#Weir and Cockerham Fst
message("\nWeir and Cockerham Fst and Fis:")
wc(coral.genind, diploid = TRUE)

# Calculate expected heterozygosity (Hs) for each pop (genetic diversity)
message("\nExpected heterozygosity (Hs) for each pop:")
Hs(coral.genind)

# G-statistic test method for overall structure of Goudet et al. (1996)
message("\nG-statistic test for overall structure:")
gstat.randtest(coral.genind.noNA, pop = NULL, method = "global", nsim = 99)

message("\nG-statistic test for overall structure - eastern pops only:")
gstat.randtest(coral.genind.eastpops.noNA, pop = NULL, method = "global", 
               nsim = 99)

# Calculate Fst per locus for the different pairwise pop comparisons
fst_compare <- per_loc_pairwise_fsts(coral.genind.bypop)
write.csv(fst_compare, file = "afra_5a_fst_compare.csv")
boxplot(fst_compare, cex = 0.5, col=c(COL_OTHER, COL_OTHER, COL_OTHER,
                                      col_shal, col_shal, col_shal, 
                                      col_deep, col_deep, col_deep,
                                      col_west, col_west, col_west))

# Calculate mean and stdev of Fst for the different pairwise comparisons
message("\nPairwise Fst when averaged over all loci (without NA):")
apply(fst_compare, 2, function(x) mean(x, na.rm = TRUE))
apply(fst_compare, 2, function(x) sd(x, na.rm = TRUE))

# Pairwise Fst (as calculated by "pairwise.fst")
message("\nPairwise Fst:")
pairwise.fst(coral.genind)

# Non-parametric Kruskal-Wallis to test for overall differences
# (between-depths vs within-deep vs within-shallow)
message("\nNon-parametric Kruskal-Wallis to test for overall differences")
fst_compare_combine <- fst_compare
colnames(fst_compare_combine) <- c("between", "between", "between", 
                                   "shallow", "shallow", "shallow", 
                                   "deep", "deep", "deep", 
                                   "west", "west", "west")
fst_melt <- melt(fst_compare_combine)
kruskal.test(fst_melt$value ~ fst_melt$X2)

# Non-parametric Wilcoxon signed rank to test for differences in Fst
message("\nNon-parametric Wilcoxon: between-depths vs within-shallow")
wilcox.test(
  fst_melt[fst_melt$X2 == "between" | fst_melt$X2 == "shallow", ]$value ~ 
  fst_melt[fst_melt$X2 == "between" | fst_melt$X2 == "shallow", ]$X2)

message("\nNon-parametric Wilcoxon: between-depths vs within-deep")
wilcox.test(
  fst_melt[fst_melt$X2 == "between" | fst_melt$X2 == "deep", ]$value ~ 
  fst_melt[fst_melt$X2 == "between" | fst_melt$X2 == "deep", ]$X2)

message("\nNon-parametric Wilcoxon: within-shallow vs within-deep")
wilcox.test(
  fst_melt[fst_melt$X2 == "shallow" | fst_melt$X2 == "deep", ]$value ~ 
  fst_melt[fst_melt$X2 == "shallow" | fst_melt$X2 == "deep", ]$X2)

message("\nNon-parametric Wilcoxon: between-westeast vs within-deep")
wilcox.test(
  fst_melt[fst_melt$X2 == "west" | fst_melt$X2 == "deep", ]$value ~ 
  fst_melt[fst_melt$X2 == "west" | fst_melt$X2 == "deep", ]$X2)

message("\nNon-parametric Wilcoxon: between-depths vs between-westeast")
wilcox.test(
  fst_melt[fst_melt$X2 == "between" | fst_melt$X2 == "west", ]$value ~ 
  fst_melt[fst_melt$X2 == "between" | fst_melt$X2 == "west", ]$X2)
