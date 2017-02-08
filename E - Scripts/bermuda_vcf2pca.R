# PCA from VCF file
#
# NOTE: specific to analysis detailed in github.com/pimbongaerts/bermuda_rad
#       PCA settings are hard-coded
#
# Authors: Cynthia Riginos, Pim Bongaerts
# Copyright (C) 2015 Cynthia Riginos & Pim Bongaerts
# License: GPL

## Dependencies ==============================================================

library("adegenet")
library("vcfR")

## Main code =================================================================

## Obtain filenames through command-line arguments
args <- commandArgs(TRUE)
vcf_filename <- args[1]
pop_filename <- args[2]

## Import with vcfR
coral.genind <- vcfR2genind(read.vcfR(vcf_filename))
popfile <- read.csv(pop_filename, header = TRUE)
pop(coral.genind) = popfile[,2]

## Set output files
filename_base <- strsplit(vcf_filename, "\\.")[[1]][[1]]
output_pca_filename <- paste(filename_base, "_pca.txt", sep = "")
output_pca_png_filename <- paste(filename_base, "_pca.png", sep = "")
output_eigenval_png_filename <- paste(filename_base, "_eigenval.png", 
                                      sep = "")

## Principal component analysis
coral.scaled <- scaleGen(coral.genind, center=TRUE, scale=TRUE, 
                         NA.method="mean")  
pca1 <- dudi.pca(coral.scaled, cent = FALSE, scale = FALSE, scannf = FALSE, 
                 nf = 2)  # 2 PCs retained
write.table(pca1$l1, file = output_pca_filename)

# Plot - eigenvalues
png(filename = output_eigenval_png_filename)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
dev.off()

# Plot - PCA
png(filename = output_pca_png_filename)
plot(pca1$l1, cex = 1.5, col = as.character(popfile$PopColor), 
     bg = as.character(popfile$PopColor), pch = popfile$pch)
dev.off()