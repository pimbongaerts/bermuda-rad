# DAPC from VCF file
#
# NOTE: specific to analysis detailed in github.com/pimbongaerts/bermuda_rad
#       DAPC settings are hard-coded
#
# Authors: Cynthia Riginos, Pim Bongaerts
# Copyright (C) 2015 Cynthia Riginos & Pim Bongaerts
# License: GPL

## Dependencies ==============================================================

suppressMessages(library("adegenet"))
suppressMessages(library("vcfR"))

## Main code =================================================================

# Obtain filenames through command-line arguments
args <- commandArgs(TRUE)
vcf_filename <- args[1]
pop_filename <- args[2]

# Import with vcfR
coral.genind <- vcfR2genind(read.vcfR(vcf_filename))
popfile <- read.csv(pop_filename, header=TRUE)
pop(coral.genind) = popfile[,2]

# Perform Discriminant Analysis of Principal Components (DAPC)
max_PCAs <- as.integer(length(coral.genind$pop) / 3) # as <= N/3 advised
dapc.pop <- dapc(coral.genind, n.pca = max_PCAs, n.da = 6)
optimum_score <- optim.a.score(dapc.pop)
dapc.pop <- dapc(coral.genind, n.pca = optimum_score$best, n.da = 6)

## Set output files
filename_base <- strsplit(vcf_filename, "\\.")[[1]][[1]]
output_indcoord_filename <- paste(filename_base, "_indcoord.txt", sep = "")
output_posterior_filename <- paste(filename_base, "_posterior.txt", sep = "")

write.table(dapc.pop$ind.coord, file = output_indcoord_filename, sep = ",")
write.table(dapc.pop$posterior, file = output_posterior_filename, sep = ",")