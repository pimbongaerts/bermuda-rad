#!/usr/bin/env python
"""
Provides a summary of outlier FDist / BayeScan outlier analyses

FDist and Bayescan output files need to be supplied in a single folder, 
respectively labeled with the FDIST_PREFIX and BAYESCAN_PREFIX (see below)

Usage: python3 outlier_summary.py vcf_filename fdist_hetfst_filename 
               outlier_folder bayescan_cutoff fdist_cutoff
E.g.:  python3 outlier_summary.py afra_2f.vcf outlier_files/fdist_overall.txt 
               outlier_files 0.05 0.95

Note: order of SNPs in VCF file and outlier output files needs to be the same
"""

__author__ = "Pim Bongaerts"
__copyright__ = "Copyright (C) 2016 Pim Bongaerts"
__license__ = "GPL"

import os
import sys
from glob import glob

VCF_HEADER = '#'
VCF_CHROM_COL = 0
VCF_POS_COL = 1
FDIST_PREFIX = 'fdist_'
FDIST_HET_COL = 1
FDIST_FST_COL = 2
FDIST_P_COL = 3
BAYESCAN_PREFIX = 'bayescan_'
BAYESCAN_Q_COL = 3
UNKNOWN_FILE = 'unknown file'
SEPARATOR = '\t'

def chrompos_from_vcf_file(vcf_filename):
    """ Extract CHROM and POS for each SNP from VCF file """
    vcf_file = open(vcf_filename, 'r')
    snps = {}
    snp_index = 0
    for line in vcf_file:
        if line[0] != VCF_HEADER:
            snp_index += 1
            snps[snp_index] = []    # Create dictionary item for each SNP
            snps[snp_index].append(str(snp_index)) # index
            snps[snp_index].append(line.split()[VCF_CHROM_COL]) # CHROM
            snps[snp_index].append(line.split()[VCF_POS_COL]) # POS
    vcf_file.close()
    return snps

def hetfst_from_fdist_file(fdist_hetfst_filename, snps):
    """ Extract Fst and Het from FDist file and add to snps dict """    
    fst_het_file = open(fdist_hetfst_filename, 'r')
    for snp_index, line in enumerate(fst_het_file):
        if snp_index > 0:
            snps[snp_index].append(line.split()[FDIST_HET_COL]) # HET
            snps[snp_index].append(line.split()[FDIST_FST_COL]) # FST
    fst_het_file.close()
    return snps

def get_outliers_from_file(outlier_filename, snps, bayescan_cutoff, 
                           fdist_cutoff):
    """ Assess for each SNP whether it is an outlier """
    outlier_file = open(outlier_filename, 'r')
    for snp_index, line in enumerate(outlier_file):
        if snp_index > 0:
            if BAYESCAN_PREFIX in outlier_filename:
                value = float(line.split()[BAYESCAN_Q_COL])
                snps[snp_index].append(str(value <= bayescan_cutoff))
            elif FDIST_PREFIX in outlier_filename:
                value = float(line.split()[FDIST_P_COL])
                snps[snp_index].append(str(value > fdist_cutoff))
            else:
                snps[snp_index].append(UNKNOWN_FILE)
    return snps
    
def main(vcf_filename, fdist_hetfst_filename, outlier_folder, bayescan_cutoff, 
         fdist_cutoff):
    ## Iterate through VCF file and obtain CHROM and POS and add to snps
    snps = chrompos_from_vcf_file(vcf_filename)
    
    ## Iterate through file with Fst/Het and add to snps
    snps = hetfst_from_fdist_file(fdist_hetfst_filename, snps)
     
    ## Iterate through the different files in outlier_folder and add to snps
    filenames = glob('{0}/*'.format(outlier_folder))
    headers = ['index','chrom', 'pos', 'het', 'fst']
    for outlier_filename in filenames:
        snps = get_outliers_from_file(outlier_filename, snps, bayescan_cutoff, 
                                      fdist_cutoff)
        headers.append(os.path.basename(outlier_filename).replace('.txt', ''))
        
    ## Output summary to STDOUT
    print(SEPARATOR.join(headers))
    for snp_index in snps:
        print(SEPARATOR.join(snps[snp_index]))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), 
         float(sys.argv[5]))