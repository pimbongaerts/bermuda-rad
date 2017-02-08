#!/usr/bin/env python
"""
Extract HWE test results from Arlequin XML output file

Goal is to screen for suspicious SNPs that show significant deviations from 
HWE across multiple pops, or that show a significant Heterozygote excess for 
at least one population.

Note: order of loci in optional VCF needs to match that of arlequin input

Usage: arlequinHWEsum.py -i <arlequin file> 
                         [-p <P-value cut-off>] 
                         [-d <number of HWE deviations cut-off>] 
                         [-v <corresponding vcf file>]
                         
NOTE: specific to analysis detailed in github.com/pimbongaerts/bermuda_rad 
      also: terrible script - in need of rewrite/refactor.
"""

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2014 Pim Bongaerts'
__license__ = 'GPL'

import sys
import re
import getopt

USAGE_MESSAGE = ('arlequinHWEsum.py -i <arlequin file> [-p <P-value cut-off>]'
                 ' [-d <number of HWE deviations cut-off>] [-v <corresponding'
                 ' vcf file>]')
LANDMARK = '== Hardy-Weinberg equilibrium'
LANDMARK_POPULATION_SEPARATOR = ':'
HEADER_VCF = '#'
IDENT_MONOMORPHIC = 'This'
TEXT_EMPTY = ''
TEXT_MONOMORPHIC = 'm'
TEXT_SIGN = 'S'
TEXT_HETEXCESS = 'e'
TEXT_HETDEFICIT = 'd'
TEXT_NOT_SIGN = ''
TEXT_OVERALL_SIGNIFICANT = 'SIGN'
TEXT_OVERALL_MONOMORPHIC = 'FIXED'
TEXT_OVERALL_HETEXCESS = 'EXCESS'
SHORT_OUTPUT_COL_LOCUS = 'LOCUS'
SHORT_OUTPUT_COL_SUMMARY = 'SUMMARY'
DETAILED_OUTPUT_COL_LOCUS = 'LOCUS'
DETAILED_OUTPUT_POP_COLS = ['LOCUS' , '#Genot' , 'Obs.Het.' , 
                            'Exp.Het. ', 'P-value' , 'Sign']
COL_LOCUS = 0
COL_NG = 1
COL_OBSHET = 2
COL_EXPHET = 3
COL_PVAL = 4
VCF_COL_CHROM = 0
VCF_COL_POS = 1
 
def represents_int(string):
    """ Return True if string can be converted into integer """
    try:
        int(string)
        return True
    except ValueError:
        return False

def get_command_line_argv(argv):
    """ Obtain variables from command-line arguments or exit if any issues """
    # Required arguments
    arlequin_file = ''
    # Optional argument defaults
    pvalue_cutoff = 0.05
    N_sign_cut_off = 1
    vcf_file = ''
    
    try:
        opts, args = getopt.getopt(argv, 'hi:p:d:v:',
                                   ['input=','pvalue=','devs=','vcf'])
    except getopt.GetoptError:
        print(USAGE_MESSAGE)
        sys.exit(2)            # Command-line syntax error (2)
    for opt, arg in opts:
        if opt == '-hi':
            print(USAGE_MESSAGE)
            sys.exit()        # Succesful termination (0)
        elif opt in ('-i', '--input'):
            arlequin_file = arg
        elif opt in ('-p', '--pvalue'):
            pvalue_cutoff = arg
        elif opt in ('-d', '--devs'):
            N_sign_cut_off = arg
        elif opt in ('-v', '--vcf'):
            vcf_file = arg
    if not arlequin_file:
        print(USAGE_MESSAGE)
        sys.exit(2)            # Command-line syntax error (2)
    return arlequin_file, pvalue_cutoff, N_sign_cut_off, vcf_file

def get_locus_dict_from_vcf(vcf_file):
    """ Creat dict with original CHROM:POS names from VCF file """
    inputFile = open(vcf_file, 'r')
    vcf_loci = {}
    loci_ID = 1    
    for line in inputFile:
        cols = line.split('\t')
        if not line[0][0] == HEADER_VCF:
            vcf_loci[loci_ID] = '{0}:{1}'.format(cols[VCF_COL_CHROM],
                                                 cols[VCF_COL_POS])
            loci_ID += 1
    return vcf_loci

def extract_population_name(landmark_line):
    """ Extract population name/code from landmark line """
    landmark_cols = landmark_line.split(LANDMARK_POPULATION_SEPARATOR)
    population = re.sub('[()]', '', landmark_cols[1])
    return population.strip()

def get_VCF_locus_name(arlequin_locus, vcf_loci):
    """ Return original (+ arlequin) locus name if VCF file is included  """    
    if vcf_loci:
        locus_name = '{0} [{1}]'.format(vcf_loci[arlequin_locus], 
                                        arlequin_locus)
    else:
        locus_name = str(arlequin_locus)
    return locus_name

def create_detailed_output_file_header(populations):
    """ Create output file header based on number of populations """
    header_columns = []
    header_columns.append(DETAILED_OUTPUT_COL_LOCUS)
    for population in populations:
        for output_col in DETAILED_OUTPUT_POP_COLS:
            header_columns.append('{0}:{1}'.format(population,output_col))
    return ','.join(header_columns)

def create_short_output_file_header(populations):
    """ Create output file header based on number of populations """
    header_columns = []
    header_columns.append(SHORT_OUTPUT_COL_LOCUS)
    for population in populations:
        header_columns.append(population)
    header_columns.append(SHORT_OUTPUT_COL_SUMMARY)
    return ','.join(header_columns)

def output_detailed_file(output_filename, hwe_info_detailed, vcf_loci):
    """ Create detailed output file """
    output_file = open(output_filename, 'w')
    header = create_detailed_output_file_header(hwe_info_detailed.keys())
    output_file.write('{0}\n'.format(header))
    cur_pop = next(iter(hwe_info_detailed.keys()))
    
    # Iterate through loci/pops (using cur_pop to retrieve loci)
    for locus in hwe_info_detailed[cur_pop]:
        # Concatenate information from one locus into a single-line string
        hwe_info_detailed_pop = []
        
        # Include VCF CHROM/POS if given
        locus_name = get_VCF_locus_name(locus, vcf_loci)
        hwe_info_detailed_pop.append(locus_name)
        
        # Iterate over populations
        for population in hwe_info_detailed:
            hwe_info_detailed_pop.append(hwe_info_detailed[population][locus])
        concat_line = ','.join(hwe_info_detailed_pop)
        output_file.write(concat_line + '\n')
    output_file.close()    

def output_short_file(output_filename, output_elimpos_filename, 
                      output_elimchrom_filename, hwe_info_short, 
                      N_sign_cut_off, vcf_loci):
    """ Create short output file  """
    # Open file and generate header
    output_file = open(output_filename, 'w')
    header = create_short_output_file_header(hwe_info_short.keys())
    output_file.write('{0}\n'.format(header))
    output_elimpos_file = open(output_elimpos_filename, 'w')
    output_elimchrom_file = open(output_elimchrom_filename, 'w')
    cur_pop = next(iter(hwe_info_short.keys()))
    
    # Initialise counters
    overall_counter = 0
    overall_significant_counter = 0
    overall_monomorphic_counter = 0
    overall_hetexcess_counter = 0
    overall_hetdeficit_counter = 0
    
    # Iterate through loci/pops (using cur_pop to retrieve loci)
    for locus in hwe_info_short[cur_pop]:
        # Concatenate information from one locus into a single-line string
        hwe_info_short_pop = []
        for population in hwe_info_short:
            hwe_info_short_pop.append(hwe_info_short[population][locus])
        concat_line = ','.join(hwe_info_short_pop)
        
        # Include VCF CHROM/POS if given
        locus_name = get_VCF_locus_name(locus,vcf_loci)
        
        # Evaluate if locus has more sign. deviations from HWE than cut-off
        if concat_line.count(TEXT_SIGN) >= int(N_sign_cut_off):
            summary_text = TEXT_OVERALL_SIGNIFICANT
            # Only count if locus does not exhibit sign. heterozygote excess
            # (as counted separately)
            if concat_line.count(TEXT_SIGN+TEXT_HETEXCESS) == 0:
                overall_significant_counter += 1
                
        # Evaluate if locus is fixed (monomorphic for all loci)
        elif concat_line.count(TEXT_MONOMORPHIC) == len(hwe_info_short):
            summary_text = TEXT_OVERALL_MONOMORPHIC
            overall_monomorphic_counter += 1
        # Otherwise
        else:
            summary_text = TEXT_NOT_SIGN
        
        # Evaluate if site exhibits a sign. heterozygote excess in a pop
        if concat_line.count(TEXT_SIGN+TEXT_HETEXCESS) > 0:
            summary_text += TEXT_OVERALL_HETEXCESS
            overall_hetexcess_counter += 1
            
          # Evaluate if site exhibits a sign. heterozygote deficit in a pop
        if concat_line.count(TEXT_SIGN+TEXT_HETDEFICIT) > 0:
            overall_hetdeficit_counter += 1      
        
        # Output entire CHROM for SNPs that have sign. excess in >= 1 pop
        if TEXT_SIGN+TEXT_HETEXCESS in concat_line:
            # Output to screen
            print('{0}\t[{1}]\t{2}'.format(locus_name, concat_line,
                                           summary_text))
            # Output to eliminate file
            if vcf_loci:
                chrom = locus_name.split()[0].split(':')[0]
                pos = locus_name.split()[0].split(':')[1]
                output_elimchrom_file.write('{0}\n'.format(chrom))
            else:
                output_elimchrom_file.write(locus_name + '\n')
        
        elif TEXT_OVERALL_SIGNIFICANT in summary_text:
            # Output to screen
            print('{0}\t[{1}]\t{2}'.format(locus_name, concat_line,
                                           summary_text))
            # Output to eliminate file
            if vcf_loci:
                chrom = locus_name.split()[0].split(':')[0]
                pos = locus_name.split()[0].split(':')[1]
                output_elimpos_file.write('{0}\t{1}\n'.format(chrom, pos))
            else:
                output_elimpos_file.write(locus_name + '\n')        
        
        output_file.write('{0},{1},{2}\n'.format(locus_name, concat_line,
                                                 summary_text))
        overall_counter += 1
    
    print('\n')
    print('Total of {0} SNPs evaluated'.format(overall_counter))
    print(('{0} SNPs monomorphic across '
           'all populations').format(overall_monomorphic_counter))
    print(('{0} SNPs with a significant deficit in heterozygosity'
           ' in at least one population').format(overall_hetdeficit_counter))
    suspic_loci_count = overall_significant_counter+overall_hetexcess_counter
    print('\nSuspicious SNPs/loci ({0}):'.format(suspic_loci_count))
    print(('{0} SNPs with a significant excess in heterozygosity in'
           ' at least one population').format(overall_hetexcess_counter))
    print(('{0} addditional SNPs significantly deviating from HWE in {1}'
           ' or more populations').format(overall_significant_counter, 
                                          N_sign_cut_off))
    output_file.close()
    output_elimpos_file.close()
    output_elimchrom_file.close()

def get_arlequin_HWE_data(arlequin_file, pvalue_cutoff):
    """ Extract data from Arlequin XML file """
    # Open input/output files
    inputFile = open(arlequin_file, 'r')

    # Set dictionaries
    hwe_info_short = {}
    hwe_info_detailed = {}
        
    # Set iteration variables
    reached_landmark = False
    reached_start_of_table = False
    cur_pop = ''
    
    # Iterate through arlequin XML and extract HWE data
    for line in inputFile:
        # Split line into columns (separated by whitespace)
        stripped_line = line.strip()
        cols = stripped_line.split()
                
        if stripped_line == '':
            # Only evaluate line if not empty (otherwise continue)
            continue
        elif line[:len(LANDMARK)] == LANDMARK:
            # Reached landmark, so set flag and extract population name
            reached_landmark = True
            cur_pop = extract_population_name(line)
            hwe_info_short[cur_pop] = {}
            hwe_info_detailed[cur_pop] = {}
        elif reached_landmark and not reached_start_of_table:
            if represents_int(cols[0].strip()):
                # Reached start of data table, so set flag
                reached_start_of_table = True
        
        if reached_start_of_table:
            if represents_int(cols[0].strip()):
                # Read data and store
                locus = int(cols[COL_LOCUS])
                if cols[COL_NG][:len(IDENT_MONOMORPHIC)] == IDENT_MONOMORPHIC:
                    # Locus monomorphic for population - store data in library
                    hwe_info_detailed[cur_pop][locus] = ('{0}{1},{2},{3},{4}'
                                                         ',{5},{6}').format(
                                                            cur_pop, locus,
                                                            TEXT_EMPTY,
                                                            TEXT_EMPTY,
                                                            TEXT_EMPTY,
                                                            TEXT_EMPTY,
                                                            TEXT_MONOMORPHIC)
                    hwe_info_short[cur_pop][locus] = TEXT_MONOMORPHIC
                else:
                    # Locus polymorphic for population, eval if sign. Pvalue
                    if float(cols[COL_PVAL]) < float(pvalue_cutoff):
                        text = TEXT_SIGN
                        # Evaluate if excessive or deficit in heterozygotes
                        if float(cols[COL_OBSHET]) > float(cols[COL_EXPHET]):
                            text += TEXT_HETEXCESS
                        else:
                            text += TEXT_HETDEFICIT
                    else:
                        text = TEXT_NOT_SIGN
                           
                    # Store data in library
                    hwe_info_detailed[cur_pop][locus] = ('{0}{1},{2},{3},{4}'
                                                         ',{5},{6}').format(
                                                            cur_pop, locus,
                                                            cols[COL_NG],
                                                            cols[COL_OBSHET],
                                                            cols[COL_EXPHET],
                                                            cols[COL_PVAL],
                                                            text)
                    hwe_info_short[cur_pop][locus] = text
            else:
                # End of table: reinitialise flags and variable
                reached_landmark = False
                reached_start_of_table = False
    inputFile.close()
    return hwe_info_short, hwe_info_detailed
    
def main(argv):
    ## Obtain command-line arguments
    (arlequin_file, pvalue_cutoff, 
     N_sign_cut_off, vcf_file) = get_command_line_argv(argv)
    
    ## Extract HWE data from arlequin file and store in two dictionaries
    ## ('short' and 'detailed')
    hwe_info_short, hwe_info_detailed = get_arlequin_HWE_data(arlequin_file, 
                                                              pvalue_cutoff)
    
    ## Read locus name (CHROM/POS) from original VCF (if given)
    vcf_loci = get_locus_dict_from_vcf(vcf_file) if vcf_file else False
    
    ## Output file with detailed information
    output_filename = arlequin_file.replace('.xml','_detail.csv')
    output_detailed_file(output_filename, hwe_info_detailed, vcf_loci)
    
    ## Output file with summarised information
    output_filename = arlequin_file.replace('.xml','_short.csv')
    output_elimpos_filename = arlequin_file.replace('.xml','_elimpos.txt')
    output_elimchrom_filename = arlequin_file.replace('.xml','_elimchrom.txt')    
    output_short_file(output_filename, output_elimpos_filename, 
                      output_elimchrom_filename, hwe_info_short, 
                      N_sign_cut_off, vcf_loci)

if __name__ == '__main__':
   main(sys.argv[1:])