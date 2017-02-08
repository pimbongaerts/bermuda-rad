#!/usr/bin/env python
"""

Characterize outliers generated with outlier_summary.py using custom criteria

"""


__author__ = "Pim Bongaerts"
__copyright__ = "Copyright (C) 2015 Pim Bongaerts"
__license__ = "GPL"

import sys

UNDETERMINED = 'UNDETERMINED'
DEPTH_OUTLIER = 'DEPTH_OUTLIER'
FDIST_AND_BAYESCAN_OUTLIER = 'FDIST_AND_BAYESCAN_OUTLIER'
FDIST_OUTLIER = 'FDIST_OUTLIER'
BAYESCAN_OUTLIER = 'BAYESCAN_OUTLIER'
NOT_OUTLIER = 'NON_OUTLIER'

TRUE_STR = 'True'
FALSE_STR = 'False'

def main(summary_filename):
    summary_file = open(summary_filename, 'r')
    for line_index, line in enumerate(summary_file):
        if line_index == 0:
            # Create header and output to STDOUT
            header = line.split()
            concat_header = '\t'.join(header)
            sys.stdout.write('{0}\toutlier_type\n'.format(concat_header))
        else:
            # Load outlier assessment outcomes in dictionary
            outlier = {}
            outlier_type = UNDETERMINED
            for col_index, col in enumerate(line.split()):
                if col_index < 5:
                    outlier[header[col_index]] = col
                else:
                    if col.strip() == TRUE_STR:
                        outlier[header[col_index]] = True
                    elif col.strip() == FALSE_STR:
                        outlier[header[col_index]] = False
                    else:
                        sys.exit('Invalid value: {0}'.format(col))

            # Evaluate type of outlier
            if (outlier['fdist_overall'] and outlier['bayescan_overall']):
                if sum([outlier['fdist_PS-PD'], outlier['fdist_JS-JD'], 
                       outlier['fdist_GS-GD']]) >= 2 and \
                   not (outlier['fdist_shallow'] or outlier['fdist_deep']):
                    outlier_type = DEPTH_OUTLIER
                else:
                    outlier_type = FDIST_AND_BAYESCAN_OUTLIER
            elif outlier['fdist_overall']:
                outlier_type = FDIST_OUTLIER
            elif outlier['bayescan_overall']:
                outlier_type = BAYESCAN_OUTLIER
            else:
                outlier_type = NOT_OUTLIER
            
            # Add outlier_type to table and output to STDOUT
            sys.stdout.write('{0}\t{1}\n'.format(line.strip('\n'), 
                                                 outlier_type))
                
if __name__ == '__main__':
    main(sys.argv[1])