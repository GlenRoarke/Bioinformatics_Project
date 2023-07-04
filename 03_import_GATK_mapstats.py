# Import mapping statistics
# author: Glen Roarke
# Date: 12/05/2023
# Purpose - To import GATK mapping statistics into a data frame for further calculations.

import pandas as pd
import argparse
import os
from datetime import datetime

# use in a folder with GATK CollectAlignmentSummaryMetrics outputs 
# module load languages/anaconda3/2021-3.9-bioconda
# python import_mapstats.py *_map.txt

# Create the argparse object and define the command-line arguments
parser = argparse.ArgumentParser(description='Import mapping statistics from GATK CollectAlignmentSummaryMetrics')
parser.add_argument('filenames', type=str, nargs='+', help='input GATK map files')


results = []

# Parse the arguments
args = parser.parse_args()

# loop through parsed filenames and import 3 rows of GATK mapping stats
for filename in args.filenames:
    with open(filename, 'r') as f:
        df = pd.read_csv(f, sep='\t', skiprows=6, nrows =3, header = 0)
        sub_df = df.iloc[-1:] # subset last row of paired ends
        sub_df = sub_df.loc[:, ['TOTAL_READS','PF_READS', 'READS_ALIGNED_IN_PAIRS','PCT_READS_ALIGNED_IN_PAIRS','PF_READS_IMPROPER_PAIRS','PCT_PF_READS_IMPROPER_PAIRS']]
        filename_base = os.path.basename(filename)
        filename_without_ext = os.path.splitext(filename_base)[0]
        sub_df['filename'] = filename_without_ext
        results.append(sub_df)
        
# Combine the results into a single data frame and reorder
merged_df = pd.concat(results)


# calculate aligned read pairs
# PCT_READS_ALIGNED_IN_PAIRS - PCT_READS_ALIGNED_IN_PAIRS 
merged_df['Mapped_read_pairs'] =  merged_df.PCT_READS_ALIGNED_IN_PAIRS - merged_df.PCT_PF_READS_IMPROPER_PAIRS

#sample id 
merged_df['sample_id'] = merged_df['filename'].str.extract(r'\d*-(.*?)_')

# reorder cols
merged_df = merged_df[['sample_id','filename', 'Mapped_read_pairs','TOTAL_READS','PF_READS','READS_ALIGNED_IN_PAIRS'
,'PCT_READS_ALIGNED_IN_PAIRS','PF_READS_IMPROPER_PAIRS','PCT_PF_READS_IMPROPER_PAIRS']]
 
# Add datetime to output filename
now = datetime.now()
dt_string = now.strftime("%Y-%m-%d_%H-%M-%S")
output_filename = f"mapping_stats_{dt_string}.csv"

merged_df.to_csv(output_filename, index=False)

print(f'The GATK mapping stats have been saved to the file {output_filename}')


