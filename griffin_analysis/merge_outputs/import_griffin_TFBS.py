# Import griffin results
# date 23/05/2023
# To merge all griffin results into one dataframe
# use python .\import_griffin_TFBS.py --folder ..\..\Results\griffin_results\75bp_Illumina_results_21-05-2023\
# purpose - to merge the individual griffin TFBS outputs into one dataframe for downstream R analysis 


import glob
import pandas as pd
import argparse
import os
from datetime import datetime


# Create the argparse object and define the command-line arguments
parser = argparse.ArgumentParser(description='Merge Griffin results for each sample')
parser.add_argument('--folder', type=str, required=True, help='Folder directory to search for GC_coverage.tsv files')

# Parse the arguments
args = parser.parse_args()

results = []

# loop through root directory folders, import griffin results into pandas
# root_dir needs a trailing slash (i.e. /root/dir/)
for filename in glob.iglob(args.folder + '**/*.GC_corrected.coverage.tsv', recursive=True):
    df = pd.read_csv(filename, sep='\t')
    df['sample_id'] = df['sample'].str.extract(r'\d*-(.*?)_')
    results.append(df)
      
# Combine the results into a single data frame
merged_df = pd.concat(results)

#########

# Add datetime to output filename
now = datetime.now()
dt_string = now.strftime("%Y-%m-%d_%H-%M-%S")
#output_filename = f"Transposed_griffin_GC.coverage_{dt_string}.csv"
merged_filename = f"Griffin GC_corrected.coverage_{dt_string}.csv"

merged_df.to_csv(merged_filename, index=True)

print(f'The Griffin GC_corrected.coverage results have been saved to the file {output_filename} & {merged_filename}')
    
