# Read_counts
# Purpose: Generate raw read count and trimmed read count from fastp json output
# auther: glen roarke
# date : 18/05/2023
# updated : 21/05/2023

import pandas as pd
import argparse
import os
from datetime import datetime
import json

# use in a folder with fastp json outputs 
# module load languages/anaconda3/2021-3.9-bioconda
# python import_mapstats.py *.json

# Create the argparse object and define the command-line arguments
parser = argparse.ArgumentParser(description='Import read statistics from fastp json outputs')
parser.add_argument('filenames', type=str, nargs='+', help='input json file list')


results = []

# Parse the arguments
args = parser.parse_args()


# loop through parsed filenames and import 3 rows of GATK mapping stats
for filename in args.filenames:
    with open(filename, 'r') as f:
        json_data = json.load(f)
        # Extract the "summary" data
        summary_data = json_data['summary']
        # Create a DataFrame from the "summary" data
        df = pd.DataFrame(summary_data).T.reset_index()
        filename_base = os.path.basename(filename)
        df['filename'] = filename_base
        results.append(df)
        
# Combine the results into a single data frame and reorder
merged_df = pd.concat(results)

# derive sample ID from the filename
merged_df['sample_id'] = merged_df['filename'].str.extract(r'\d*-(.*?)_')

# reorder cols
merged_df = merged_df[['filename', 'index','sample_id','total_reads','total_bases','q20_bases','q30_bases','q20_rate','q30_rate',
'read1_mean_length','read2_mean_length','gc_content']]

#group by sample and before and after trimming (old index)
agg_df  = merged_df.groupby(['sample_id','index'])[['total_reads']].sum() / 2


# Pivot the DataFrame
df_pivot = agg_df.pivot_table(index='sample_id' , columns='index', values='total_reads').reset_index()

df_pivot = df_pivot[['sample_id', 'before_filtering', 'after_filtering']]

# Add datetime to output filename
now = datetime.now()
dt_string = now.strftime("%Y-%m-%d_%H-%M-%S")
output_filename = f"read_stats_{dt_string}.csv"

print(f' Read statisitcs have been saved to the file {output_filename}.')

df_pivot.to_csv(output_filename, index=False)




