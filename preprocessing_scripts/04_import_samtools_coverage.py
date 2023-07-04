import pandas as pd
import argparse
import os
from datetime import datetime
import glob

# Create the argparse object and define the command-line arguments
parser = argparse.ArgumentParser(description='Import results from samtools coverage')
parser.add_argument('input_directory', type=str, help='input directory path containing *.txt files')

results = []

# Parse the arguments
args = parser.parse_args()

# Retrieve *.txt files from the input directory
file_paths = glob.glob(os.path.join(args.input_directory, '*.txt'))

# loop through file paths and import 
for file_path in file_paths:
    with open(file_path, 'r') as f:
        df = pd.read_csv(f, sep='\t')
        # remove mitochondrial DNA, chromosome variations, and sex chrs
        df = df[~df['#rname'].str.contains("random|alt|chrM|chrUn_|chrY|chrX", regex=True)]
        filename_base = os.path.basename(file_path)
        filename_without_ext = os.path.splitext(filename_base)[0]
        df['filename'] = filename_without_ext
        df['sample_id'] = df['filename'].str.extract(r'\d*-(.*?)_')
        agg_df = df.groupby(['sample_id', 'filename'])[['meandepth', 'coverage']].mean()
        results.append(agg_df)

# Combine the results into a single data frame and reorder
merged_df = pd.concat(results)

# Add datetime to output filename
now = datetime.now()
dt_string = now.strftime("%Y-%m-%d_%H-%M-%S")
output_filename = f"samtools_stats_{dt_string}.csv"

merged_df.to_csv(output_filename, index=True)

print(f'The samtools mapping stats have been saved to the file {output_filename}')
