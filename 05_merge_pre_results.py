import pandas as pd
from datetime import datetime
import glob
import argparse

# Create the argparse object and define the command-line arguments
parser = argparse.ArgumentParser(description='Merge CSV files based on sample ID')
parser.add_argument('--folder', type=str, required=True, help='Folder directory to search for CSV files')

# Parse the arguments
args = parser.parse_args()

# Perform wildcard search for CSV files in the specified folder
csv_files = glob.glob(args.folder + '/*.csv')

# Perform wildcard search for CSV files
map_file = glob.glob(args.folder + '/mapping_stats_*.csv')
read_file = glob.glob(args.folder + '/read_stats_*.csv')
sam_file =  glob.glob(args.folder + '/samtools_stats_*.csv')

# Read the three CSV files assume there is only one of each
df1 = pd.read_csv(read_file[0])
df2 = pd.read_csv(map_file[0])
df3 = pd.read_csv(sam_file[0])


# Join the DataFrames based on 'sample ID'
merged_df = pd.merge(df1, df2, on='sample_id', how='inner')
merged_df = pd.merge(merged_df, df3, on='sample_id', how='inner')

merged_df = merged_df[['sample_id','before_filtering','after_filtering','Mapped_read_pairs','meandepth','coverage']]

# Format 'before_filtering' and 'after_filtering' columns
merged_df['before_filtering'] = merged_df['before_filtering'].astype(int)
merged_df['after_filtering'] = merged_df['after_filtering'].astype(int)

# Add datetime to output filename
now = datetime.now()
dt_string = now.strftime("%Y-%m-%d_%H-%M-%S")
output_filename = f"preprocessing_stats_{dt_string}.csv"

merged_df.to_csv(output_filename, index=False)

print(f'The preprocessing stats have been saved to the file {output_filename}')


