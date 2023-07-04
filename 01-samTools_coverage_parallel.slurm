#!/bin/bash
# create samtools statistics

#SBATCH --job-name=samtools_coverage
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --account=panm024922
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=fh22528@bristol.ac.uk

module load apps/samtools/1.12

# Check if the input directory is provided
if [ -z "$1" ]; then
    echo "Please provide the input directory as an argument."
    exit 1
fi

# Directory for coverage files
coverage_dir="coverage_stats"

# Create the directory if it doesn't exist
if [ ! -d "$coverage_dir" ]; then
    mkdir "$coverage_dir"
    echo "Created directory: $coverage_dir"
fi

# Export the coverage directory variable for parallel execution
export coverage_dir

# Function to calculate coverage for a BAM file
calculate_coverage() {
    bam_file="$1"
    file_name=$(basename "$bam_file" .bam)
    echo "Calculating coverage for $bam_file"
    samtools coverage "$bam_file" > "$coverage_dir/${file_name}_cov.txt"
}

# Export the function for parallel execution
export -f calculate_coverage

# Run the parallel command to process BAM files
parallel calculate_coverage ::: "$1"/*.bam

echo "The analysis has completed successfully"

