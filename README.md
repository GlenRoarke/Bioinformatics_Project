# Preprocessing of Circulating tumour DNA

## Authors
Glen Roarke

This set of python and shell scripts preprocess the outputs of bioinformatics preprocessing scripts.

## Description

## Getting started

ctDNA prep script must be run first on BlueCrystal HPC. Folders and file paths need to be defined.

### Dependencies

01-samTools_coverage_parallel.slurm
This script runs the samtools coverage command in parallel on each ctDNA prepared fastq file.

02_read_fastP_stats.py
This script uses pandas to import fastP json results into a tabular format.

03_import_GATK_mapstats.py
This script imports the results of the GATK outputs into a tabular format.

04_import_samtools_coverage.py
This script imports the results from 01-samTools_coverage_parallel.slurm into a tabular format.
      
05_merge_pre_results.py
This script uses pandas to combine all the preprocessing outputs into one table. Each file has been assigned a filename from a sample.
               
### Executing programs.
All the python scripts make use of the argparse module and have instructions for inputs.
