# Preprocessing of Circulating tumour DNA

## Authors
Glen Roarke & Francisca Segars

This set of python and shell scripts preprocess the outputs of bioinformatics preprocessing scripts.

## Description

This workflow uses NHS patient information so no data files are included in the repo.

## Getting started

The ctDNA preparation step is run on the raw fastq reads, to trim low quality reads, map them to the human reference genome and recalibrate the bam files.

Once this is completed the preprocessing scripts can be run in order to summarise the results in a table.

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
