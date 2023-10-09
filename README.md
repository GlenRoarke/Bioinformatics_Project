# ctDNA nucleosome analysis using Griffin

## Authors
Glen Roarke , Francisca Segers & Adam Chambers

## Description
Exploration of the colorectal cancer genome by liquid biopsies and Griffin nucleosome profiling of ctDNA


## Getting started

The analysis work is unavailable pending future publication.

Currently the reusable programs i developed for data manipulation to improve repeatability of the bioinformatic pipeline have been included. 

Pending publication, Markdown files exist for R code used in this project, providing examples of multivariate analysis such as PCA , complex heatmaps and correlations.

This analysis used patient information so no raw data is available to ensure privacy.

## Folder Structure

A summary of each folders purpose is summarised below.


### Preprocessing
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
  

### griffin_analysis (pending)
The main folder of analysis for griffin TFBS and TSS configurations. markdwon files with patient personal information have not been included. 

### reports (pending)
A folder containing a liquid biopsy and cfDNA review and my thesis.

### Executing programs
preprocessing_scripts can be used on the same respective outputs but may need editing based on the output.


