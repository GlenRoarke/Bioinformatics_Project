# ctDNA nucleosome analysis using Griffin

## Authors
Glen Roarke , Francisca Segers & Adam Chambers

## Description
Exploration of the colorectal cancer genome by liquid biopsies and Griffin nucleosome profiling of ctDNA


## Getting started
Markdown files exist for R code used in this project, providing examples of multivariate analysis such as PCA , complex heatmaps and correlations.

This analysis used patient information so no raw data is available to ensure privacy.

## Folder Structure

A summary of each folders purpose is summarised below.

### ctDNA_prep
The ctDNA preparation step is run on the raw fastq reads, to trim low quality reads, map them to the human reference genome and recalibrate the bam files.
Once this is completed the preprocessing scripts can be run in order to summarise the results in a table.

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

### Griffin_workflow
Example of the open source software griffin and the custom snakemake developed griffin_nucleosome_profiling_tss.snakefile for Transcription start site analysis.   

### griffin_analysis
The main folder of analysis for griffin TFBS and TSS configurations. markdwon files with patient personal information have not been included. 

### reports 
A folder containing a liquid biopsy and cfDNA review and my thesis.

### Executing programs.
preprocessing_scripts can be used on the same respective outputs but may need editing based on the output.


