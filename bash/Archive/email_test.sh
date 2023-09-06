#!/bin/bash

#SBATCH --job-name=ctdna_prep
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --account=panm024922
#SBATCH --mail-type=BEGIN, FAIL, END
#SBATCH --mail-user=fh22528@bristol.ac.uk

echo "$(date "+%Y-%m-%d %H:%M:%S"): $SLURM_JOB_NAME start id=$SLURM_JOB_ID"

echo "This is a test script - to confirm if email notifications work" > email_test.txt

cat JOB${SLURM_JOB_ID}.out | mail -s "$SLURM_JOB_NAME Ended id=$SLURM_JOB_ID"