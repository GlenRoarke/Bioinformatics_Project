
# initialise conda
source ~/initConda.sh

# run ctDNA preperation

# example
sbatch ctdna_prep.slurm /user/work/ec20449/ctDNA/fastq Resources

# 75bp - batch 2
sbatch ctdna_prep_v2.slurm /user/work/fh22528/Illumina_75bp_paired_end/fastq Resources
# job id v2 - 11497758

#150bp - batch 1
sbatch ctdna_prep_v2.slurm /user/work/fh22528/Illumina_150bp_paired_end/fastq Resources
#job id - 11495941

#batch 3 - 75bp 
sbatch ctdna_prep_v2.slurm /user/work/fh22528/b3_Illumina_75bp_paired_end_May2023/fastq Resources
#job id - 11518942

sbatch ctdna_prep_v2.slurm /user/work/fh22528/b3_Illumina_75bp_paired_end_May2023/Fastq Resources
# job id - 11519069



