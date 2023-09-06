
curl --remote-name --remote-time ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz

# initialise conda
source ~/initConda.sh

# run ctDNA preperation

sbatch ctdna_prep.slurm /user/work/fh22528/Illumina_75bp_paired_end/fastq Resources
# job id - 11456109

sbatch ctdna_prep.slurm /user/work/fh22528/Illumina_150bp_paired_end/fastq Resources
# job id -

# example
sbatch ctdna_prep.slurm /user/work/ec20449/ctDNA/fastq Resources

# recalibration step - error
sbatch 10-05-ctdna_prep_recal.slurm /user/work/fh22528/Illumina_75bp_paired_end/fastq Resources
# job id - 11461551

#create a list of bam files in a variable

# loop through this list and print the name.

#!/bin/bash
# try to loop through corrected fasta files process through transdecoder.
# we could look to move each file do its own folder not sure if needed











