#!/bin/bash

#SBATCH --job-name=ctdna_prep
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --account=panm024922

#Use: sbatch ctdna_prep.slurm <fastq_file_location> <path_to_Resources_dir>

BEGIN=$(date)

module load apps/fastp/0.21.0
module load apps/bwa/0.7.17
module load apps/samtools/1.12
module load apps/gatk/4.2.3.0

#Get location of folder with fastq files
fastq_files=$1

#Get location of folder with genome resources
resources=$2

#Set location where output folders for large files are placed to one directory up from the fastq_files
work_dir=$1/..

#Make fastp output directories and run program
mkdir ${work_dir}/trimmed_reads

mkdir fastp_reports
for i in $(ls ${fastq_files}/*R1*)
do
    r1=$(basename $i)
    extension=${r1##*R1}
    name=$(basename $i _R1${extension})
    r2=${name}_R2${extension}
    fastp -w 12 --cut_front --cut_tail --trim_front1 1 --trim_front2 1 \
    --html fastp_reports/${name}.html --json fastp_reports/${name}.json \
    -i ${fastq_files}/${r1} -I ${fastq_files}/${r2} \
    -o ${work_dir}/trimmed_reads/${r1} -O ${work_dir}/trimmed_reads/${r2}
done

#Make output directory for read alignments and map reads
#A read group is added with the "-R" flag because this is a requirement for base recalibration
#For more info on read groups see: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
mkdir ${work_dir}/bam_files
for i in $(ls ${work_dir}/trimmed_reads/*R1*)
do
    r1=$(basename $i)
    extension=${r1##*R1}
    sample=${r1%_L00*}
    name=$(basename $i _R1${extension})
    r2=${name}_R2${extension}
#Get flowcell and lane number for read group from the r1 file:
    FC=$(zcat $i | head -1 | cut -d':' -f3)
    L=$(zcat $i | head -1 | cut -d':' -f4)
#Generate read group string:
    RG=@RG\\tID:${FC}.${L}\\tLB:LIB-${sample}\\tPL:ILLUMINA\\tSM:${sample}
#Map and pipe output into samtools to sort the .bam files by coordinate
    bwa mem -t 12 ${resources}/hg38.fa \
    -R $RG \
    ${work_dir}/trimmed_reads/${r1} \
    ${work_dir}/trimmed_reads/${r2} | \
    samtools sort -@ 11 -o ${work_dir}/bam_files/${name}.bam -
done

#Make output directories for deduplicated alignment files and map statistics
mkdir ${work_dir}/dedup_bam_files
mkdir map_statistics

#Create a list with sample names from the files in the "bam_files" directory:
sample_list=$(ls ${work_dir}/bam_files | sed 's:^.*/::' | awk -F '_L' '{print $1}' | sort | uniq)
#echo $sample_list

#Remove duplicated reads from .bam files
for s in $sample_list
do
    input=$(ls ${work_dir}/bam_files/${s}* | sed 's:^/:I=/:g')
    gatk MarkDuplicates REMOVE_DUPLICATES=TRUE \
    $input O=${work_dir}/dedup_bam_files/${s}.bam METRICS_FILE=map_statistics/${s}_dedup.txt
    gatk CollectAlignmentSummaryMetrics \
    R=${resources}/hg38.fa I=${work_dir}/dedup_bam_files/${s}.bam O=map_statistics/${s}_map.txt
done

#Make output directories for recalibration tables and recalibrated alignments
mkdir ${work_dir}/recal_bam_files

#Recalibrate the bases in the .bam files:

parallel --jobs 12 \
gatk --java-options "-Xmx4g" BaseRecalibrator \
-R ${resources}/hg38.fa -I {} \
--known-sites ${resources}/Mills_and_1000G_gold_standard.indels.hg38.vcf \
--known-sites ${resources}/All_20180418.vcf \
-O {.}_recal.txt ::: ${work_dir}/dedup_bam_files/*.bam

parallel --jobs 12 \
gatk --java-options "-Xmx4g" ApplyBQSR \
-R ${resources}/hg38.fa -I {} \
-bqsr {.}_recal.txt -O {.}_recal.bam ::: ${work_dir}/dedup_bam_files/*.bam

mv ${work_dir}/dedup_bam_files/*recal* ${work_dir}/recal_bam_files

#Index the recalibrated bam files for Griffin and for wig file creation with hmmcopy_utils/bin/readCounter
for i in $(ls ${work_dir}/recal_bam_files/*.bam)
do
    samtools index -@ 11 $i
done

END=$(date)

echo Begin:$BEGIN, End:$END > time.txt

# rename time file with date to stop overwriting GR - 05/05/2023
mv time.txt time_$(date +%Y-%m-%d_%H%M%S).txt
