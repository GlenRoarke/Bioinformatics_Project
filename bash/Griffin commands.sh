### Griffin

source ~/initConda.sh
conda activate griffin

# Griffin GC and mappability correction 
# preprocessing - file management - GC mappability correction 

# within current analysis cp required files and folders.
cd analysis # input -variable

#move or copy recal bam files
mv recal_bam_files/ ../griffin_workflow/analyis_15-05-2023/ 

# copy griffin GC mappability files 
cp -r ../snakemakes/griffin_GC_and_mappability_correction .
cp ../create_sample_file.py .
cp ../GC_and_mappability_correction.slurm griffin_GC_and_mappability_correction

# create sample file - could change output
python create_sample_file.py /user/work/fh22528/griffin_workflow/analyis_15-05-2023/recal_bam_files 
python create_sample_file.py /user/work/fh22528/griffin_workflow/analysis_b1_24-05-2023/recal_bam_files

python create_sample_file.py /user/work/fh22528/griffin_workflow/analysis_b2_30-05-2023/recal_bam_files

cp samples.yaml griffin_GC_and_mappability_correction/config

#execute
sbatch GC_and_mappability_correction.slurm
# job-id - 11524757

### nucleosome profiling ####

#copy nucleosome profiling snakemake
cp -r ../snakemakes/griffin_nucleosome_profiling .

# copy the GC yaml output into config
cp griffin_GC_and_mappability_correction/results/samples.GC.yaml \
griffin_nucleosome_profiling/config


#create sites
#move into griffin workflow ref
cd ../Ref/


python create_sites_file.py /user/work/fh22528/griffin_workflow/Ref/30000_unfiltered_sites_CIS_BP_v2

cp sites.yaml ../analysis_b2_30-05-2023/griffin_nucleosome_profiling/config/

# copy the slurm script 
cd ../Analysis1/griffin_nucleosome_profiling

#copy the slurm script into this folder
cp ../../nucleosome_profiling_v2.slurm .

#execute nucleosome profiling

source ~/initConda.sh
conda activate griffin


sbatch nucleosome_profiling_v2.slurm
# job id 11525452