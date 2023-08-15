#!/bin/bash
#SBATCH --job-name=impute_CHR_DIR
#SBATCH --output=logs/impute_CHR_DIR.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --array=1-NJOBS

module load singularity

config=../data/out/CHR/imputation_data/CHR_DIR/config.txt

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${config})

cd ../data/out/CHR/imputation_data/CHR_DIR/${sample}
echo "Processing ${sample}"
singularity exec /projects/chesler-lab/phenome/snp_grid/v2/code/sif/python36.sif python haplohmm.py

