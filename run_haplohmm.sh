#!/usr/bin/env bash
#SBATCH --job-name=STRAIN
#SBATCH --output=STRAIN.log
#SBATCH --mem=8GB
#SBATCH -n 2
#SBATCH --time=1:00:00

cd STRAIN/
module load singularity
singularity exec /projects/chesler-lab/phenome/snp_grid/v2/code/sif/python36.sif python haplohmm.py 

