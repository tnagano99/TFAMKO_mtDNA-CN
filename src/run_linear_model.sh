#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ccastel
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH -n 16

echo "change directory to script"
cd ./projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/src

echo "Load R module"
module load R

echo "Run linear mixed model script"
Rscript linear_model_mtDNACN.R
