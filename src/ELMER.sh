#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --account=def-ccastel
#SBATCH --nodes=1
#SBATCH --mem=16G

echo "change directory to script"
cd ./projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/src

echo "Load R module"
module load R

echo "Run Elmer Analysis script"
Rscript 'ELMER-Analysis-DMP-Cont.R'