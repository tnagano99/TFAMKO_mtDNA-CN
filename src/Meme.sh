#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=def-ccastel
#SBATCH --mem=50G
#SBATCH --nodes=1


echo "Run Meme Suite Analysis script"

cd ../data/fasta

streme -oc /home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/streme_out_500 -objfun de -dna -minw 8 -maxw 20 -p sig_genes_upstream.fa 