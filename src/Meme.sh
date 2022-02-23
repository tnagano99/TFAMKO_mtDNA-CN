#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ccastel
#SBATCH --mem=240G
#SBATCH --nodes=1


echo "Run Meme Suite Analysis script"
meme -oc /home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/results/meme_out -neg /home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/data/fasta/neg_controls.fa -objfun de -rna -mod oops -minw 5 -maxw 20 /home/tnagano/projects/def-ccastel/tnagano/TFAMKO_mtDNA-CN/data/fasta/tfam_ko.fa