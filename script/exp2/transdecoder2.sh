#!/bin/bash
#SBATCH --job-name=transdecoder2
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

module load TransDecoder

#make the final prediction about ORFs
TransDecoder.Predict -t /scistor/guest/yna300/LymPG/Trinity_output/trinity_combine.fasta \
        --retain_pfam_hits pfam.domtblout \
        --cpu 16

