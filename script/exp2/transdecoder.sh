#!/bin/bash
#SBATCH --job-name=transdecoder
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

module load TransDecoder #TransDecoder/5.5.0

#identify long ORFs (min 100 AA long)
TransDecoder.LongOrfs -t ../Trinity_output/trinity_combine.fasta



