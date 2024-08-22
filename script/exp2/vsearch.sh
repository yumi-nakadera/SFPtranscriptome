#!/bin/bash
#SBATCH --job-name=vsearch
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
 
module load vsearch

vsearch --threads 16 --log LOGFile \
	--cluster_fast ../01_CodingRegions/trinity_combine.fasta.transdecoder.cds \
	--id 0.90 \
	--centroids centroids.fasta \
	--uc clusters.uc

