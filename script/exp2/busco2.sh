#!/bin/bash
#SBATCH --job-name=busco2
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

singularity exec --bind /scratch:/scratch /cm/shared/package/singularity/cpu/busco-v5.0.0_cv1.sif \
busco -i ../02_Clustering/centroids.fasta -o busco2.output -m transcriptome -l metazoa_odb10
