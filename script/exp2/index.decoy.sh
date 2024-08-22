#!/bin/bash
#SBATCH --job-name=index.decoy
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 12
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

module load salmon

grep "^>" </scistor/guest/yna300/LymPG/Lymnaea_stagnalis_v1.0.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat /scistor/guest/yna300/LymPG/de.novo/02_Clustering/centroids.fasta /scistor/guest/yna300/LymPG/Lymnaea_stagnalis_v1.0.fa > ref.trans.genome.fasta

salmon index -t ref.trans.genome.fasta -d decoys.txt -p 12 -i salmon_index --gencode
