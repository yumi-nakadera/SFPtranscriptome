#!/bin/bash

module load salmon

grep "^>" </scistor/guest/yna300/LymPG/Lymnaea_stagnalis_v1.0.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat /scistor/guest/yna300/LymPG/de.novo/02_Clustering/centroids.fasta /scistor/guest/yna300/LymPG/Lymnaea_stagnalis_v1.0.fa > ref.trans.genome.fasta

salmon index -t ref.trans.genome.fasta -d decoys.txt -p 12 -i salmon_index --gencode
