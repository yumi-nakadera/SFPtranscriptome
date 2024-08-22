#!/bin/bash

 
module load vsearch

vsearch --threads 16 --log LOGFile \
	--cluster_fast ../01_CodingRegions/trinity_combine.fasta.transdecoder.cds \
	--id 0.90 \
	--centroids centroids.fasta \
	--uc clusters.uc

