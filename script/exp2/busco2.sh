#!/bin/bash

singularity exec --bind /scratch:/scratch /cm/shared/package/singularity/cpu/busco-v5.0.0_cv1.sif \
busco -i ../02_Clustering/centroids.fasta -o busco2.output -m transcriptome -l metazoa_odb10
