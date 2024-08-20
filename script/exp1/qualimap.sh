#!/bin/bash

module load site-apps qualimap

qualimap bamqc -bam ../HISAT2/PG_F.sorted.bam -outdir qualimap_results
qualimap bamqc -bam ../HISAT2/PG_M.sorted.bam -outdir qualimap_results
qualimap bamqc -bam ../HISAT2/AG_F.sorted.bam -outdir qualimap_results
qualimap bamqc -bam ../HISAT2/AG_M.sorted.bam -outdir qualimap_results
