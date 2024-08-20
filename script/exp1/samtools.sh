#!/bin/bash

samtools view -S -b AG_F.sam > AG_F.bam
samtools view -S -b AG_M.sam > AG_M.bam

#merge .bam
samtools merge PG.merged.bam PG_F.bam PG_M.bam
samtools merge AG.merged.bam AG_F.bam AG_M.bam
#samtools merge PG.merged.bam PG_F.bam PG_M.bam
#samtools merge AG.merged.bam AG_F.bam AG_M.bam

#sort .bam 
samtools sort PG.merged.bam -o PG.merged.bam
samtools sort AG.merged.bam -o AG.merged.bam

samtools sort PG_F.bam -o PG_F.sorted.bam
samtools sort PG_M.bam -o PG_M.sorted.bam
samtools sort AG_F.bam -o AG_F.sorted.bam
samtools sort AG_M.bam -o AG_M.sorted.bam
