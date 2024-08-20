#!/bin/bash

#load HISAT2
module load site-apps HISAT2

#Building an index
hisat2-build centroids.fasta ref.PG

#Aligning
##PG_F
hisat2 -q -x ref.PG -1 ../raw.data/BYD_AIOSRB_7_1_HJLK5BBXX.12BA019_noribo_clean.fastq -2 ../raw.data/BYD_AIOSRB_7_2_HJLK5BBXX.12BA019_noribo_clean.fastq -S PG_F.sam

##PG_M
hisat2 -q -x ref.PG -1 ../raw.data/BYD_ALOSRB_7_1_HJLK5BBXX.12BA021_noribo_clean.fastq -2 ../raw.data/BYD_ALOSRB_7_2_HJLK5BBXX.12BA021_noribo_clean.fastq -S PG_M.sam

##AG_F
hisat2 -q -x ref.PG -1 ../raw.data/BYD_AKOSRB_7_2_HJLK5BBXX.12BA020_noribo_clean.fastq -2 ../raw.data/BYD_AKOSRB_7_1_HJLK5BBXX.12BA020_noribo_clean.fastq -S AG_F.sam

##AG_M
hisat2 -q -x ref.PG -1 ../raw.data/BYD_AMOSRB_7_1_HJLK5BBXX.12BA022_noribo_clean.fastq -2 ../raw.data/BYD_AMOSRB_7_2_HJLK5BBXX.12BA022_noribo_clean.fastq -S AG_M.sam
