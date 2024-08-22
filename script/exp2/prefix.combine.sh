#!/bin/bash

##give prefix for all the Trintiy files

SAM=01
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=02
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=03
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=04
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=05
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=06
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=07
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=08
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=09
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=10
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=11
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=12
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=13
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=14
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta
SAM=15
sed "s/>/>${SAM}_/g" ../Trinity_output/trinity_${SAM}.Trinity.fasta > ../Trinity_output/trinity_prefix_${SAM}.Trinity.fasta

#combined all the renamed files and place it the working dir
cat ../Trinity_output/trinity_prefix_* > ../Trinity_output/trinity_combine.fasta
