#!/bin/bash


module load TransDecoder

#make the final prediction about ORFs
TransDecoder.Predict -t /scistor/guest/yna300/LymPG/Trinity_output/trinity_combine.fasta \
        --retain_pfam_hits pfam.domtblout \
        --cpu 16

