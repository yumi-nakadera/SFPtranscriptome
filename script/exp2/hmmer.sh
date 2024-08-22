#!/bin/bash

module load hmmer #hmmer/3.3.2

hmmscan --cpu 16 \
       --domtblout $SLURM_SUBMIT_DIR/pfam.domtblout /scistor/guest/yna300/hmmdb/Pfam-A.hmm \
       $SLURM_SUBMIT_DIR/trinity_combine.fasta.transdecoder_dir/longest_orfs.pep
