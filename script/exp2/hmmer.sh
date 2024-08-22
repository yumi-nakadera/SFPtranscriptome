#!/bin/bash
#SBATCH --job-name=hmmer
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

module load hmmer #hmmer/3.3.2

hmmscan --cpu 16 \
       --domtblout $SLURM_SUBMIT_DIR/pfam.domtblout /scistor/guest/yna300/hmmdb/Pfam-A.hmm \
       $SLURM_SUBMIT_DIR/trinity_combine.fasta.transdecoder_dir/longest_orfs.pep
