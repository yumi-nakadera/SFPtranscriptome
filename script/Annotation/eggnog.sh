#!/bin/bash
#SBATCH --job-name=emapper
#SBATCH -p defq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=y.nakadera@vu.nl
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

#move to TMPDIR
cd $TMPDIR

# copy the data - protein seq fasta file (N=640)
cp $SLURM_SUBMIT_DIR/target.sfp .

#load eggnog-mapper/2.1.4-foss-2020b
module load site-apps
module load eggnog-mapper

# run eggNog mapper 
emapper.py -i target.sfp --itype proteins -o target.sfp_eggnog --output_dir $SLURM_SUBMIT_DIR --excel
