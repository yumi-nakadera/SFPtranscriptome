#!/bin/bash

LIST=($(echo 01 02 03))

SAM=${LIST[$SLURM_ARRAY_TASK_ID]}

cd $TMPDIR

#path to Trinity
export TRINITY=/cm/shared/package/singularity/cpu/trinityrnaseq.v2.11.0.simg

#local setting issue - language
export LC_ALL=en_US.UTF-8

# copy input
cp $SLURM_SUBMIT_DIR/Ls_${SAM}.1_trimmed.fastq.gz .
cp $SLURM_SUBMIT_DIR/Ls_${SAM}.2_trimmed.fastq.gz .

singularity exec --bind /scratch:/scratch $TRINITY Trinity \
 --seqType fq --left $TMPDIR/Ls_${SAM}.1_trimmed.fastq.gz --right $TMPDIR/Ls_${SAM}.2_trimmed.fastq.gz \
 --max_memory 120G --CPU 16 --output $SLURM_SUBMIT_DIR/Trinity_output/trinity_${SAM} --full_cleanup
