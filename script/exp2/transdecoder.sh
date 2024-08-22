#!/bin/bash

module load TransDecoder #TransDecoder/5.5.0

#identify long ORFs (min 100 AA long)
TransDecoder.LongOrfs -t ../Trinity_output/trinity_combine.fasta



