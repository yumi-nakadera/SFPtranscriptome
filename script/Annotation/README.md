# Annotation
Based on both experiments, we carried out functional annotation to select candidate SFPs. 

# Data
PG specific transcripts from Exp1 

# Methods 
1. run Trinotate, SignalP and TargetP (Trinotate https://github.com/Trinotate/Trinotate.github.io/wiki)
2. select the target transcripts based on signal peptides and tissue specificity 
3. run InterProScan on target proteins
```
docker run --rm \
--platform linux/amd64 \
-v $PWD/interproscan-5.67-99.0/data:/opt/interproscan/data \
-v $PWD/output:/output \
-v $PWD/temp:/temp \
-v $PWD/input:/input \
interpro/interproscan:5.67-99.0 \
--input /input/target_sfp.fasta \
--disable-precalc \
--output-dir /output \
--tempdir /temp \
--cpu 8
```
4. run eggNOG-mapper on target proteins (eggnog.sh)

    
