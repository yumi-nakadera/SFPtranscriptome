# Exp 2
We reared the snails isolated, paired and grouped with five individuals. At the end of experiment, we collected their PG to screen differentially expressed transcripts.

## Data 
RNAseq Deposted in NCBI
The reference genome (Koene et al. in subm., NCBI)

## Methods 
### Making the reference PG transcriptome 
1. Quality control using FastQC and MultiQC
```
for file in *.fastq
do
fastqc --outdir . ../01_QC/ $file 
done

multiqc .
```

2. Run Trimmomatic to remove adaptors and run FastQC and MultiQC again 
```
for R1 in *1.fastq*
do
        R2=${R1//1.fastq/2.fastq}
        R1trimmed=${R1//.fastq/_trimmed.fastq.gz}
        R1singles=${R1//.fastq/_singles.fastq.gz}
        R2trimmed=${R2//.fastq/_trimmed.fastq.gz}
        R2singles=${R2//.fastq/_singles.fastq.gz}
        Output=${R1//1.fastq/trimmed.output.txt}
        java -jar ~/Software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 -phred33 $R1 $R2 $R1trimmed $R1singles $R2trimmed $R2singles ILLUMINACLIP:/Users/ynakadera/Software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads 
        LEADING:3 TRAILING:3 MINLEN:36 &> $Output
done 
```
3. Run Trinity for assembly (trinity.sh)
4. Add the sample-specific prefix and merged all the data (prefix.combine.sh)
5. Run transDecoder and hmmer for finding ORFs and proteins domains (xxx).
6. remove the dupulicated transcripts
7. Run BUSCO

### Screening differentially expressed transcripts between different mating conditions 
1. run Salmon for quantification
2. run DEseq2 to calculate LFCs
3. Select the trascripts which increased the expression when they had a mate(s) (xxx)
