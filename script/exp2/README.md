# Exp 2
We reared the snails isolated, paired and grouped with five individuals. At the end of experiment, we collected their PG to screen differentially expressed transcripts.

## Data 

## Methods 
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
4. 
