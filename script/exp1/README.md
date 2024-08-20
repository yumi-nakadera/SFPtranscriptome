## Aim
By comparing the transcriptome of AG and PG, we screen PG specific transcripts

## Data
The data of transcriptomes are obtained from https://www.genoscope.cns.fr/sadc/projet_BYD/ 

## Sample names
BYD_AIOSRB_7_1_HJLK5BBXX.12BA019_noribo_clean.fastq	7_prostate gland female
BYD_AIOSRB_7_2_HJLK5BBXX.12BA019_noribo_clean.fastq	7_prostate gland female
BYD_ALOSRB_7_1_HJLK5BBXX.12BA021_noribo_clean.fastq	9_prostate gland male
BYD_ALOSRB_7_2_HJLK5BBXX.12BA021_noribo_clean.fastq	9_prostate gland male
BYD_AMOSRB_7_1_HJLK5BBXX.12BA022_noribo_clean.fastq	10_albumen gland male
BYD_AMOSRB_7_2_HJLK5BBXX.12BA022_noribo_clean.fastq	10_albumen gland male
BYD_AKOSRB_7_2_HJLK5BBXX.12BA020_noribo_clean.fastq	8_albumen gland female
BYD_AKOSRB_7_1_HJLK5BBXX.12BA020_noribo_clean.fastq	8_albumen gland female

We also used the referenec transcriptome from exp2. 
centroids.fasta

## File name change
Since these files have .gz but not compressed, I changed their names by removing .gz.
```
find . -depth -name "*.gz" -exec sh -c 'f="{}"; mv -- "$f" "${f%.gz}"' \;
```
## Methods
1. Align using HISAT2 (hisat.sh)
2. convert and sort using samtool (samtools.sh)
3. Run QualiMap (qualimap.sh)
4. Calculate LFCs using DEseq2 (DE_PGspecific.R)

## Fig
Fig.2.R
