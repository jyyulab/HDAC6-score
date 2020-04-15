#!usr/bin/bash 

# Coded by Qingfei Pan (Qingfei.Pan@stjude.org)

# FastQC: 0.11.5
# cutadapt: 1.13
# Salmon: 0.9.1
# Reference genome: hg38 (GENCODE, v28)

# This script runs on files in the current working directory

hg38=../../Database/References/hg38/Gencode/Salmon/index
tr2gene=../../Database/References/hg38/Gencode/gencode.v28.transcript2geneSymbol.txt

# 1 QC by FastQC
for sample in pt03 pt04 pt05 pt06 pt07 pt08 pt09 pt13 pt14 pt15rep1 pt15rep2
do
	fastqc ./$sample/$sample.L001_R1_001.fastq.gz -o ./$sample
done

# 2 adapter trimming and low-quality sequence removal
for sample in pt03 pt04 pt05 pt06 pt07 pt08 pt09 pt13 pt14 pt15rep1 pt15rep2
do
	cutadapt -a AGATCGGAAGAG --trim-n --max-n=0.5 --quality-base=33 -q 30 -m 30 -o ./$sample/$sample.cutadapt.fq.gz ./$sample/$sample.L001_R1_001.fastq.gz
	salmon quant -i $hg38 -l A -p 8 -g $tr2gene -o ./$sample/salmon -r ./$sample/$sample.cutadapt.fq.gz
done