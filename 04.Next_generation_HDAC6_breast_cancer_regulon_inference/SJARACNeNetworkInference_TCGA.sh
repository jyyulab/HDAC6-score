#!usr/bin/bash 

# Coded by Qingfei Pan (Qingfei.Pan@stjude.org)

# python: 3.6.1
# scMINER: https://github.com/jyyulab/scMINER

# This script runs on files in the current working directory


for sample in HRpositive HER2positive tripleNegative normalSamples LumA LumB Her2 Basal Normal
do
	python scMINER SJARACNE $sample ./TCGA/$sample/log2FPKM_0.1/GeneExpression.txt ./TCGA/$sample/log2FPKM_0.1/sig/HubGene_sig.txt ./TCGA/$sample/log2FPKM_0.1/sig/
	python scMINER SJARACNE $sample ./TCGA/$sample/log2FPKM_0.1/GeneExpression.txt ./TCGA/$sample/log2FPKM_0.1/tf/HubGene_tf.txt ./TCGA/$sample/log2FPKM_0.1/tf/
done
