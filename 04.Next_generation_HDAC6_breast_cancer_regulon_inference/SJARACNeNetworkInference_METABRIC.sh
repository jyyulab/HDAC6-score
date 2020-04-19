#!usr/bin/bash 

# Coded by Qingfei Pan (Qingfei.Pan@stjude.org)

# python: 3.6.1
# scMINER: https://github.com/jyyulab/scMINER

# This script runs on files in the current working directory


for sample in HRpositive_HP HRpositive_LP HER2positive tripleNegative LumA LumB Her2 Basal Claudin_Low Normal
do
	python scMINER SJARACNE $sample ./METABRIC/$sample/log2Intensity/GeneExpression.txt ./METABRIC/$sample/log2Intensity/sig/HubGene_sig.txt ./METABRIC/$sample/log2Intensity/sig/
	python scMINER SJARACNE $sample ./METABRIC/$sample/log2Intensity/GeneExpression.txt ./METABRIC/$sample/log2Intensity/tf/HubGene_tf.txt ./METABRIC/$sample/log2Intensity/tf/
done
