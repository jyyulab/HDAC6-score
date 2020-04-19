#!usr/bin/bash 

# Coded by Qingfei Pan (Qingfei.Pan@stjude.org)

# python: 3.6.1
# scMINER: https://github.com/jyyulab/scMINER

# This script runs on files in the current working directory


for sample in ACC_T BLCA_T CESC_T CHOL_T COAD_T DLBC_T ESCA_T GBM_T HNSC_T KICH_T KIRC_T KIRP_T LAML_T LGG_T LIHC_T LUAD_T LUSC_T MESO_T OV_T PAAD_T PCPG_T PRAD_T READ_T SARC_T SKCM_T STAD_T TGCT_T THCA_T THYM_T UCEC_T UCS_T UVM_T
do
	python scMINER SJARACNE $sample ./$sample/log2FPKM_0.1/GeneExpression.txt ./$sample/log2FPKM_0.1/sig/HubGene_sig.txt ./$sample/log2FPKM_0.1/sig/
	python scMINER SJARACNE $sample ./$sample/log2FPKM_0.1/GeneExpression.txt ./$sample/log2FPKM_0.1/tf/HubGene_tf.txt ./$sample/log2FPKM_0.1/tf/
done
