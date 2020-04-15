#BSUB -P HDAC6
#BSUB -n 1
#BSUB -M 64000
#BSUB -oo std.out -eo std.err
#BSUB -J Test
#BSUB -q priority

##module load python/3.6.1
##export PYTHON_PATH=/research/rgs01/applications/hpcf/apps/python/install/3.6.1/bin/python3
##export SJARACNE_PATH=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/SJARACNe-master
##export scMINER_PATH=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/scMINER/scMINER-master/
python3=/research/rgs01/applications/hpcf/apps/python/install/3.6.1/bin/python3
sjaracne=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/SJARACNe-master/generate_pipeline.py
scMINER=/research/projects/yu3grp/scRNASeq/yu3grp/qpan/Software/scMINER/scMINER-master/scMINER.py
dir=/research/projects/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/02_subNetworks/02_METABRIC/01_sjaracne

##python3 $scMINER SJARACNE ALL_Tumor $dir/ALL_Tumor/logIntensity/GeneExpression.txt $dir/ALL_Tumor/logIntensity/sig/HubGene_sig.txt $dir/ALL_Tumor/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE ALL_Tumor $dir/ALL_Tumor/logIntensity/GeneExpression.txt $dir/ALL_Tumor/logIntensity/tf/HubGene_tf.txt $dir/ALL_Tumor/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_HER2positive $dir/IHC_HER2positive/logIntensity/GeneExpression.txt $dir/IHC_HER2positive/logIntensity/sig/HubGene_sig.txt $dir/IHC_HER2positive/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_HER2positive $dir/IHC_HER2positive/logIntensity/GeneExpression.txt $dir/IHC_HER2positive/logIntensity/tf/HubGene_tf.txt $dir/IHC_HER2positive/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_HRpositive_HP $dir/IHC_HRpositive_HP/logIntensity/GeneExpression.txt $dir/IHC_HRpositive_HP/logIntensity/sig/HubGene_sig.txt $dir/IHC_HRpositive_HP/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_HRpositive_HP $dir/IHC_HRpositive_HP/logIntensity/GeneExpression.txt $dir/IHC_HRpositive_HP/logIntensity/tf/HubGene_tf.txt $dir/IHC_HRpositive_HP/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_HRpositive_LP $dir/IHC_HRpositive_LP/logIntensity/GeneExpression.txt $dir/IHC_HRpositive_LP/logIntensity/sig/HubGene_sig.txt $dir/IHC_HRpositive_LP/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_HRpositive_LP $dir/IHC_HRpositive_LP/logIntensity/GeneExpression.txt $dir/IHC_HRpositive_LP/logIntensity/tf/HubGene_tf.txt $dir/IHC_HRpositive_LP/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_tripleNegative $dir/IHC_tripleNegative/logIntensity/GeneExpression.txt $dir/IHC_tripleNegative/logIntensity/sig/HubGene_sig.txt $dir/IHC_tripleNegative/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE IHC_tripleNegative $dir/IHC_tripleNegative/logIntensity/GeneExpression.txt $dir/IHC_tripleNegative/logIntensity/tf/HubGene_tf.txt $dir/IHC_tripleNegative/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Basal $dir/PM50_Basal/logIntensity/GeneExpression.txt $dir/PM50_Basal/logIntensity/sig/HubGene_sig.txt $dir/PM50_Basal/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Basal $dir/PM50_Basal/logIntensity/GeneExpression.txt $dir/PM50_Basal/logIntensity/tf/HubGene_tf.txt $dir/PM50_Basal/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Claudin_Low $dir/PM50_Claudin_Low/logIntensity/GeneExpression.txt $dir/PM50_Claudin_Low/logIntensity/sig/HubGene_sig.txt $dir/PM50_Claudin_Low/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Claudin_Low $dir/PM50_Claudin_Low/logIntensity/GeneExpression.txt $dir/PM50_Claudin_Low/logIntensity/tf/HubGene_tf.txt $dir/PM50_Claudin_Low/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Her2 $dir/PM50_Her2/logIntensity/GeneExpression.txt $dir/PM50_Her2/logIntensity/sig/HubGene_sig.txt $dir/PM50_Her2/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Her2 $dir/PM50_Her2/logIntensity/GeneExpression.txt $dir/PM50_Her2/logIntensity/tf/HubGene_tf.txt $dir/PM50_Her2/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_LumA $dir/PM50_LumA/logIntensity/GeneExpression.txt $dir/PM50_LumA/logIntensity/sig/HubGene_sig.txt $dir/PM50_LumA/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_LumA $dir/PM50_LumA/logIntensity/GeneExpression.txt $dir/PM50_LumA/logIntensity/tf/HubGene_tf.txt $dir/PM50_LumA/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_LumB $dir/PM50_LumB/logIntensity/GeneExpression.txt $dir/PM50_LumB/logIntensity/sig/HubGene_sig.txt $dir/PM50_LumB/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_LumB $dir/PM50_LumB/logIntensity/GeneExpression.txt $dir/PM50_LumB/logIntensity/tf/HubGene_tf.txt $dir/PM50_LumB/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Normal $dir/PM50_Normal/logIntensity/GeneExpression.txt $dir/PM50_Normal/logIntensity/sig/HubGene_sig.txt $dir/PM50_Normal/logIntensity/sig/ --resource 4000 64000 64000 64000 --queue compbio
##python3 $scMINER SJARACNE PM50_Normal $dir/PM50_Normal/logIntensity/GeneExpression.txt $dir/PM50_Normal/logIntensity/tf/HubGene_tf.txt $dir/PM50_Normal/logIntensity/tf/ --resource 4000 64000 64000 64000 --queue compbio

