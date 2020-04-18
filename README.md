# HDAC6-score
Network-based assessment of HDAC6 activity is highly predictive of pre-clinical and clinical response to the HDAC6 inhibitor

## Transcriptomics data collection and processing
We made a comprehensive collection of breast cancer transcriptomics data of cell lines, mouse models and primary patients from five sources. We introduced normalization (if not done yet), quality assessment and removed the non-informative genes and outlier samples, if any. The expression set (eset) file was created for each of them, and were used for subsequent analysis.
* See scripts: 01.Transcriptomics_data_collection_and_processing/esetGeneration_\*.R

## RNAseq analysis of the clinical trial
We employed the RNA-SEQ/RiboZero strategy to profile the gene expression patterns of 11 FFPE biopsy samples from 10 ductal breast cancer patients who were treated with the HDAC6 inhibitor ricolinostat. We first trimmed the adpters and low-quality sequences from the raw FastQ data. We then used Salmon to perform the trancriptome quantification analysis. We normalized the libraries by calculating the CPM, and evaluated the quality of the data. Batch effect was observed and removed. Expressoin set file was generated for subsequent analysis.
* See scripts: 02.RNAseq_analysis_of_the_clinical_trial/01_quantBySalmon.sh and 02_esetGeneration.R

## Integration of human breast cancer transcriptomics data from different sources
In order to evalue the performance of HDAC6 regulon, we merged the gene expression profiles of four human breast cancer datasets, including TCGA, METABRIC, IBC, and RNA-seq data from our clinical trial, by keeping the genes shared across all four datasets. Batch effects were removed by limma.
* See script: 03.Integration_of_human_breast_cancer_transcriptomics_data_from_different_sources.R

## Next-generation HDAC6 breast cancer regulon inference
To reconstruct HDAC6 regulon with high accuracy and applicability across all breast cancer subtypes, we used SJARACNe to generate the subtype-specific regulons of 9 breast cancer subtypes with over 3,000 primary breast cancer samples from TCGA and METABRIC. We finally got a new HDAC6 breast cancer regulon with 294 targets.
* See scripts: 04.Next_generation_HDAC6_breast_cancer_regulon_inference/\*R and \*.sh
* See the new HDAC6 breast cancer regulon: 04.Next_generation_HDAC6_breast_cancer_regulon_inference/HDAC6_Breast_Cancer_Regulon.txt

## HDAC6 cancer type–specific regulon inference
We used the SJARACNe to generate the cancer type-specific regulons of 32 human primary cancer types from TCGA dataset.
* See scripts: 05.HDAC6_cancer_type_specific_regulon_inference/\*R and \*.sh

## HDAC6 score inference by NetBID
We used the gene level expression data to calculate the HDAC6 score with the function of "cal.Activity" from NetBID. The HDAC6 score was calculated through all datasets, including cell lines, mouse models, primary cancer patiens (TCGA and METABRIC), inflamatory breast cancers, clinical trial samples and the integaratd human breast cancer transcriptomics data.
* See scripts: 06.HDAC6_score_inference_by_NetBID/\*.R

## HDAC6 score inference by VIPER
We evaluated the HDAC6 relative protein activity based on the new HDAC6 breast cancer regulon with the VIPER algorithm, a powerful tool which has been aproved by the NYS Department of Health CLIA/CLEP Validatioin Unit. The VIPER-based HDAC6 score ROC analysis was performed by DarwinHealth Inc.

## Correlative analysis of HDAC6 score in breast cancer cell lines
We used ggplot2 to draw the correlaton scatter plot and do the curve fitting between HDAC6 and IC50 among cell lines.
* See script: 08.Correlative_analysis_of_HDAC6_score_in_breast_cancer_cell_lines.R

## HDAC6 score distribution across contexts
We used ggplot2 to visulize the distribution of HDAC6 score in each context. The two-tailed p-values were calculated by limma.
* See scripts: 09.HDAC6_score_distribution_accross_contexts/\*.R

## Overlaps of cancer type–specific HDAC6 regulons
The p-value of Fisher's Exact Test and overlap statistics were calculated and used to scale the similarities of HDAC6 regulon among different cancer types. The scatter plot was made by ggplots.
* See script: 10.Overlaps_of_cancer_type_specific_HDAC6_regulons.R

## Progression free survival (PFS) analysis against HDAC6 score in the clinical trial
We used the lubridate package to calculate the PFS from the dates of trial assignment and disease progression/last follow-up. The PFS analysis of the clinical trial against HDAC6 score, including the COX model fitting and  Kaplan–Meier plot, was performed by using the R package survminer. We used -0.36 as the cutoff of HDAC6 score from the ROC analysis to define HDAC6 high and low patients.
* See script: 11.PFS_analysis_against_HDAC6_score_in_the_clinical_trial.R

## Receiver operating characteristic (ROC) curve analysis of HDAC6 score in the clinical trial
We introduced the ROC curve analysis to evaluate the performance of HDAC6 score in predicting clinical response of patients. It was performed by using the R package pROC. 
* See script: 12.ROC_curve_analysis_of_HDAC6_score_in_the_clinical_trial.R

## References
1.	Tallarida RJ. Quantitative methods for assessing drug synergism. Genes Cancer 2011;2(11):1003-8 doi 10.1177/1947601912440575.
2.	Chou TC. Drug combination studies and their synergy quantification using the Chou-Talalay method. Cancer Res 2010;70(2):440-6 doi 10.1158/0008-5472.CAN-09-1947.
3.	Cancer Cell Line Encyclopedia C, Genomics of Drug Sensitivity in Cancer C. Pharmacogenomic agreement between two cancer cell line data sets. Nature 2015;528(7580):84-7 doi 10.1038/nature15736.
4.	Barretina J, Caponigro G, Stransky N, Venkatesan K, Margolin AA, Kim S, et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 2012;483(7391):603-7 doi 10.1038/nature11003.
5.	Meyers RM, Bryan JG, McFarland JM, Weir BA, Sizemore AE, Xu H, et al. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nat Genet 2017;49(12):1779-84 doi 10.1038/ng.3984.
6.	Du X, Wen J, Wang Y, Karmaus PWF, Khatamian A, Tan H, et al. Hippo/Mst signalling couples metabolic state and immune function of CD8α+ dendritic cells. Nature 2018;558(7708):141-5 doi 10.1038/s41586-018-0177-0.
7.	Dong X, Wang X, Ding L, Liu J, Yu J. NetBID2: https://jyyulab.github.io/NetBID. 2020.
8.	Pfefferle AD, Herschkowitz JI, Usary J, Harrell JC, Spike BT, Adams JR, et al. Transcriptomic classification of genetically engineered mouse models of breast cancer identifies human subtype counterparts. Genome Biol 2013;14(11):R125 doi 10.1186/gb-2013-14-11-r125.
9.	OmicSoft C. OmicSoft® ArraySuite® software. omicsoftArraySuite 2015;8.
10.	Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature 2012;486(7403):346-52 doi 10.1038/nature10983.
11.	Bekhouche I, Finetti P, Adelaide J, Ferrari A, Tarpin C, Charafe-Jauffret E, et al. High-resolution comparative genomic hybridization of inflammatory breast cancer and identification of candidate genes. PLoS One 2011;6(2):e16950 doi 10.1371/journal.pone.0016950.
12.	Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. 2011 2011;17(1):3 doi 10.14806/ej.17.1.200.
13.	Andrews S. FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/. 2010.
14.	Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods 2017;14(4):417-9 doi 10.1038/nmeth.4197.
15.	Frankish A, Diekhans M, Ferreira AM, Johnson R, Jungreis I, Loveland J, et al. GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Res 2019;47(D1):D766-D73 doi 10.1093/nar/gky955.
16.	Anders S, Huber W. Differential expression analysis for sequence count data. Genome Biol 2010;11(10):R106 doi 10.1186/gb-2010-11-10-r106.
17.	Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res 2015;43(7):e47 doi 10.1093/nar/gkv007.
18.	Prat A, Parker JS, Karginova O, Fan C, Livasy C, Herschkowitz JI, et al. Phenotypic and molecular characterization of the claudin-low intrinsic subtype of breast cancer. Breast Cancer Res 2010;12(5):R68 doi 10.1186/bcr2635.
19.	Khatamian A, Paull EO, Califano A, Yu J. SJARACNe: a scalable software tool for gene network reverse engineering from big data. Bioinformatics 2019;35(12):2165-6 doi 10.1093/bioinformatics/bty907.
20.	Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 2005;102(43):15545-50 doi 10.1073/pnas.0506580102.
21.	Alvarez MJ, Shen Y, Giorgi FM, Lachmann A, Ding BB, Ye BH, et al. Functional characterization of somatic mutations in cancer using network-based inference of protein activity. Nat Genet 2016;48(8):838-47 doi 10.1038/ng.3593.
22.	Love M, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 2014;15(12):550.
23.	Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York; 2016.
24.	Grolemund G, Wickham H. Dates and Times Made Easy with lubridate. J Stat Softw 2011;40(3):1-25.
25.	survminer: https://rpkgs.datanovia.com/survminer/index.html.
26.	Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 2011;12:77 doi 10.1186/1471-2105-12-77.
