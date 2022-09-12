# Supplemental information
Ariadna Colmenero Cobo de Guzman | MSc Omics Data Analysis

This branch contains a series of files, either outputs or inputs of different functions, so that you can see an **example** of them. This includes:
* *Example_Copy_Number_FlatList_LymphGen.txt*: Example input for Lymphgen in case you have information for the Copy Number.
* *Example_ENTREZID_CopyGeneList.txt*: Example input of the gene list for Lymphgen in case you have information for the Copy Number.
* *Example_MutationFlat_LymphGen.txt*: Example of a mutation flat needed for the prediction via LymphGen.
* *Example_PROVEAN_results.tsv*: Output retrieved from PROVEAN (when it still worked) and needed for the filtering steps applied to mutational variants.
* *Example_SAMPLE_ANNOTATION_LymphGen.txt*: Input which refers to BCL2, BCL6 translocations and CN availability for LymphGen prediction. 
* *Example_VCF.txt*: Output retrieved from Mutect2. This is the one used for the implementation of SNPEff/SNPSift.
* *Example_VCF_annotatedviaSNPEff.txt*: Information provided by SNPEff/SNPSift once the VCF has been annotated.  
* *Example_insert_size_metrics.txt*: Information concerning the insert size of each sample. This file can optionally be used to apply MultiQC.
* *Example_mark_duplicates.txt*: Information concerning the duplicates marked before applying the variant caller (Mutect2). This information can optionally serve as input to the MultiQC report.
* *Example_matrix_01_Oncoprint.xlsx*: Matrix obtained once we have applied the filtering to the mutational variants and obtained the potential driver mutations. It is used as input to perform the Oncoprint.
* *Example_output_hs_metrics.txt*: Example of a quality control metrics report created using HsCollectMetrics via Picard.
* *Example_subset_VCF_annotated.txt*: Example of subset information of the columns of interest for later application of the potential driver mutation filters. These columns are obtained from reports such as the *Example_VCF_annotatedviaSNPEff*.
* *Example_FINAL_CLASSIFICATION_LEARNER_TABLE_WITHOUT_MYC_REARRANGEMENT*:This file belongs to the dataset with which the proposed Bagged trees model has been trained. That is to say, they are the 104 samples and 84 predictors proposed in a first instance.

In addition, examples of the interactive reports created using **MultiQC** have also been uploaded.
* *multiqc_report_HGBCL_samples.html*
* *multiqc_report_WES_KIEL.html*
