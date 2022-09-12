# SCRIPTS
Ariadna Colmenero Cobo de Guzmán | MSc Omics Data Analysis

In this branch of the repository, you can find the most important **codes** necessary for the realisation of the project. Thus, they are divided into folders according to which part they belong to. If they are sorted in order of execution, we find:

**1)** *FROM RAW READS TO BAM*: In this folder you will find the SH, PY and TXT needed to batch the job concerning the mapping of raw reads together with the marking of PCR duplicates and the retrieve of insert size metrics.

**2)** *QC METRICS*: It contains both the commands for transforming a BED file into an interval list, how to apply HsCollectMetrics via Picard, and instructions on how to use MultiQC for quality control reporting. 

**3)** *VARIANT CALLER MUTECT2*: In this folder you will find the scripts concerning the application of the variant caller, as well as the implementation of BaseRecalibrator as ApplyBQSR.

**4)** *SNPEFF SNPSift*: The files refer to both the annotation of the VCFs resulting from the variant caller and their further processing. That is, the selection of the columns of interest necessary for subsequent filtering.

**5)** *FILTRATION STEPS*: Filters that have been applied to the mutational variants once they have been annotated. Thus, the aim of this code is to obtain the potential driver mutations for each of the samples analysed.

**6)** *LYMPHGEN*: Elaboration of the input files for this predictor. In addition to the mandatory inputs, the code necessary for the elaboration of files to further adjust the prediction is also shown.

**7)** *CLASSIFICATION LEARNER*: Code used for the elaboration of the chosen classification model. In addition, some outputs are also shown in order to better explain the procedure.

**8)** *OTHER SCRIPTS*: Scripts concerning the elaboration of figures and statistical tests. The justification for the application of each of them can be found in the code.
