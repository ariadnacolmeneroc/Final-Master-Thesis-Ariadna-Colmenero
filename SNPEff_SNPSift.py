# -*- coding: utf-8 -*-
"""
SNPEff and SNPSift code
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
"""

############################### SNPEff and SNPSift code #####################################

"""
The follwoing lines of code have been executed for the step: Genetic variant nnotation of the 
filtered Variant Calling File (VCF) resulting from the script VARIANT_CALLER_MUTECT2.py. The 
functional annotation and prediction is performed through SNPEff and SNPSift v.5.1. In this way, 
another VCF is obtained with the corresponding annotation.

IMPORTANT CONSIDERATION: As this code had to be run from a local machine without the use of the
                         Starlife server because it was down, a code that can be launched directly 
                         in Ubuntu on Windows (18.04.5) is used.
"""

# First the code to launch one of the jobs is shown. That is, to get the annotation of one of the filtered VCF:
    
java -jar SnpSift.jar dbnsfp -v -db db/dbNSFP4.1a.txt.gz D4669_filt.vcf -m -a | 
java -Xmx12g -jar SnpSift.jar annotate -info 'CLNDISDB,CLNDISDBINCL,CLNDN,CLNSIG,CLNREVSTAT' CLINVAR/clinvar.vcf.gz /dev/stdin -a |
java -Xmx12g -jar snpEff.jar -canon GRCh37.75 > D4669_filt_ann_c.vcf

# However, to save computation time, we launch all the jobs in parallel, so we get the resulting VCF for all the 
# files in the directory.
