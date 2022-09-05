
############################### SNPEff and SNPSift code #####################################

"""
SNPEff and SNPSift code
Ariadna Colmenero Cobo de Guzman
Master's Final Thesis
2021-2022
"""

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

import subprocess # Call python scripts from python.

for f in glob.glob("*_filt.vcf"):
      avinput = f.replace("_filt.vcf", ".avinput")
      final = f.replace("_filt.vcf", ".txt")
      bashArguments = "perl /media/ramis/Elements/SCRIPTS/annovar/convert2annovar.pl -format vcf4 "+f+" -outfile "+avinput+" -includeinfo -withzyg"
      subprocess.call(bashArguments, shell=True)
      
      bashArguments = "perl /media/ramis/Elements/SCRIPTS/annovar/table_annovar.pl "+avinput+" /media/ramis/Elements/annovar/humandb/ -buildver hg19 -out "+final+" -remove -protocol refGene,cytoBand,exac03,gnomad211_exome,gnomad211_genome,dbnsfp33a -operation g,r,f,f,f,f -nastring . -otherinfo -thread 8"
      subprocess.call(bashArguments, shell=True)
