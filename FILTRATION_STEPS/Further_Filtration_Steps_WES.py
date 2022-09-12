#!/usr/bin/env python
# coding: utf-8

# # $\color {#042823} {\text {WES ADDITIONAL FILTERS AND CONSIDERATIONS}}$

# Ariadna Colmenero Cobo de Guzmán | TFM | MSc Omics Data Analysis

# In the case of the WES data, and as information on the Panel of Normals (PON) was available, an additional filtering is performed for this parameter. Furthermore, in the case of Target Sequencing, there were no Immunoglobulins (IGs) in the corresponding panel and, therefore, no potential driver mutation belonged to them. However, in this case, as we are working with whole exome, it would be necessary to eliminate them. 
# This is because they already present a high mutation rate in their heavy chain, and, therefore, we could confuse it with a hotspot and associate it with a potential driver mutation (Li et al., 2019).

# ### $\color {#25756A} {\text {E) FILTER FOR THE NORMALS}}$

# In this section we are removing those mutations which are at least in two normals. For that reason, we are only selecting those which are 0 or 1.

# In[87]:
def PON_filter(df):
    PON_filter = df[(df.PON == 0) | (df.PON == 1)]
    return PON_filter

# In[88]:
df7 = PON_filter(df6)
df7


# ### $\color {#25756A} {\text {F) REMOVAL OF IGs}}$

# Finally, we remove all those rows which have an immunoglobulin as a gene name. For that reason, we are dropping those which start with IG.

# In[89]:
discard = ["IG"]
df8 = df7[~df7.GENE_NAME.str.contains('|'.join(discard))]
df8


# Consideration:
# In addition, when this filtering was being carried out, the website https://www.jcvi.org/research/provean?jobid=1331373500324697, allowed us to enter the mutational variants online and returned an output with their prediction: whether they were deleterious, beningnesian, etc. However, after the withdrawal of this website, the PROVEAN prediction has continued to be carried out via SNPSift.

# BIBLIOGRAPHY
# Li, X., Wu, N., & Li, B. (2019). A high mutation rate of immunoglobulin heavy chain variable region gene associates with a poor survival and chemotherapy response of mantle cell lymphoma patients. Medicine, 98(22), e15811. https://doi.org/10.1097/MD.0000000000015811
