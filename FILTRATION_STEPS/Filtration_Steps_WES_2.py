#!/usr/bin/env python
# coding: utf-8

# # $\color {#042823} {\text {WES FILTER VARIANT EFFECT PREDICTOR RESULTS}}$

# Ariadna Colmenero Cobo de Guzmán | TFM |
# *This code is used prior to perform the Oncoprint. So, from the results for the VEP predictor, we get a table of potential driver mutations to analyze.*

# ## $\color {#25756A} {\text {FIRTS CONSIDERATIONS}}$

# In this section, we are working with the results provided by **SNPEff/SNPSift**. Thus, our main objective is to filter the **16748** results obtained and, in turn, to establish which mutations are **potential drivers**. This will serve as an imput for the next steps of the work: oncoprinting, machine learning or prediction using LymphGene2.

# ## $\color {#25756A} {\text {1. IMPORT LIBRARIES AND .csv}}$

# First, we must import **numpy** and **pandas** in order to use their implementations in this code. In addition, we assign them an alias that will be used to call the different functions. Once we have read our .csv (with the results obtained with VEP (*Ensembl Variant Effect Predictor* well described), we can view it (df). With this, we must consider that the names of the columns (if they are composed of more than one word) must be separated by _. 

# In[59]:


import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns


# In[60]:


df = pd.read_csv('final.csv', low_memory=False) # Read our csv of interest. It is important to have all the columns well defined. 
df # Take a look at our df 


# In[61]:


pd.options.display.max_columns = None # We use this line of code to be able to see all the columns while working in Jupiter.


# Thus, we can see how we work with a data frame containing 16748 rows and 48 columns.

# ## $\color {#25756A} {\text {2. APPLICATION OF THE DIFFERENT FILTERS}}$

# In this section, the filters to be applied to the selected data are defined. With this, consider that we should not delete the df that is created in each step. In this way, we can export it, see what results we are obtaining and assess how to proceed according to the established criteria. 

# The first important thing and following the analysis performed by *Karube et. al*, we are assigniing to all those **TRUNCATING** mutations, that they are **DRIVER** mutations. With that, we are selecting them and excluding them from further filters.

# In[62]:


truncating = ["stop_gained", "frameshift_variant", "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&start_lost", "frameshift_variant&stop_gained", "frameshift_variant&stop_lost", "frameshift_variant&stop_lost", "plice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant", "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",  "splice_acceptor_variant&intron_variant", "splice_acceptor_variant&splice_region_variant&intron_variant", "splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant", "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant", "splice_donor_variant&intron_variant", "splice_donor_variant&intron_variant"]
df_truncating = df[df.ANNOTATION.str.contains('|'.join(truncating))]


# In[63]:


df_truncating.to_excel('TRUNCATING_MUTATIONS_WES_KIEL.xlsx')


# Thus, we see that **606 mutations** can be considered as **DRIVERS**. The remaining ones will go through a **filtering** and **selection process** that follows different criteria. Thus, we can remove these results from our main dataframe. Hence, we are using the *df_filter* df.

# In[64]:


df_filter = df.loc[~((df.ANNOTATION.isin(df_truncating['ANNOTATION']))),:]
df_filter # We can see that we have removed the 606 truncating mutations.
df_filter.to_excel('df_WITHOUT_TRUNCATING.xlsx') # Let's make sure that we have done it correctly.


# ### $\color {#25756A} {\text {A) FILTER FOR THE NUMBER OF CASES}}$

# Next, and considering the **NUMBER OF CASES** column (which is only based on **position** and **Chr**), we are interested in those with a value less than or equal to 3, i.e. that the mutational variant is found in 3 or fewer cases. Thus, we eliminate the most recurrent mutations. 
# * Considerations:
#     * Despite eliminating those mutations present in more than 3 cases, we must take into account that, perhaps, some of them, although present in more cases could be of interest to us (for example, if they were part of a **hot spot**).
#     * For that reason, we are comparing those mutations which are in more than 3 cases with the list of **114 genes** used in the **Lymph2gene** predictor (*Wright et al....*) which is used for this kind of prediction.

# In[65]:


def NUMBER_OF_CASES(df): # We define a function which is selecting those mutations which are in more than 3 cases.
    number_of_cases = df[(df.NUMBER_OF_CASES_2 > 3)]
    return number_of_cases


# In[66]:


df_CASES = NUMBER_OF_CASES(df_filter) # We obtain a df with the mutations which are in more than 3 cases (in the df that we have removed the truncating ones).
df_CASES


# Next, we select all the genes which contain recurrent mutations and import the ones from LymphGen.

# In[67]:


genes1 = df_CASES['GENE_NAME'] # We set all the gene names from the df which contains the recurrent mutations
genes1 = set(genes1) # We stablish that it has to be a set.


# In[70]:


genes2 = pd.read_csv('Genes_Lymph2.csv', low_memory=False) # We import those genes used for the prediction by Lymphgen2.
genes2 = set(genes2['Gene']) # We stablish it as a set.
len(genes2)


# In[71]:


hotspots = set(genes1).intersection(genes2) # We perform the intersection between both gene sets (ours and from LymphGene2)
print(f"We have to consider that these genes could be hotspots: {hotspots}, based on the genes used in the prediction of LymphGen2.")


# In[72]:


hotspot_list = list(hotspots)
len(hotspot_list)


# Con ello, vemos que existen 6 genes que aunque quedarían eliminados por ser recurrentes, resultan interesantes en pertenecer a aquellos utilizados para la predicción.

# In[73]:


df1 = df_filter[(df_filter.GENE_NAME.isin(hotspot_list))|(df_filter.NUMBER_OF_CASES_2 <= 3)]
df1
df1.to_excel('df_CASES+HOTSPOTS.xlsx') #Saving the csv allows us to check whether it has been performed correctly.


# Thus, although we could have removed these genes because they were present in more than 3 cases, we have to mantain them as they are crutial used in the Lymph2.

# In[74]:


count_row3 = df1.shape[0] 
print(f"If select those mutations which are in less than 3 cases we get {count_row3} rows")


# ### $\color {#25756A} {\text {B) FILTER FOR DP}}$

# In[75]:


def DP_filter(df):
    DP_filter_A = df[(df.DP >= 10)]
    return DP_filter_A


# In[76]:


df2 = DP_filter(df1)
df2


# In[77]:


count_row1 = df2.shape[0] 
print(f"If we only select those rows which have a value higher than 10 in the DP column, we get {count_row1} mutations.")


# In[78]:


def DP_filter2(df):
    DP_filter_B = df[(df.DP_ALT >= 3)]
    return DP_filter_B


# In[79]:


df3 = DP_filter2(df2)
df3


# In[80]:


count_row2 = df3.shape[0] 
print(f"If further select those rows which have a value higher or equal to 3 in the DP column, we get {count_row2} rows")


# ### $\color {#25756A} {\text {C) FILTER FOR THE ANNOTATION}}$

# Although in this case, this filter is already applied in the results from which we started, it is important to take it into account in future analyses.

# In[81]:


def db_filter(df):
    db = df = df[~(df[['dbNSFP_1000Gp3_AF_simplified','dbNSFP_ExAC_AF_simplified', 'dbNSFP_gnomAD_exomes_AF_simplified']] > 0.01).any(axis=1)]
    return(db)   


# In[24]:


df4 = db_filter(df3)
df4


# ### $\color {#25756A} {\text {D) FILTER FOR THE ANNOTATION IMPACT}}$

# The first step is to remove all mutations with an annotation impact that is a **MODIFIER**. In addition, from those that are **LOW**, we will remove the rows with an annotation containing **5_prime_UTR_premature_start_codon_gain_variant**.

# In[82]:


def annotation_filter(df):
    annotation = df[(df.ANNOTATION_IMPACT != 'MODIFIER')]
    return annotation


# In[83]:


df5 = annotation_filter(df4)
df5


# In[84]:


count_row4 = df5.shape[0] 
print(f"If further remove those mutations which are MODIFIER we get {count_row4} results")


# In[85]:


impact_low = df5[(df5.ANNOTATION_IMPACT == 'LOW')]
remove = impact_low[(impact_low.ANNOTATION == '5_prime_UTR_premature_start_codon_gain_variant')]
df6 = df5.loc[~((df5.ANNOTATION.isin(remove['ANNOTATION']))&(df5.ANNOTATION_IMPACT.isin(remove['ANNOTATION_IMPACT']))),:]


# In[86]:


count_row5 = df6.shape[0] 
print(f"If further remove those LOW impact mutations which are 5_prime_UTR_premature_start_codon_gain_variant  we get {count_row5} results")


# Además, también deberíamos eliminar todas las mutaciones que presenten '3_prime_UTR_variant' excepto en aquellos casos pertenecientes a NOTCH. Aún así, cómo en éste caso, todas las '3_prime_UTR_variant
# ', son MODIFIERS y ninguna de ellas pertenece a NOTCH, ya quedan eliminadas. 
# 

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


# ## $\color {#25756A} {\text {3. MUTATIONAL PREDICTION}}$

# ### 3.1. PROVEAN for INDELS

# Next, we must apply the prediction of mutations. Thus, the potential drivers of the remaining mutations were selected based on the functional prediction established by the **OncodriveCLUST Mutation Assessor (MA)** and the **SIFT** algorithm. Thus, those that do not have a result for MA, the result of SIFT is used.

# In[90]:


import warnings
warnings.filterwarnings('ignore')


# In[91]:


df8['PREDICTION'] = df8['dbNSFP_MutationAssessor_pred'].mask(df8['dbNSFP_MutationAssessor_pred'].isna(), df8['dbNSFP_SIFT_pred'])
#df8.to_excel('FINAL_FILTERED_WES_KIEL.xlsx')


# One of the things that we have noticed is that INDELs do not have information about their annotation. Thus, we cannot conclude whether or not they are potential driver mutations. For that reason, we can predict them using **PROVEAN**. With that, we do not loose information about more than 600 mutations.

# In[35]:


sns.heatmap(df8.isnull(), cbar=False)


# Then, we have to import our df_filtered data frame with the additional information based on PROVEAN results. With that we obtain:

# In[48]:


df_FINAL = pd.read_csv('FINAL_FILTERED+PROVEAN.csv', low_memory=False) # Read our csv of interest. It is important to have all the columns well defined. 
df_FINAL


# We can further see how many mutations belong to each level of the prediction:

# In[49]:


df_FINAL['PREDICTION'].value_counts() # The ones already named as DELETERIOUS or NEUTRAL come from PROVEAN (and, thus, from INDELS).


# In[38]:


my_dictionary = {'N':'Neutral',
                             'L':'Low',
                             'M':'Medium',
                             'H':'High',
                             'D': 'Deleterious',
                             'T': 'Tolerated',
                             'T,.': 'Tolerated',
                             'N,.': 'Neutral',
                             'Neutral': 'Neutral',
                             '.,T': 'Tolerated',
                             'L,.': 'Low',
                             '.,D': 'Deleterious',
                             '.,L': 'Low',
                             'D,.,D': 'Deleterious',
                             'D,.':'Deleterious',
                             'T,.,.': 'Tolerated',
                             '.,N': 'Neutral',
                             'N,.,.': 'Neutral',
                             'N,N': 'Neutral',
                             'T,T,.,.,.,T,T,T': 'Tolerated',
                             '.,M': 'Medium',
                             'M,.': 'Medium',
                             'D,D,.,D': 'Deleterious',
                             'L,L,L,L': 'Low',
                             'L,L': 'Low',           
                             'D,D': 'Deleterious',
                             '.,D,D,D': 'Deleterious',
                             'D,.,D,D': 'Deleterious',           
                             'T,.,.,.': 'Tolerated'}


# However, although these predictions can be entered into the dictionary through the results obtained, there are some that must be evaluated manually. Thus, if we enter our prediction in the following code, it will tell us whether we really need to look for more information about it (for example in *Cosmic*). For instance, for the following ones:
# * D,.,T,T             
# * D,T                  
# * T,D                  
# * .,T,D                
# * D,T,.                
# * D,D,T  

# In[77]:





# In[39]:


try:  
  print (my_dictionary["H"])  
except:  
  print ("key error: This prediction is not in the dictionary! Please, further evaluate this mutation.") 


# In[40]:


gene_evaluate = df8[df8['PREDICTION'] == "D,D,T"]['GENE_NAME']
print(f"You have to further evaluate the following mutation to get its prediction (for instance on COSMIC): {gene_evaluate}")


# In[50]:


df_FINAL['PREDICTION'] = df_FINAL['PREDICTION'].map({'N':'Neutral',
                             'L':'Low',
                             'M':'Medium',
                             'H':'High',
                             'D': 'Deleterious',
                             'T': 'Tolerated',
                             'T,.': 'Tolerated',
                             'N,.': 'Neutral',
                             '.,T': 'Tolerated',
                             'L,.': 'Low',
                             '.,D': 'Deleterious',
                             '.,L': 'Low',
                             'D,.,D': 'Deleterious',
                             'D,.':'Deleterious',
                             'T,.,.': 'Tolerated',
                             '.,N': 'Neutral',
                             'N,.,.': 'Neutral',
                             'N,N': 'Neutral',
                             'T,T,.,.,.,T,T,T': 'Tolerated',
                             '.,M': 'Medium',
                             'M,.': 'Medium',
                             'D,D,.,D': 'Deleterious',
                             'L,L,L,L': 'Low',
                             'L,L': 'Low',           
                             'D,D': 'Deleterious',
                             '.,D,D,D': 'Deleterious',
                             'D,.,D,D': 'Deleterious',           
                             'T,.,.,.': 'Tolerated'},                                  
                             na_action=None)
#df_FINAL.to_excel('final_oncoprint.xlsx')


# In[51]:


df_FINAL['PREDICTION'].value_counts()


# In[53]:


prediction_drivers = list(df_FINAL[(df_FINAL['PREDICTION'] == "High") | (df_FINAL['PREDICTION'] == "Medium")|(df_FINAL['PREDICTION'] == "Deleterious")]['GENE_NAME'])
num_drivers = len(prediction_drivers)
print(f"The following genes are predicted as DRIVERS: {prediction_drivers}. We can see, we have {num_drivers}. ")


# Now, we have to select the driver mutations, and so, we are going to append a new column for them.

# In[55]:


df_FINAL_DRIVERS = df_FINAL[(df_FINAL['PREDICTION'] == "High") | (df_FINAL['PREDICTION'] == "Medium")|(df_FINAL['PREDICTION'] == "Deleterious")]
df_FINAL_DRIVERS


# In[57]:


FINAL_DRIVER_MUTATIONS = df_FINAL_DRIVERS.append(df_truncating, ignore_index=True)
len(FINAL_DRIVER_MUTATIONS)
#FINAL_DRIVER_MUTATIONS.to_excel('FINAL_DRIVER_WES_KIEL.xlsx')


# In[58]:


sns.heatmap(FINAL_DRIVER_MUTATIONS.isnull(), cbar=False)


# ## Exploratory analysis

# A parte de explorar cuáles de ellos son mutational potential si hacemos una tabla con los más recurrentes, más allá de su predicción, encontramos:

# In[48]:


genes = df8['GENE_NAME'].value_counts()
df_genes = pd.DataFrame(genes)
genes_1 = df_genes[df_genes.GENE_NAME > 3] # Hay 57 genes que se encuentran mutados en más de 3 casos.
genes_1


# From the selected mutations, we have to look at which of them are present in the list of **96 genes** used for Machine Learning, so we create the imput table for the following analysis.
