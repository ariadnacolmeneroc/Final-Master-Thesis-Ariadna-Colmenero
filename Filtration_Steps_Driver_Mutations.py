#!/usr/bin/env python
# coding: utf-8

# # $\color {#042823} {\text {WES FILTER VARIANT EFFECT PREDICTOR RESULTS}}$

# Ariadna Colmenero Cobo de Guzmán | TFM |
# *This code is used prior to perform the Oncoprint. So, from the results for the SNPEff/SNPSift procedure.

# ## $\color {#25756A} {\text {0. FIRTS CONSIDERATIONS}}$

# In this section, we are working with the results provided by **SNPEff/SNPSift**. Thus, our main objective is to filter the results obtained and, in turn, to establish which mutations are **potential drivers**. This will serve as an imput for the next steps of the work: **oncoprinting** or using **LymphGene2**. In fact, this particular case of code serves for all 23 samples of the High-grades analysed.

# ## $\color {#25756A} {\text {1. IMPORT LIBRARIES AND .csv}}$

# First, we must import **numpy** and **pandas** in order to use their implementations in this code. In addition, we assign them an alias that will be used to call the different functions. Once we have read our CSV, we can view it (df). With this, we must consider that the names of the columns (if they are composed of more than one word) must be separated by _. 

# In[1]:


import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns


# In[4]:


df = pd.read_csv('final_ann_HG_x17_081722.csv', low_memory=False) # Read our csv of interest. It is important to have all the columns well defined. 
df # Take a look at our df 


# In[5]:


pd.options.display.max_columns = None # We use this line of code to be able to see all the columns while working in Jupiter.


# From this point on, it is very important to export the different data frames that are created and to look carefully at which mutations are being eliminated, which groups they belong to, etc. In this way, filters can be determined on the basis of the biological significance of the results.

# ## $\color {#25756A} {\text {2. APPLICATION OF THE DIFFERENT FILTERS}}$

# In this section, the filters to be applied to the selected data are defined. With this, consider that we should not delete the df that is created in each step. In this way, we can export it, see what results we are obtaining and assess how to proceed according to the established criteria. 

# The first important thing and following the analysis performed by (*Karube et al., 2018*) we are assigniing to all those **TRUNCATING** mutations, that they are **DRIVER** mutations. With that, we are selecting them and excluding them from further filters.

# In[6]:


truncating = ["stop_gained", "disruptive_inframe_deletion", "frameshift_variant", "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&start_lost", "frameshift_variant&stop_gained", "frameshift_variant&stop_lost", "frameshift_variant&stop_lost", "splice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant", "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",  "splice_acceptor_variant&intron_variant", "splice_acceptor_variant&splice_region_variant&intron_variant", "splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant", "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant", "splice_donor_variant&intron_variant", "splice_donor_variant&intron_variant", "conservative_inframe_deletion"]


# In[7]:


df_truncating = df[df.Annotation.str.contains('|'.join(truncating))]
df_truncating


# Para cada una de éstas alteraciones, indicaremos en una nueva columna que sean TRUNCATING, así podran ser consideradas, tras su filtraje, cómo drivers directamente sin necesidad de tener que consultar ninguna base de datos.

# In[8]:


df["PREDICTION"] = df.apply(lambda x: "TRUNCATION" if df_truncating["Annotation"].isin(x).any() else " ",axis=1)
print(df)


# It can be seen that there are 1381 mutations which can be considered as **TRUNCATING** and are thus directly selected as **POTENTIAL DRIVER MUTATIONS**.

# In[19]:


df.to_excel('TRUNCATING_MUTATIONS_HG.xlsx') # Exporing the truncating mutations as an XLSX file.


# ### $\color {#25756A} {\text {A) FILTER FOR THE NUMBER OF CASES}}$

# Next, and considering the **NUMBER OF CASES** column (which is only based on **position** and **Chr**), we are interested in those with a value less than or equal to 3, i.e. that the mutational variant is found in 3 or fewer cases. Thus, we eliminate the most recurrent mutations. 
# * Considerations:
#     * Despite eliminating those mutations present in more than 3 cases, we must take into account that, perhaps, some of them, although present in more cases could be of interest to us (for example, if they were part of a **hot spot**).
#     * For that reason, we are comparing those mutations which are in more than 3 cases with the list of **114 genes** used in the **Lymph2gene** predictor (*Wright et al., 2020.*) which is used for this kind of prediction.

# In[9]:


def NUMBER_OF_CASES(df): # We define a function which is selecting those mutations which are in more than 3 cases.
    number_of_cases = df[(df.Number_Cases > 3)]
    return number_of_cases


# In[10]:


df_CASES = NUMBER_OF_CASES(df) # We obtain a df with the mutations which are in more than 3 cases (in the df that we have removed the truncating ones).
df_CASES


# Next, we select all the genes which contain recurrent mutations and import the 114 ones from **LymphGen2**.

# In[11]:


genes1 = df_CASES['Gene_name'] # We set all the gene names from the df which contains the recurrent mutations
genes1 = set(genes1) # We stablish that it has to be a set.


# In[12]:


genes2 = pd.read_csv('Genes_Lymph2.csv', low_memory=False) # We import those genes used for the prediction by Lymphgen2.
genes2 = set(genes2['Gene']) # We stablish it as a set.
len(genes2)


# In[13]:


hotspots = set(genes1).intersection(genes2) # We perform the intersection between both gene sets (ours and from LymphGene2)
print(f"We have to consider that the following genes could be hotspots: {hotspots}, based on the genes used in the prediction of LymphGen2.")


# In[14]:


hotspot_list = list(hotspots)
len(hotspot_list)


# With this, although there are 57 genes that would be eliminated because they are recurrent, they are interesting in that they belong to those used for prediction.

# In[15]:


df1 = df[(df.Gene_name.isin(hotspot_list))|(df.Number_Cases <= 3)] # Mutation selection based on both criteria.
df1


# In[64]:


df1.to_excel('df_CASES+LYMPH.xlsx') #Saving the csv allows us to check whether it has been performed correctly.


# In[16]:


count_row3 = df1.shape[0] 
print(f"If select those mutations which are in less than 3 cases we get {count_row3} rows")


# ### $\color {#25756A} {\text {B) FILTER FOR DP}}$

# In[16]:


def DP_filter(df):
    DP_filter_A = df[(df.totalDepth >= 10)]
    return DP_filter_A


# In[17]:


df2 = DP_filter(df1)
df2


# In[18]:


count_row1 = df2.shape[0] 
print(f"If we only select those rows which have a value higher than 10 in the DP column, we get {count_row1} mutations.")


# In[19]:


def DP_filter2(df):
    DP_filter_B = df[(df.altDepth >= 5)]
    return DP_filter_B


# In[20]:


df3 = DP_filter2(df2)
df3


# In[21]:


count_row2 = df3.shape[0] 
print(f"If further select those rows which have a value higher or equal to 3 in the DP column, we get {count_row2} rows")


# ### $\color {#25756A} {\text {D) FILTER FOR THE ANNOTATION IMPACT}}$

# The first step is to remove all mutations with an annotation impact that is a **MODIFIER**. In addition, from those that are **LOW**, we will remove the rows with an annotation containing **5_prime_UTR_premature_start_codon_gain_variant**.

# In[22]:


df_notch = df3[(df3.Annotation_impact == "MODIFIER") & (df3.Gene_name == "NOTCH1") & (df3.Annotation == "3_prime_UTR_variant")]


# In[23]:


def annotation_filter(df):
    annotation = df[(df.Annotation_impact != 'MODIFIER')]
    return annotation


# In[24]:


df4 = annotation_filter(df3)
df4


# In[25]:


df5 = df4.append(df_notch, ignore_index=True)


# In[26]:


count_row4 = df5.shape[0] 
print(f"If further remove those mutations which are MODIFIER we get {count_row4} results")


# In[27]:


impact_low = df5[(df5.Annotation_impact == 'LOW')]
remove = impact_low[(impact_low.Annotation == '5_prime_UTR_premature_start_codon_gain_variant')]
df6 = df5.loc[~((df5.Annotation.isin(remove['Annotation']))&(df5.Annotation_impact.isin(remove['Annotation_impact']))),:]


# In[28]:


count_row5 = df6.shape[0] 
print(f"If further remove those LOW impact mutations which are 5_prime_UTR_premature_start_codon_gain_variant  we get {count_row5} results")


# In[29]:


def synonymous_filter(df):
    annotation = df[(df.Annotation != 'synonymous_variant')]
    return annotation


# In[30]:


df6 = synonymous_filter(df6)
df6


# Además, también deberíamos eliminar todas las mutaciones que presenten '3_prime_UTR_variant' excepto en aquellos casos pertenecientes a NOTCH. Aún así, cómo en éste caso, todas las '3_prime_UTR_variant
# ', son MODIFIERS y ninguna de ellas pertenece a NOTCH, ya quedan eliminadas. 
# 

# ### $\color {#25756A} {\text {F) FILTER FOR QUALITY}}$

# In[31]:


def QUALITY_filter(df):
    annotation = df[(df.AS_FilterStatus == 'SITE')]
    return annotation


# In[32]:


df7 = QUALITY_filter(df6)
df7


# ### $\color {#25756A} {\text {F) FILTER FOR SNPs}}$

# In[33]:


df7['dbNSFP_ExAC_AF'] = df7['dbNSFP_ExAC_AF'].replace(np.nan, 0)
df7['dbNSFP_ExAC_AF'] 


# In[131]:


#df.drop(df.loc[df['dbNSFP_ExAC_AF'] >= 0.01].index, inplace=True)


# In[34]:


def db_filter(df):
    db = df[~(df['dbNSFP_ExAC_AF'] >= 0.01)]
    return(db)  


# In[35]:


df9 = db_filter(df7)
df9


# In[107]:


df9.to_excel('df9.xlsx')


# In[36]:


df_truncating = df9[df9.Annotation.str.contains('|'.join(truncating))]
df_prediction = df9.loc[~((df9.Annotation.isin(df_truncating['Annotation']))),:]
df_truncating


# ## $\color {#25756A} {\text {3. MUTATIONAL PREDICTION}}$

# ### 3.1. PROVEAN for INDELS

# Next, we must apply the prediction of mutations. Thus, the potential drivers of the remaining mutations were selected based on the functional prediction established by the **OncodriveCLUST Mutation Assessor (MA)** and the **SIFT** algorithm. Thus, those that do not have a result for MA, the result of SIFT is used.

# In[37]:


import warnings
warnings.filterwarnings('ignore')


# In[38]:


df_prediction['dbNSFP_MutationAssessor_pred'] = df_prediction['dbNSFP_MutationAssessor_pred'].replace(".", np.nan)
df_prediction['dbNSFP_SIFT_pred'] = df_prediction['dbNSFP_SIFT_pred'].replace(".", np.nan)
df_prediction['dbNSFP_PROVEAN_pred'] = df_prediction['dbNSFP_PROVEAN_pred'].replace(".", np.nan)
df_prediction['dbNSFP_Polyphen2_HVAR_pred'] = df_prediction['dbNSFP_Polyphen2_HVAR_pred'].replace(".", np.nan)


# In[39]:


df_prediction['PREDICTION']= df_prediction['dbNSFP_MutationAssessor_pred'].mask(df_prediction['dbNSFP_MutationAssessor_pred'].isna(), df_prediction['dbNSFP_SIFT_pred'])
df_prediction


# In[40]:


df_prediction.PREDICTION.fillna(df_prediction.dbNSFP_PROVEAN_pred, inplace = True)


# In[41]:


df_prediction.PREDICTION.fillna(df_prediction.dbNSFP_Polyphen2_HVAR_pred, inplace = True)


# In[42]:


df_prediction


# In[101]:


df_prediction.to_excel('FINAL_FILTERED_HG.xlsx')


# In[43]:


df_HG_prediction = df_prediction.append(df_truncating, ignore_index=True)
df_HG_prediction
#df_HG_prediction.to_excel('df_HG_prediction.xlsx')


# One of the things that we have noticed is that INDELs do not have information about their annotation. Thus, we cannot conclude whether or not they are potential driver mutations. For that reason, we can predict them using **PROVEAN**. With that, we do not loose information about more than 600 mutations.

# Then, we have to import our df_filtered data frame with the additional information based on PROVEAN results. With that we obtain:

# We can further see how many mutations belong to each level of the prediction:

# In[44]:


df_HG_prediction['PREDICTION'].value_counts() # The ones already named as DELETERIOUS or NEUTRAL come from PROVEAN (and, thus, from INDELS).


# In[45]:


my_dictionary = {'N':'Neutral',
                 'T': 'Tolerated',
                 'M': 'Medium',
                 'L': 'Low',
                 'D': 'Deleterious',
                 'H': 'High',
                 'B': 'Benign',
                 'P': 'Probably deleterious'}


# However, although these predictions can be entered into the dictionary through the results obtained, there are some that must be evaluated manually. Thus, if we enter our prediction in the following code, it will tell us whether we really need to look for more information about it (for example in *Cosmic*). For instance, for the following ones:
# * D,.,T,T             
# * D,T                  
# * T,D                  
# * .,T,D                
# * D,T,.                
# * D,D,T  

# In[107]:


try:  
  print (my_dictionary["H"])  
except:  
  print ("key error: This prediction is not in the dictionary! Please, further evaluate this mutation.") 


# In[46]:


gene_evaluate = df8[df8['PREDICTION'] == "D,D,T"]['GENE_NAME']
print(f"You have to further evaluate the following mutation to get its prediction (for instance on COSMIC): {gene_evaluate}")


# In[47]:


df_HG_prediction['PREDICTION'] = df_HG_prediction['PREDICTION'].map({'N':'Neutral',
                 'T': 'Tolerated',
                 'M': 'Medium',
                 'L': 'Low',
                 'D': 'Deleterious',
                 'H': 'High',
                 'B': 'Benign',
                 'P': 'Probably deleterious',
                 'TRUNCATION': 'Truncation',                                                  },                                  
                             na_action=None)
#df_FINAL.to_excel('final_oncoprint.xlsx')


# In[48]:


df_HG_prediction['PREDICTION'].value_counts()


# lOOK AT THEM IN COSMIC:

# In[49]:


df_HG_prediction[(df_HG_prediction['PREDICTION'] == "Probably deleterious")]


# In[50]:


prediction_drivers = list(df_HG_prediction[(df_HG_prediction['PREDICTION'] == "High") | (df_HG_prediction['PREDICTION'] == "Medium")|(df_HG_prediction['PREDICTION'] == "Deleterious")|(df_HG_prediction['PREDICTION'] == "Truncation")]['Gene_name'])
num_drivers = len(prediction_drivers)
print(f"The following genes are predicted as DRIVERS: {prediction_drivers}. We can see, we have {num_drivers}. ")


# In[51]:


drivers = ["High", "Medium", "Deleterious", "Truncation"]

df_HG_prediction['PREDICTION'] = df_HG_prediction['PREDICTION'].replace(np.nan, "Unknown")


# In[52]:


df_HG_prediction


# In[ ]:





# In[53]:


df_drivers = df_HG_prediction[df_HG_prediction.PREDICTION.str.contains('|'.join(drivers))]


# In[54]:


df_HG_prediction["DRIVERS"] = df_HG_prediction.apply(lambda x: "DRIVER" if df_drivers["PREDICTION"].isin(x).any() else " ",axis=1)
print(df_HG_prediction)


# In[65]:


df_HG_prediction.to_excel('FINAL_POTENTIALDRIVER_MUTATIONS_HG_x17_ARI_81722.xlsx')


# Now, we have to select the driver mutations, and so, we are going to append a new column for them.

# In[55]:


df_FINAL_DRIVERS = df_HG_prediction[(df_HG_prediction['PREDICTION'] == "High") | (df_HG_prediction['PREDICTION'] == "Medium")|(df_HG_prediction['PREDICTION'] == "Deleterious")|(df_HG_prediction['PREDICTION'] == "Truncation")]
df_FINAL_DRIVERS


# In[52]:


df_FINAL_DRIVERS.to_excel('FINAL_POTENTIALDRIVER_MUTATIONS_HG_X18_081722.xlsx')

