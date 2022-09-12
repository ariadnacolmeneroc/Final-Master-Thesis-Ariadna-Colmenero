
# Consideration of filtering steps
Ariadna Colmenero Cobo de Guzman | TFM MSc Omics Data Analysis

All the necessary steps for the **filtering** of the mutational variants obtained by the variant caller (Mutect2) and annotated by SNPEff/SNPSift, can be found in this folder. Even so, and as indicated in the report, an extra file is added, with the steps or changes that must be made if instead of working with **Target Sequencing** data, you work with **WES** data. We then recommend not to apply a general filtering directly to your data, but to study it case by case and, mutation by mutation, so that it is as accurate as possible. In this way, we will also be able to take into account all relevant exceptions and examine whether the predictions of the databases make sense with each other.

Thus, this filtering has been applied to the 17 HGBCL, NOS, the 5 DLBCL samples from own data and the 25 DLBCL samples from Reddy (*Reddy et al., 2017*) and collaborators.

BIBLIOGRAPHY
Reddy, A. et al. Genetic and Functional Drivers of Diffuse Large B Cell Lymphoma. Cell 171, (2017).
