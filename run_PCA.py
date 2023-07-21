#!/usr/bin/env python
# coding: utf-8

# In[25]:


import allel
import sys
from sklearn.decomposition import PCA
import numpy as np
#create a scikit-allele VCF loading function
def extract_genotypes_allel(my_path):
    my_GT = []
    callset = allel.read_vcf(my_path, tabix=tabix_path,alt_number=1)
    my_GT = callset['calldata/GT']
    my_GT = allel.GenotypeArray(my_GT,dtype='i1')
    my_GT = my_GT.to_n_alt(fill=-1)
    
    return my_GT


# In[12]:


#input_path = "/mnt/stsi/stsi0/raqueld/VMV_VCF_Extractions/chr22/HRC.r1-1.EGA.GRCh37.chr22.haplotypes.17274081-17382360.vcf.VMV1"
input_path = sys.argv[1]
#tabix required by scikit-allele, otherwise loading may be slower
tabix_path='/opt/applications/samtools/1.9/gnu/bin/tabix'


# In[13]:


dosages = extract_genotypes_allel(input_path).T


# In[18]:


#Bayesian method to determine the dimensionality of PCA
#https://tminka.github.io/papers/pca/minka-pca.pdf
#my_pca = PCA(n_components = "mle", svd_solver ="auto")
my_pca = PCA(n_components = 100, svd_solver ="auto")
result = my_pca.fit_transform(dosages)


# In[16]:


cumulative_explained_ratio=0
n_components_needed=0
for ratio in my_pca.explained_variance_ratio_:
    cumulative_explained_ratio+=ratio
    n_components_needed+=1
    if cumulative_explained_ratio >=0.9:
        break


# In[23]:


out_file = open(input_path+".NCOMP", "wt")
n = out_file.write(str(n_components_needed)+"\n")
out_file.close()


# In[30]:


np.savetxt(input_path+".PCA_EXP_RATIO", [my_pca.explained_variance_ratio_], delimiter="\t")


# In[ ]:


print("Outputs at:")
print(input_path+".NCOMP")
print(input_path+".PCA_EXP_RATIO")

