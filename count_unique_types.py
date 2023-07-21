#!/usr/bin/env python
# coding: utf-8

# In[32]:


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
    #my_GT = my_GT.to_n_alt(fill=-1)
    
    return my_GT


# In[34]:


#input_path = "/mnt/stsi/stsi0/raqueld/VMV_VCF_Extractions/chr22/HRC.r1-1.EGA.GRCh37.chr22.haplotypes.17274081-17382360.vcf.VMV1"
input_path = sys.argv[1]
#tabix required by scikit-allele, otherwise loading may be slower
tabix_path='/opt/applications/samtools/1.9/gnu/bin/tabix'


# In[54]:


my_GT = extract_genotypes_allel(input_path)
my_DS = my_GT.to_n_alt(fill=-1).T
my_GT = my_GT.transpose((1,0,2))
_, unique_indexes_DS = np.unique(my_DS, axis=0, return_index=True)
_, unique_indexes_GT = np.unique(my_GT, axis=0, return_index=True)
my_GT = np.reshape(my_GT, [my_GT.shape[0], my_GT.shape[1]*2, -1])
my_HAP = np.vstack((my_GT[:,0::2],my_GT[:,1::2]))
_, unique_indexes_HAP = np.unique(my_HAP, axis=0, return_index=True)


# In[63]:


NDS = len(unique_indexes_DS)
NDIP = len(unique_indexes_GT)
NHAP = len(unique_indexes_HAP)
NHET = np.count_nonzero(my_DS == 1)
result_array = np.array([NDS,NDIP,NHAP,NHET])


# In[64]:


np.savetxt(input_path+".NDS_NDIP_NHAP_NHET", [result_array], delimiter="\t")


# In[65]:


print("Outputs at:")
print(input_path+".NDS_NDIP_NHAP_NHET")

