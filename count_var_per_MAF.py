#!/usr/bin/env python
# coding: utf-8

# In[109]:


import os
import allel
import numpy as np
import datetime
import allel
import time
import pandas as pd
import sys

# %%
#UTILITIES

def cmd_exists(cmd):
    for path in os.environ["PATH"].split(os.pathsep):
        if(os.access(os.path.join(path, cmd), os.X_OK) == True):
            return path+"/"+cmd
    return False

tabix_path = cmd_exists('tabix')
if(tabix_path == False):
    tabix_path='/opt/applications/samtools/1.9/gnu/bin/tabix'
    #tabix_path='/home/rdias/bin/tabix'
    if(os.path.isfile(tabix_path)==False):
        print("WARNING!!! Tabix not found, VCF extraction may be slow!!!")
        tabix_path = None
    else:
        print ("Using tabix at:", tabix_path)
else:
    print ("Using tabix at:", tabix_path)


# In[114]:


def extract_genotypes_allel(my_path, region=None):
    my_GT = []
    callset = allel.read_vcf(my_path, region=region, tabix=tabix_path,alt_number=2)
    my_GT = callset['calldata/GT']
    coord = [str(i)+':'+str(j)+'_'+','.join(filter(None,x))+'_'+','.join(filter(None,y)) for i,j,x,y in zip(callset['variants/CHROM'],callset['variants/POS'],callset['variants/REF'],callset['variants/ALT'])]
    sample_ids = callset['samples']
    sample_ids = [i.split('_')[0] for i in sample_ids]
    my_GT = allel.GenotypeArray(my_GT,dtype='i1')
    ac = my_GT.count_alleles()
    my_GT = my_GT.to_n_alt(fill=-1)
    results = dict(zip(sample_ids, my_GT.T))
    #REF,ALT frequecnies
    MAF = ac.to_frequencies()
    #ALT only frequecnies
    MAF = np.round(np.amin(MAF, axis=1),6)
    return results, MAF, coord


# In[115]:


#my_input = '/mnt/stsi/stsi0/raqueld/VMV_VCF_Extractions/chr22/HRC.r1-1.EGA.GRCh37.chr22.haplotypes.17274081-17382360.vcf.VMV1'
my_input = sys.argv[1]
name = os.path.basename(my_input)
my_output = sys.argv[1]+"_counts_per_MAF.csv"


# In[117]:


gt, maf, coord = extract_genotypes_allel(my_input)
cut_bins = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
maf_df = pd.DataFrame(maf, columns=['MAF'])
maf_df[name] = pd.cut(maf_df['MAF'], bins=cut_bins)
counts = maf_df[name].value_counts(normalize=True).sort_index().to_frame().transpose().round(8)
counts.to_csv(my_output, encoding='utf-8', sep='\t', index=True)
print("output at:",my_output)


# In[122]:


#get_ipython().system('jupyter nbconvert VCF_count_per_MAF.ipynb --to script')

