from cmdstanpy import model
import numpy as np
import scanpy as sc
import json
import pandas as pd
import os
from cmdstanpy import cmdstan_path, CmdStanModel
import scipy.stats as stats
import matplotlib.pyplot as plt
import random

from cmdstanpy import cmdstan_path, set_cmdstan_path


set_cmdstan_path('/Users/alexander/cmdstan')


infer=pd.read_csv("/Users/alexander/classes/Ml prob prog/project/N_infer_full_10.16.21",index_col=0)

adata=sc.read_h5ad("/Users/alexander/classes/Ml prob prog/project/N_ribas310_clones.h5")
X = adata.layers['counts']


#Joy's subsampling code

def subsample(adata,key,index,n):
    '''
    adata: AnnData object to subsample from
    key: key in adata.obs that contains metadata of interest
    index: value of key of interest, i.e. keep if adata.obs[key]==index
    n: number of cells to select from cells for which adata.obs[key]==index
    '''
    keep=[i for i in range(adata.shape[0]) if adata.obs[key][i]==index]
    keep_subsampled=np.random.choice(keep,size=n,replace=False)
    return adata[keep_subsampled,:]


#filter adata genes down to the genes that also appear in the inferCNV data
keep=[]
for i in range(adata.shape[1]):
    if adata.var_names[i] in infer.columns:
        keep.append(i)
adata=adata[:,keep]

sc.pp.highly_variable_genes(adata,n_top_genes=800,subset=True)

clonal_ad=[]
for i in np.unique(adata.obs["inferCNV_clones"]):
     clonal_ad.append(subsample(adata,"inferCNV_clones",i,300))

adata_subsampled=clonal_ad[0].concatenate(clonal_ad[1:])


infer_subsampled=np.zeros((adata_subsampled.shape))
for i in range(adata_subsampled.shape[0]):
    for j in range(adata_subsampled.shape[1]):
        barcode_ad=adata_subsampled.obs_names[i]
        barcode_infer=barcode_ad.split("-")[0]+"-"+barcode_ad.split("-")[1]+"-"+barcode_ad.split("-")[2]
        infer_subsampled[i][j]=infer.loc[barcode_infer][adata_subsampled.var_names[j]]

small_infer=infer_subsampled.tolist()

def get_sampled_infer():
    return small_infer

def get_sampled_adata():
    return adata_subsampled
