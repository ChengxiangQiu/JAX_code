
#################################
### Section - 7, Birth series ###
#################################

#################################################################################################################
### To systematically identify which cell types exhibit abrupt transcriptional changes before vs. after birth ###
#################################################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import gc
import time

WORK_PATH = "./"

adata_1 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_1.h5ad'))
adata_2 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_2.h5ad'))
adata_3 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_3.h5ad'))
adata_4 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_4.h5ad'))

adata = adata_1.concatenate(adata_2, adata_3, adata_4)
del adata_1, adata_2, adata_3, adata_4
gc.collect()

day_include = ["E16.0", "E16.25", "E16.5", "E16.75", "E17.0", "E17.25", "E17.5", 
               "E17.75", "E18.0", "E18.25", "E18.5", "E18.75", "P0"]

celltype_list = {}
file = open(os.path.join(WORK_PATH, "celltype_include.txt"))
for line in file:
    l = line.rstrip().split('\t')
    celltype_list[l[0]] = l[1]
file.close()

for celltype_i in celltype_list:

    adata = adata_i[adata_i.obs["celltype_update"] == celltype_i]
    adata = adata_i[adata_i.obs["day"].isin(day_include)]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2500)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

    sc.tl.umap(adata, min_dist=0.3, n_components=3)
    adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
    adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
    adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

    sc.tl.umap(adata, min_dist=0.3, n_components=2)
    adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
    adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])

    adata.obs.to_csv(os.path.join(WORK_PATH, '%s_adata_scale.obs.csv'%celltype_list[celltype_i]))

    X = adata.obsm['X_pca']
    print(X.shape)
    np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%celltype_list[celltype_i]), X, delimiter=",", fmt='%1.3f')

