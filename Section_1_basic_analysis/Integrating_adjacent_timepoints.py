
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
from annoy import AnnoyIndex

WORK_PATH = './'

file = open(os.path.join(WORK_PATH, "batch_list.txt"))
example_list = [line.rstrip() for line in file]
file.close()

for example_i in example_list:
    print(example_i)

    adata = sc.read_mtx(os.path.join(WORK_PATH, '%s.gene_count.mtx'%example_i))
    fdata = pd.read_csv(os.path.join(WORK_PATH, '%s.df_gene.csv'%example_i), index_col = 0)
    pdata = pd.read_csv(os.path.join(WORK_PATH, '%s.df_cell.csv'%example_i), index_col = 0)
    adata.obs = pdata
    adata.var = fdata

    print(adata.shape)

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

    adata.obs.to_csv(os.path.join(WORK_PATH, '%s_adata_scale.obs.csv'%example_i))

    adata.write(os.path.join(WORK_PATH, '%s_adata_scale.h5ad'%example_i), compression="gzip")

    X = adata.obsm['X_pca']
    print(X.shape)
    np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%example_i), X, delimiter=",", fmt='%1.3f')









