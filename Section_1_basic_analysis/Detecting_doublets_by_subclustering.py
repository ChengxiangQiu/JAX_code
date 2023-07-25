
### We performed two rounds of clustering and used the doublet annotations to identify subclusters that are enriched in doublets

import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import os, sys

os.mkdir("doublet_cluster")

fdata = pd.read_csv("df_gene.csv", index_col = 0)

adata1 = sc.read_mtx('gene_count_1.mtx')
pdata1 = pd.read_csv('df_cell_1.csv', index_col = 0)
adata1.obs = pdata1
adata1.var = fdata

adata2 = sc.read_mtx('gene_count_2.mtx')
pdata2 = pd.read_csv('df_cell_2.csv', index_col = 0)
adata2.obs = pdata2
adata2.var = fdata

adata3 = sc.read_mtx('gene_count_3.mtx')
pdata3 = pd.read_csv('df_cell_3.csv', index_col = 0)
adata3.obs = pdata3
adata3.var = fdata

adata4 = sc.read_mtx('gene_count_4.mtx')
pdata4 = pd.read_csv('df_cell_4.csv', index_col = 0)
adata4.obs = pdata4
adata4.var = fdata

adata5 = sc.read_mtx('gene_count_5.mtx')
pdata5 = pd.read_csv('df_cell_5.csv', index_col = 0)
adata5.obs = pdata5
adata5.var = fdata

adata6 = sc.read_mtx('gene_count_6.mtx')
pdata6 = pd.read_csv('df_cell_6.csv', index_col = 0)
adata6.obs = pdata6
adata6.var = fdata

adata = adata1.concatenate(adata2, adata3, adata4, adata5, adata6)

adata_orig = adata

### remove sex genes
adata = adata_orig[:, ~adata_orig.var['chr'].isin(['chrX', 'chrY'])]
### high variable genes
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
filter_genes = list(adata.var.loc[adata.var['highly_variable'] == True, 'gene_id'])

### 
adata = adata_orig[:, adata_orig.var['gene_id'].isin(filter_genes)]
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)
sc.pp.scale(adata)
###
sc.tl.pca(adata, svd_solver='arpack', n_comps = 30)
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)
sc.tl.louvain(adata)
sc.tl.umap(adata, min_dist=0.1)

adata.obs['umap_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['umap_2'] = list(adata.obsm['X_umap'][:,1])
name = "global.csv"
adata.obs.to_csv(os.path.join('doublet_cluster', name))

obs_all = adata.obs
obs_all['louvain'].value_counts()
cluster_list = list(set(list(obs_all['louvain'])))

xx = 0
for cnt in range(len(cluster_list)):
    xx += 1
    print('Processing: ' + str(xx) + '/' + str(len(cluster_list)))
    cluster = cluster_list[cnt]
    include_cell = list(obs_all.loc[obs_all['louvain'] == cluster, 'sample'])
    adata = adata_orig[adata_orig.obs['sample'].isin(include_cell)]
    adata = adata[:, ~adata.var['chr'].isin(['chrX', 'chrY'])]
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    filter_genes = list(adata.var.loc[adata.var['highly_variable'] == True, 'gene_id'])
    adata = adata_orig[adata_orig.obs['sample'].isin(include_cell)]
    adata = adata[:, adata.var['gene_id'].isin(filter_genes)]
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack', n_comps = 30)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
    sc.tl.louvain(adata, resolution = 3)
    sc.tl.umap(adata, min_dist=0.1)
    adata.obs['umap_1'] = list(adata.obsm['X_umap'][:,0])
    adata.obs['umap_2'] = list(adata.obsm['X_umap'][:,1])
    name = 'adata.obs.louvain_' + cluster + '.csv'
    adata.obs.to_csv(os.path.join('doublet_cluster', name))




