

### pull out Fgf17 positive cells

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys

work_path = '/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/Neurons'
example_id = "Cerebellar_Purkinje"

adata = sc.read_mtx(os.path.join(work_path, '%s.mtx'%example_id))
pdata = pd.read_csv(os.path.join(work_path, '%s.df_cell.csv'%example_id), index_col = 0)
fdata = pd.read_csv("/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mtx/df_gene.csv", index_col = 0)
adata.obs = pdata
adata.var = fdata

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

sc.tl.leiden(adata, resolution=1, n_iterations=2)
adata.obs['leiden_res_1'] = adata.obs['leiden']

adata.obs.to_csv(os.path.join(work_path, '%s_adata_scale.obs.csv'%example_id))

emb = adata.obsm['X_pca']
np.savetxt(os.path.join(work_path, '%s_adata_scale.PCs.csv'%example_id), emb, delimiter=",", fmt='%1.3f')
