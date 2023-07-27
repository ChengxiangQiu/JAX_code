
############################################################################################################
### Here, we peformed basic analysis (normalization, dimension reuction, and clustering) on renal subset ###
############################################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import time
import gc

start_time = time.time()

WORK_PATH = './'

example_id = "LPM_somite_5_20"
print(example_id)


### First, only including backbone cells

adata = sc.read_mtx(os.path.join(WORK_PATH, '%s_backbone.gene_count.mtx'%example_id))
pdata = pd.read_csv(os.path.join(WORK_PATH, '%s_backbone.df_cell.csv'%example_id), index_col = 0)
fdata = pd.read_csv(os.path.join(WORK_PATH, "df_gene.csv"), index_col = 0)
adata.obs = pdata
adata.var = fdata

print("Done reading data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.normalize_total(adata, target_sum=1e4)
print("Done normalization by total counts ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.log1p(adata)
print("Done log transformation ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.highly_variable_genes(adata, n_top_genes=2500)
print("Done finding highly variable genes ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata = adata[:, adata.var.highly_variable]
print("Done filtering in highly variable genes ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.scale(adata, max_value=10)
print("Done scaling data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')
### done with regress_out and scale ###

sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
print("Done performing PCA ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
print("Done computing neighborhood graph ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.umap(adata, min_dist=0.3, n_components=3)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.leiden(adata, resolution=1, n_iterations=2)
adata.obs['leiden_res_1'] = adata.obs['leiden']
print("Done clustering using res = 1 ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

sc.tl.umap(adata, min_dist=0.3, n_components=2)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])

adata.obs.to_csv(os.path.join(WORK_PATH, '%s_backbone_adata_scale.obs.csv'%example_id))

adata.write(os.path.join(WORK_PATH, '%s_backbone_adata_scale_processed.h5ad'%example_id), compression="gzip")
print("Done writing data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

emb = adata.obsm['X_pca']
print(emb.shape)
np.savetxt(os.path.join(WORK_PATH, '%s_backbone_adata_scale.PCs.csv'%example_id), emb, delimiter=",", fmt='%1.3f')



### Next, including both backbone cells and derivatives

adata = sc.read_mtx(os.path.join(WORK_PATH, '%s_all.gene_count.mtx'%example_id))
pdata = pd.read_csv(os.path.join(WORK_PATH, '%s_all.df_cell.csv'%example_id), index_col = 0)
fdata = pd.read_csv(os.path.join(WORK_PATH, "df_gene.csv"), index_col = 0)
adata.obs = pdata
adata.var = fdata

print("Done reading data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.normalize_total(adata, target_sum=1e4)
print("Done normalization by total counts ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.log1p(adata)
print("Done log transformation ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.highly_variable_genes(adata, n_top_genes=2500)
print("Done finding highly variable genes ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata = adata[:, adata.var.highly_variable]
print("Done filtering in highly variable genes ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.scale(adata, max_value=10)
print("Done scaling data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')
### done with regress_out and scale ###

sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
print("Done performing PCA ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
print("Done computing neighborhood graph ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.umap(adata, min_dist=0.3, n_components=3)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.leiden(adata, resolution=1, n_iterations=2)
adata.obs['leiden_res_1'] = adata.obs['leiden']
print("Done clustering using res = 1 ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

sc.tl.umap(adata, min_dist=0.3, n_components=2)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])

adata.obs.to_csv(os.path.join(WORK_PATH, '%s_all_adata_scale.obs.csv'%example_id))

adata.write(os.path.join(WORK_PATH, '%s_all_adata_scale_processed.h5ad'%example_id), compression="gzip")
print("Done writing data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

emb = adata.obsm['X_pca']
print(emb.shape)
np.savetxt(os.path.join(WORK_PATH, '%s_all_adata_scale.PCs.csv'%example_id), emb, delimiter=",", fmt='%1.3f')

