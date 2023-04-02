
##########################################
### running Scanpy to perform analysis ###
##########################################

import scanpy as sc
import pandas as pd
import numpy as np
import os
import time
import gc
import sys

start_time = time.time()

WORK_PATH = '/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mtx'

trajectory_list = ["CNS_neurons", "Brain_and_spinal_cord"]

day_list = ["E85","E9","E95","E10","E11","E12"]

celltype_list_exclude = ["Amacrine cells", "Amacrine/Horizontal precursor cells", "Cholinergic amacrine cells", "Horizontal cells", "PV-containing retinal ganglion cells", "Retinal ganglion cells"]

example_id = "Neurons_plus"
print(example_id)

fdata = pd.read_csv("/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mtx/df_gene.csv", index_col = 0)

adata_1 = sc.read_h5ad(os.path.join(WORK_PATH, 'E85_adata_scale.h5ad'))
adata_2 = sc.read_h5ad(os.path.join(WORK_PATH, 'E9_adata_scale.h5ad'))
adata_3 = sc.read_h5ad(os.path.join(WORK_PATH, 'E95_adata_scale.h5ad'))
adata_4 = sc.read_h5ad(os.path.join(WORK_PATH, 'E10_adata_scale.h5ad'))
adata_5 = sc.read_h5ad(os.path.join(WORK_PATH, 'E11_adata_scale.h5ad'))
adata_6 = sc.read_h5ad(os.path.join(WORK_PATH, 'E12_adata_scale.h5ad'))
adata_7 = sc.read_h5ad(os.path.join(WORK_PATH, 'E13_adata_scale.h5ad'))
adata_8 = sc.read_h5ad(os.path.join(WORK_PATH, 'E14_adata_scale.h5ad'))
adata_9 = sc.read_h5ad(os.path.join(WORK_PATH, 'E15_adata_scale.h5ad'))
adata_10 = sc.read_h5ad(os.path.join(WORK_PATH, 'E16_adata_scale.h5ad'))
adata_11 = sc.read_h5ad(os.path.join(WORK_PATH, 'E17_adata_scale.h5ad'))
adata_12 = sc.read_h5ad(os.path.join(WORK_PATH, 'E18_adata_scale.h5ad'))
adata_13 = sc.read_h5ad(os.path.join(WORK_PATH, 'P0_adata_scale.h5ad'))

adata = adata_1.concatenate(adata_2, adata_3, adata_4, adata_5, adata_6, adata_7, adata_8, adata_9, adata_10, adata_11, adata_12, adata_13)

del adata_1, adata_2, adata_3, adata_4, adata_5, adata_6, adata_7, adata_8, adata_9, adata_10, adata_11, adata_12, adata_13
gc.collect()

pdata = pd.read_csv(os.path.join(WORK_PATH, 'backup','adata_scale.full.obs.csv'), index_col = 0)
adata.obs = pdata
adata = adata[adata.obs["major_trajectory"].isin(trajectory_list)]
adata = adata[adata.obs["group"].isin(day_list)]
adata = adata[~adata.obs["celltype_update"].isin(celltype_list_exclude)]

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

### performing regress_out and scale ###
### sc.pp.regress_out(adata, ['log2_umi', 'S.Score', 'G2M.Score'])
### print("Done regressing out co-factors ...")
### print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

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

sc.tl.leiden(adata, resolution=5, n_iterations=2)
adata.obs['leiden_res_5'] = adata.obs['leiden']
print("Done clustering using res = 5 ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.leiden(adata, resolution=20, n_iterations=2)
adata.obs['leiden_res_20'] = adata.obs['leiden']
print("Done clustering using res = 20 ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

sc.tl.umap(adata, min_dist=0.3, n_components=2)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])

adata.obs.to_csv(os.path.join(WORK_PATH, 'example', '%s_adata_scale.obs.csv'%example_id))

adata.write(os.path.join(WORK_PATH, 'example', '%s_adata_scale.h5ad'%example_id))
print("Done writing data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

emb = adata.obsm['X_pca']
print(emb.shape)
np.savetxt(os.path.join(WORK_PATH, 'example', '%s_adata_scale.PCs.csv'%example_id), emb, delimiter=",", fmt='%1.3f')


