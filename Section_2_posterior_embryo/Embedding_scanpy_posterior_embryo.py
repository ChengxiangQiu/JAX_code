
###################################################################################################################################################
### Here, we peformed basic analysis (normalization, dimension reuction, and clustering) on NMPs, gut, and notochord during early somitogenesis ###
###################################################################################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import time
import gc

start_time = time.time()

WORK_PATH = './'

adata_1 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_1.h5ad'))
adata_2 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_2.h5ad'))
adata_3 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_3.h5ad'))
adata_4 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_4.h5ad'))

adata = adata_1.concatenate(adata_2, adata_3, adata_4)
del adata_1, adata_2, adata_3, adata_4
gc.collect()

pdata = pd.read_csv(os.path.join(WORK_PATH, 'df_cell.csv'), index_col = 0)
adata.obs = pdata

day_include = ["E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0"]
celltype_include = ["Notochord", "Nodal cilia", "Neuromesodermal progenitors", "Gut", "Mesodermal progenitors (Tbx6+)"]

example_id = "posterior_embryo"
print(example_id)

adata = adata[adata.obs["day"].isin(day_include)]
adata = adata[adata.obs["celltype_update"].isin(celltype_include)]

adata.write(os.path.join(WORK_PATH, '%s_adata_scale.h5ad'%example_id), compression="gzip")

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

adata.obs.to_csv(os.path.join(WORK_PATH, '%s_adata_scale.obs.replicate.csv'%example_id))

adata.write(os.path.join(WORK_PATH, '%s_adata_scale_processed.h5ad'%example_id), compression="gzip")
print("Done writing data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

emb = adata.obsm['X_pca']
print(emb.shape)
np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%example_id), emb, delimiter=",", fmt='%1.3f')


#########################################################
### Perform subclustering on three major trajectories ###
#########################################################

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys

work_path = './'
example_id = "posterior_embryo"

adata_all = sc.read_h5ad(os.path.join(work_path, '%s_adata_scale.h5ad'%example_id))
pdata = pd.read_csv(os.path.join(work_path, '%s_adata_scale.obs.csv'%example_id), index_col = 0)
adata_all.obs = pdata

subcluster_list = ["NMP_Mesoderm", "Notochord", "Gut"]

for i in subcluster_list:

    adata = adata_all[adata_all.obs["cluster_tmp"] == i]
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

    sc.tl.leiden(adata, resolution=1, n_iterations=2)
    adata.obs['leiden_res_1'] = adata.obs['leiden']

    adata.obs.to_csv(os.path.join(work_path, '%s_adata_scale.%s.obs.csv'%(example_id, i)))
    
    emb = adata.obsm['X_pca']
    np.savetxt(os.path.join(work_path, '%s_adata_scale.%s.PCs.csv'%(example_id, i)), emb, delimiter=",", fmt='%1.3f')


