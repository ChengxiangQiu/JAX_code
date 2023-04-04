
#############################################################################################################
### Here, we peformed basic analysis (normalization, dimension reuction, and clustering) on 11.4M dataset ###
#############################################################################################################

### Of note, this data should include 24,552 genes and 11,441,407 cells
### Hint: I suggest to request >500GB for the following analysis.

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

sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
print("Done performing PCA ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)
print("Done computing neighborhood graph ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.umap(adata, min_dist=0.1, n_components=3)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

sc.tl.umap(adata, min_dist=0.1, n_components=2)
print("Done UMAP ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])

sc.tl.leiden(adata, resolution=1, n_iterations=2)
adata.obs['leiden_res_1'] = adata.obs['leiden']
print("Done clustering using res = 1 ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

sc.tl.leiden(adata, resolution=2, n_iterations=2)
adata.obs['leiden_res_2'] = adata.obs['leiden']
print("Done clustering using res = 2 ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

adata.obs.to_csv(os.path.join(WORK_PATH, 'adata_scale.obs.csv'))
pd.DataFrame(adata.obsm['X_pca']).to_csv(os.path.join(WORK_PATH, 'adata_scale.pca.csv'))

adata.write(os.path.join(WORK_PATH, 'adata_scale.h5ad'), compression="gzip")
print("Done writing data ...")
print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')


#################################################################################################
### For each major trajectories, we further performed sub-clustering to get higher resolution ###
#################################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import os
import time
import sys

start_time = time.time()

WORK_PATH = './'

adata_1 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_1.h5ad'))
adata_2 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_2.h5ad'))
adata_3 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_3.h5ad'))
adata_4 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_4.h5ad'))

adata_orig = adata_1.concatenate(adata_2, adata_3, adata_4)
del adata_1, adata_2, adata_3, adata_4
gc.collect()

### Of note, please read df_cell.rds and then write it to df_cell.csv in R
pdata = pd.read_csv(os.path.join(WORK_PATH, 'df_cell.csv'), index_col = 0)
adata_orig.obs = pdata

trajectory_list = ["Neuroectoderm_and_glia",
"Intermediate_neuronal_progenitors",
"Eye_and_other",
"Ependymal_cells",
"CNS_neurons",
"Mesoderm",
"Definitive_erythroid",
"Epithelium",
"Endothelium",
"Muscle_cells",
"Hepatocytes",
"White_blood_cells",
"Neural_crest_PNS_glia",
"Adipocytes",
"Primitive_erythroid",
"Neural_crest_PNS_neurons",
"T_cells",
"Lung_and_airway",
"Intestine",
"B_cells",
"Olfactory_sensory_neurons",
"Cardiomyocytes",
"Oligodendrocytes",
"Mast_cells",
"Megakaryocytes",
"Testis_and_adrenal"]

for trajectory_id in trajectory_list:
    print("Processing: %s"%trajectory_id)

    adata = adata_orig[adata_orig.obs["major_trajectory"] == trajectory_id]

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

    adata.obs['UMAP_1'] = list(adata.obsm['X_umap'][:,0])
    adata.obs['UMAP_2'] = list(adata.obsm['X_umap'][:,1])
    adata.obs['UMAP_3'] = list(adata.obsm['X_umap'][:,2])

    sc.tl.umap(adata, min_dist=0.1, n_components=2)
    print("Done UMAP ...")
    print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

    adata.obs['UMAP_2d_1'] = list(adata.obsm['X_umap'][:,0])
    adata.obs['UMAP_2d_2'] = list(adata.obsm['X_umap'][:,1])

    sc.tl.leiden(adata, resolution=1, n_iterations=2)
    adata.obs['leiden_res_1'] = adata.obs['leiden']
    print("Done clustering using res = 1 ...")
    print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

    sc.tl.leiden(adata, resolution=5, n_iterations=2)
    adata.obs['leiden_res_5'] = adata.obs['leiden']
    print("Done clustering using res = 5 ...")
    print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')

    adata.obs.to_csv(os.path.join(WORK_PATH, '%s_adata_scale.obs.csv'%trajectory_id))
    pd.DataFrame(adata.obsm['X_pca']).to_csv(os.path.join(WORK_PATH, '%s_adata_scale.pca.csv'%trajectory_id))

    adata.write(os.path.join(WORK_PATH, '%s_adata_scale.h5ad'%trajectory_id), compression="gzip")
    print("Done writing data ...")
    print(str(format((time.time() - start_time)/3600, '.4f')) + 'hours')






