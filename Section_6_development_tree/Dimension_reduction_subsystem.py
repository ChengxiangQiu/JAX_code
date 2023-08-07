
##########################################################################################
### First, we manually split all the cell types, from the organogenesis & fetal development, into 12 systems, 
### to perform dimension reducting using Scanpy, followed by identifying the kNNs across cells using annoy in Python
##########################################################################################


import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import time
import gc
from annoy import AnnoyIndex

start_time = time.time()

WORK_PATH = './'

adata_1 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_1.h5ad'))
adata_2 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_2.h5ad'))
adata_3 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_3.h5ad'))
adata_4 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_4.h5ad'))

adata_orig = adata_1.concatenate(adata_2, adata_3, adata_4)
del adata_1, adata_2, adata_3, adata_4
gc.collect()

##############################################
### Of note, please read df_cell_graph.rds and then write it to df_cell_graph.csv in R

### >>> dat = readRDS("df_cell_graph.rds")
### >>> rownames(dat) = as.vector(dat$cell_id)
### >>> write.csv(dat, "df_cell_graph.csv")

pdata = pd.read_csv(os.path.join(WORK_PATH, 'df_cell_graph.csv'), index_col = 0)
adata_orig.obs = pdata

system_list = ["Endothelium",
"Epithelial_cells",
"Eye",
"Gut",
"Notochord",
"PNS_glia",
"PNS_neurons",
"Renal",
"Lateral_plate_mesoderm",
"Blood",
"Brain_spinal_cord",
"Mesoderm"]

for system_i in trajectory_list:

    print("Processing: %s"%system_i)

    adata = adata_orig[adata_orig.obs["system"] == system_i]
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

    sc.tl.leiden(adata, resolution=2, n_iterations=2)
    adata.obs['leiden_res_2'] = adata.obs['leiden']

    sc.tl.leiden(adata, resolution=5, n_iterations=2)
    adata.obs['leiden_res_5'] = adata.obs['leiden']

    sc.tl.leiden(adata, resolution=10, n_iterations=2)
    adata.obs['leiden_res_10'] = adata.obs['leiden']

    adata.obs.to_csv(os.path.join(WORK_PATH, '%s_adata_scale.obs.csv'%system_i))

    adata.write(os.path.join(WORK_PATH, '%s_adata_scale.h5ad'%system_i), compression="gzip")

    X = adata.obsm['X_pca']
    print(X.shape)
    np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%system_i), X, delimiter=",", fmt='%1.3f')

    ### calculating kNN using annoy, this is much faster than using R

    dist_metric = 'euclidean'
    k = 15
    ### Here, why we use 15? because the log2(mean cell number across cell types) is around 15.

    npc = X.shape[1]
    ncell = X.shape[0]
    annoy_index = AnnoyIndex(npc, metric=dist_metric)

    for i in range(ncell):
        annoy_index.add_item(i, list(X[i,:]))
    annoy_index.build(15) ### bigger number will make the result more accurate

    knn = []
    for iCell in range(ncell):
        knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
    knn = np.array(knn, dtype=int)

    np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.kNN_15.csv'%system_i), knn, delimiter=",", fmt='%s')



######################################################
### Second, we found the neuroectoderm is too complex, so we subset the patterned neuroectoderm cells to perform embedding


import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import time
import gc
from annoy import AnnoyIndex

start_time = time.time()

WORK_PATH = './'

adata_1 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_1.h5ad'))
adata_2 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_2.h5ad'))
adata_3 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_3.h5ad'))
adata_4 = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_JAX_dataset_4.h5ad'))

adata_orig = adata_1.concatenate(adata_2, adata_3, adata_4)
del adata_1, adata_2, adata_3, adata_4
gc.collect()

##############################################
### Of note, please read df_cell_graph.rds and then write it to df_cell_graph.csv in R

### >>> dat = readRDS("df_cell_graph.rds")
### >>> rownames(dat) = as.vector(dat$cell_id)
### >>> write.csv(dat, "df_cell_graph.csv")

pdata = pd.read_csv(os.path.join(WORK_PATH, 'df_cell_graph.csv'), index_col = 0)
adata_orig.obs = pdata

patterned_neuroectoderm = ["Anterior floor plate",
"Diencephalon",
"Floorplate and p3 domain",
"Hypothalamus",
"Midbrain",
"Posterior roof plate",
"Telencephalon",
"Anterior roof plate",
"Dorsal telencephalon",
"Hindbrain",
"Hypothalamus (Sim1+)",
"Midbrain-hindbrain boundary",
"Spinal cord/r7/r8"]

system_i = "Neuroectoderm"
print("Processing: %s"%system_i)

adata = adata_orig[adata_orig.obs["celltype_update"].isin(patterned_neuroectoderm)]
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

sc.tl.leiden(adata, resolution=2, n_iterations=2)
adata.obs['leiden_res_2'] = adata.obs['leiden']

sc.tl.leiden(adata, resolution=5, n_iterations=2)
adata.obs['leiden_res_5'] = adata.obs['leiden']

sc.tl.leiden(adata, resolution=10, n_iterations=2)
adata.obs['leiden_res_10'] = adata.obs['leiden']

adata.obs.to_csv(os.path.join(WORK_PATH, '%s_adata_scale.obs.csv'%system_i))

adata.write(os.path.join(WORK_PATH, '%s_adata_scale.h5ad'%system_i), compression="gzip")

X = adata.obsm['X_pca']
print(X.shape)
np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%system_i), X, delimiter=",", fmt='%1.3f')

### calculating kNN using annoy, this is much faster than using R

dist_metric = 'euclidean'
k = 15
### Here, why we use 15? because the log2(mean cell number across cell types) is around 15.

npc = X.shape[1]
ncell = X.shape[0]
annoy_index = AnnoyIndex(npc, metric=dist_metric)

for i in range(ncell):
    annoy_index.add_item(i, list(X[i,:]))
annoy_index.build(15) ### bigger number will make the result more accurate

knn = []
for iCell in range(ncell):
    knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
knn = np.array(knn, dtype=int)

np.savetxt(os.path.join(WORK_PATH, '%s_adata_scale.kNN_15.csv'%system_i), knn, delimiter=",", fmt='%s')


