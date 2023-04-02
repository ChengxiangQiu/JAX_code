import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
import tangram as tg

mosta_data_path = "/net/shendure/vol10/projects/cxqiu/nobackup/data/stereo_seq"
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mosta"

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print('Using device:', device)
print(torch.cuda)
print(torch.version.cuda)
print(torch.cuda.is_available())

file = open(os.path.join(work_path, "file_list.txt"))
file_list = [line.rstrip().replace(".MOSTA.h5ad", "") for line in file]
file.close()


for spatial_id in file_list[18:47]:

    print(spatial_id)

    mosta = sc.read(os.path.join(mosta_data_path, spatial_id + '.MOSTA.h5ad'))

    day_id = spatial_id.split('_')[0].replace('.','')
    adata = sc.read(os.path.join(work_path, 'sc_data', '%s.sc_data.h5ad'%day_id))

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2500)

    if mosta.shape[0] > 90000:
        sc.pp.subsample(mosta, n_obs = 90000)

    var_genes = adata.var.index[adata.var.highly_variable]
    tg.pp_adatas(adata, mosta, genes=var_genes)

    ad_map = tg.map_cells_to_space(
        adata_sc=adata,
        adata_sp=mosta,
        device='cuda:0'
    )

    tg.project_cell_annotations(ad_map, mosta, annotation='celltype')
    annotation_list = list(pd.unique(adata.obs['celltype']))

    colnames = ','.join(list(mosta.obsm['tangram_ct_pred'].columns))

    np.savetxt(os.path.join(work_path, 'result', '%s.result.csv'%spatial_id), mosta.obsm['tangram_ct_pred'], delimiter=",", fmt='%.8e', header=colnames)
    np.savetxt(os.path.join(work_path, 'result', '%s.result.coor.csv'%spatial_id), mosta.obsm['spatial'], delimiter=",")

