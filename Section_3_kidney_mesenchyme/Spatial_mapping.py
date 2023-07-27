
#######################################################################################
### Part-1: extracting spatial coordinates and annotations from individual sections ###
#######################################################################################

import os, sys
import numpy as np
import pandas as pd
import scanpy as sc

WORK_PATH = "./"

newpath = os.path.join(WORK_PATH, 'annotation')
if not os.path.exists(newpath):
    os.makedirs(newpath)

file = open(os.path.join(WORK_PATH, "Mosta_file_list.txt"))
file_list = [line.rstrip().replace(".MOSTA.h5ad", "") for line in file]
file.close()

for file_id in file_list:

    mosta = sc.read_h5ad(os.path.join(WORK_PATH, '%s.MOSTA.h5ad'%file_id))
    mosta_meta = mosta.obs

    np.savetxt(os.path.join(WORK_PATH, 'annotation', '%s.spatial_color.csv'%file_id), mosta.uns['annotation_colors'], delimiter=",", fmt='%s')
    np.savetxt(os.path.join(WORK_PATH, 'annotation', '%s.spatial_color_id.csv'%file_id), mosta.obs['annotation'].cat.categories, delimiter=",", fmt='%s')
    np.savetxt(os.path.join(WORK_PATH, 'annotation', '%s.spatial_coor.csv'%file_id), mosta.obsm['spatial'], delimiter=",")
    np.savetxt(os.path.join(WORK_PATH, 'annotation', '%s.spatial_anno.csv'%file_id), mosta.obs['annotation'], delimiter=",", fmt='%s')


##################################################################
### Part-2: generating h5ad format profile for sc-RNA-seq data ###
##################################################################

import os, sys
import numpy as np
import pandas as pd
import scanpy as sc

WORK_PATH = "./"

day_list = ["E95","E105","E115","E125","E135","E145","E155","E165"]

for day_id in day_list:

    adata = sc.read_mtx(os.path.join(WORK_PATH, 'sc_data', '%s.gene_count.mtx'%day_id))
    pdata = pd.read_csv(os.path.join(WORK_PATH, 'sc_data', '%s.df_cell.csv'%day_id), index_col = 0)
    fdata = pd.read_csv(os.path.join(WORK_PATH, 'sc_data', '%s.df_gene.csv'%day_id), index_col = 0)
    adata.obs = pdata
    adata.var = fdata

    adata.write(os.path.join(WORK_PATH, 'sc_data', '%s.sc_data.h5ad'%day_id)) 


########################################################
### Part-3: performing spatial mapping using Tangram ###
########################################################

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
import tangram as tg

WORK_PATH = "./"

newpath = os.path.join(WORK_PATH, 'result')
if not os.path.exists(newpath):
    os.makedirs(newpath)

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print('Using device:', device)
print(torch.cuda)
print(torch.version.cuda)
print(torch.cuda.is_available())

file = open(os.path.join(work_path, "MOSTA_file_list.txt"))
file_list = [line.rstrip().replace(".MOSTA.h5ad", "") for line in file]
file.close()


for spatial_id in file_list:

    print(spatial_id)

    mosta = sc.read(os.path.join(WORK_PATH, spatial_id + '.MOSTA.h5ad'))

    day_id = spatial_id.split('_')[0].replace('.','')
    adata = sc.read(os.path.join(WORK_PATH, 'sc_data', '%s.sc_data.h5ad'%day_id))

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

    np.savetxt(os.path.join(WORK_PATH, 'result', '%s.result.csv'%spatial_id), mosta.obsm['tangram_ct_pred'], delimiter=",", fmt='%.8e', header=colnames)
    np.savetxt(os.path.join(WORK_PATH, 'result', '%s.result.coor.csv'%spatial_id), mosta.obsm['spatial'], delimiter=",")




















