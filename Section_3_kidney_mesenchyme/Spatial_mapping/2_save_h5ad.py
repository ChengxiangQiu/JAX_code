import os, sys
import numpy as np
import pandas as pd
import scanpy as sc

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mosta"
data_path = "/net/shendure/vol10/projects/cxqiu/nobackup/data/stereo_seq"

day_list = ["E95","E105","E115","E125","E135","E145","E155","E165"]

day_id = day_list[int(sys.argv[1]) - 1]
adata = sc.read_mtx(os.path.join(work_path, 'sc_data', '%s.gene_count.mtx'%day_id))
pdata = pd.read_csv(os.path.join(work_path, 'sc_data', '%s.df_cell.csv'%day_id), index_col = 0)
fdata = pd.read_csv(os.path.join(work_path, 'sc_data', '%s.df_gene.csv'%day_id), index_col = 0)
adata.obs = pdata
adata.var = fdata

adata.write(os.path.join(work_path, 'sc_data', '%s.sc_data.h5ad'%day_id)) 
