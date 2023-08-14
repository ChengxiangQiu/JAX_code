
################################################
### run geosketch to downsample to 10% cells ###
################################################

import numpy as np
import pandas as pd
import scanpy as sc
import os, sys
from time import time
from geosketch import gs

start_time = time.time()

WORK_PATH = './'

### Of note, this is the adata object after running Dimension_reduction_scanpy.py in Section-1
### We need the PCA features to perform geosketch

adata = sc.read_h5ad(os.path.join(WORK_PATH, 'adata_scale.h5ad'))

X_dimred = adata.obsm['X_pca']

start = time()

N = 1144141 # Number of samples to obtain from the data set.
sketch_index = gs(X_dimred, N, replace=False)

np.savetxt(os.path.join(WORK_PATH, "adata_scale_geosketch_downsample.csv"), np.array(sketch_index)+1, delimiter=",", fmt='%s')
adata.obs.to_csv(os.path.join(WORK_PATH, 'adata_scale.obs.csv'))

end = time()
print(end - start)


