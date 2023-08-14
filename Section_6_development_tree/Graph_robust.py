
################################################################################
### Second, to assess the robustness of MNNs to cell sampling, we randomly subsampled 80% 
### of cells from each developmental system during organogenesis & fetal development

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import time
import gc
from annoy import AnnoyIndex

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
"Mesoderm",
"Neuroectoderm"]

WORK_PATH = "./"

for system_i in system_list:

    ### PC features were calculated by Dimension_reduction_subsystem.py

    X = pd.read_csv(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%system_i), index_col = False, header=None)
    X = pd.DataFrame.to_numpy(X)

    original_size = X.shape[0]
    subset_size = int(original_size * 0.8)
    npc = X.shape[1]

    dist_metric = 'euclidean'
    k = 15
    ### Here, why we use 15? because the log2(mean cell number across cell types) is around 15.

    for cnt in range(100):
        
        idx = np.random.choice(original_size, size=subset_size, replace=False)
        X_sub = X[idx,:]

        ncell = X_sub.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)

        for i in range(ncell):
            annoy_index.add_item(i, list(X_sub[i,:]))
        annoy_index.build(15) ### bigger number will make the result more accurate

        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

        np.savetxt(os.path.join(WORK_PATH, '%s_knn_%s.csv'%(system_i, str(cnt+1))), knn, delimiter=",", fmt='%s')
        np.savetxt(os.path.join(WORK_PATH, '%s_idx_%s.csv'%(system_i, str(cnt+1))), idx, delimiter=",", fmt='%s')




##################################################################################
### Third, to determine the effect of k parameter choice on the MNNs identified 
### between cell types, we examined different k values (k = 5, 10, 20, 30, 40, 50)
##################################################################################


import scanpy as sc
import pandas as pd
import numpy as np
import os, sys
import time
import gc
from annoy import AnnoyIndex

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
"Mesoderm",
"Neuroectoderm"]

WORK_PATH = "./"

for system_i in system_list:

    ### PC features were calculated by Dimension_reduction_subsystem.py

    X = pd.read_csv(os.path.join(WORK_PATH, '%s_adata_scale.PCs.csv'%system_i), index_col = False, header=None)
    X = pd.DataFrame.to_numpy(X)

    ncell = X.shape[0]
    npc = X.shape[1]
    dist_metric = 'euclidean'

    annoy_index = AnnoyIndex(npc, metric=dist_metric)

    for i in range(ncell):
        annoy_index.add_item(i, list(X[i,:]))
    annoy_index.build(15) ### bigger number will make the result more accurate

    for k in [5,10,20,30,40,50]:

        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

        np.savetxt(os.path.join(WORK_PATH, '%s_knn_%s.csv'%(system_i, str(k))), knn, delimiter=",", fmt='%s')



















