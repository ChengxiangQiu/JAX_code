import os, sys
import numpy as np
import pandas as pd
import scanpy as sc

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mosta"
data_path = "/net/shendure/vol10/projects/cxqiu/nobackup/data/stereo_seq"

file = open(os.path.join(work_path, "file_list.txt"))
file_list = [line.rstrip().replace(".MOSTA.h5ad", "") for line in file]
file.close()

file_id = file_list[int(sys.argv[1]) - 1]
print(file_id)

mosta = sc.read_h5ad(os.path.join(data_path, '%s.MOSTA.h5ad'%file_id))
mosta_meta = mosta.obs

np.savetxt(os.path.join(work_path, 'annotation', '%s.spatial_color.csv'%file_id), mosta.uns['annotation_colors'], delimiter=",", fmt='%s')
np.savetxt(os.path.join(work_path, 'annotation', '%s.spatial_color_id.csv'%file_id), mosta.obs['annotation'].cat.categories, delimiter=",", fmt='%s')
np.savetxt(os.path.join(work_path, 'annotation', '%s.spatial_coor.csv'%file_id), mosta.obsm['spatial'], delimiter=",")
np.savetxt(os.path.join(work_path, 'annotation', '%s.spatial_anno.csv'%file_id), mosta.obs['annotation'], delimiter=",", fmt='%s')


