
##################################
### Section - 5, Neuroectoderm ###
##################################

########################################################################
### Mapping patterned neuroectoderm (backbone) and their derivatives ###
########################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Neuroectoderm_derivative"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
pd_x = pd %>% left_join(pd_all[,c("cell_id", "neurons_sub_clustering")])
pd_x = pd_x[!is.na(pd_x$neurons_sub_clustering),]

### manually merging some cell types

pd_1 = pd[pd$cell_id %in% pd_x$cell_id,]
pd_1 = pd_1[as.vector(pd_x$cell_id),]
pd_1$celltype_update = as.vector(pd_x$neurons_sub_clustering)

pd_2 = pd[pd$major_trajectory == "Intermediate_neuronal_progenitors",]
pd_2$celltype_update = "Intermediate neuronal progenitors"

pd_3 = pd[pd$major_trajectory == "Ependymal_cells",]
pd_3$celltype_update = "Choroid plexus"

pd_4 = pd[pd$major_trajectory == "Oligodendrocytes",]
pd_4$celltype_update = "Oligodendrocytes"

pd_5 = pd[(!pd$cell_id %in% pd_x$cell_id) & 
              (!pd$major_trajectory %in% c("Intermediate_neuronal_progenitors", "Ependymal_cells", "Oligodendrocytes")),]

pd_y = rbind(pd_1, pd_2, pd_3, pd_4, pd_5)
pd_y = pd_y[rownames(pd),]
pd = pd_y

### Based on Knn identified in PCA space (n=30 dimensions), find the nearest neighboring cells for individual derivatives

emb = readRDS(paste0(work_path, example_i, "_adata_scale.PCs.rds"))

pd$celltype = as.vector(pd$celltype_update)

pd_integration = pd
PC_integration = emb

day_value = as.vector(pd_integration$somite_count)
day_value[pd_integration$day == "E10.25"] = "35 somites"
day_value[pd_integration$day == "E10.5"]  = "36 somites"
day_value[pd_integration$day == "E10.75"] = "37 somites"
day_value[pd_integration$day == "E11.0"]  = "38 somites"
day_value[pd_integration$day == "E11.25"] = "39 somites"
day_value[pd_integration$day == "E11.5"]  = "40 somites"
day_value[pd_integration$day == "E11.75"] = "41 somites"
day_value[pd_integration$day == "E12.0"]  = "42 somites"
day_value[pd_integration$day == "E12.25"] = "43 somites"
day_value[pd_integration$day == "E12.5"]  = "44 somites"
day_value[pd_integration$day == "E12.75"] = "45 somites"
pd_integration$somite_count = as.vector(day_value)

day_x = data.frame(somite_count = paste0(c(0, 2:12, 14:18, 20:45), " somites"),
                   day_value = c(1:length(c(0, 2:12, 14:18, 20:45))))
day_y = pd_integration %>% left_join(day_x, by = "somite_count")
pd_integration$day_value = as.vector(day_y$day_value)

### cell type of each column (brain + spinal cord, the origins)
col_cell_type = c("Anterior floor plate",
                  "Diencephalon",
                  "Hypothalamus",
                  "Midbrain",
                  "Floorplate and p3 domain",
                  "Telencephalon",
                  "Anterior roof plate",
                  "Dorsal telencephalon",
                  "Hindbrain",
                  "Hypothalamus (Sim1+)",
                  "Midbrain-hindbrain boundary",
                  "Posterior roof plate",
                  "Spinal cord/r7/r8")

### cell type of each row (neurons, the derivatives)
row_cell_type = c("Astrocytes",
                  "Intermediate neuronal progenitors",
                  "Choroid plexus",
                  "Oligodendrocytes",
                  "Cajal-Retzius cells",
                  "Cerebellar Purkinje cells",
                  "Cranial motor neurons",
                  "Di/mesencephalon GABAergic neurons",
                  "Di/mesencephalon glutamatergic neurons",
                  "GABAergic cortical interneurons",
                  "Hypothalamic Sim1 neurons",
                  "Midbrain dopaminergic neurons",
                  "Neural progenitor cells (Neurod1+)",
                  "Neural progenitor cells (Ror1+)",
                  "Neurons (Slc17a8+)",
                  "Spinal cord motor neuron progenitors",
                  "Precerebellar neurons",
                  "Spinal cord motor neurons",
                  "Spinal dI1 interneurons",
                  "Spinal dI2 interneurons",
                  "Spinal dI3 interneurons",
                  "Spinal dI4 interneurons",
                  "Spinal dI5 interneurons",
                  "Spinal dI6 interneurons",
                  "Spinal V0 interneurons",
                  "Spinal V1 interneurons",
                  "Spinal V2a interneurons",
                  "Spinal V2b interneurons",
                  "Spinal V3 interneurons",
                  "Striatal projection neurons",
                  "Thalamic neuronal precursors")

pd_1 = pd_integration[pd_integration$celltype %in% row_cell_type,]
pd_2 = pd_integration[pd_integration$celltype %in% col_cell_type,]

pd_1_x = pd_1 %>% group_by(celltype) %>% slice_min(order_by = day_value, n = 500, with_ties = F)
pd_1_y = pd_1_x %>% group_by(celltype) %>% tally() %>% filter(n >= 500) %>% pull(celltype)
pd_1_x = pd_1_x %>% filter(celltype %in% pd_1_y)

pd_1 = pd_integration[pd_integration$cell_id %in% pd_1_x$cell_id,]

pd_2_x = pd_2
pd_2 = pd_integration[pd_integration$cell_id %in% pd_2_x$cell_id,]

emb_1 = PC_integration[pd_integration$cell_id %in% pd_1_x$cell_id,]
emb_2 = PC_integration[pd_integration$cell_id %in% pd_2_x$cell_id,]

### mutually nearest neighbors

k.param = 10; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = emb_2,
    query = emb_1,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked

k.param = 10; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = emb_1,
    query = emb_2,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix_2 = nn.ranked

### save this result because it will be used for step6_Key_TFs.R
saveRDS(list(pd_1 = pd_1, pd_2 = pd_2, 
             nn_matrix = nn_matrix, nn_matrix_2 = nn_matrix_2), paste0(work_path, example_i, ".MNN.rds"))

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(pd_2$cell_id)[as.vector(nn_matrix[,i])])
}
datA = data.frame(A = rep(pd_1$cell_id, k.param),
                  B = c(resultA), 
                  value = 1, stringsAsFactors = F)

resultB = NULL
for(i in 1:k.param){
    print(i)
    resultB = cbind(resultB, as.vector(pd_1$cell_id)[as.vector(nn_matrix_2[,i])])
}
datB = data.frame(A = c(resultB),
                  B = rep(pd_2$cell_id, k.param), 
                  value = 2, stringsAsFactors = F)

dat = datB %>% left_join(datA, by = c("A" = "A", "B" = "B")) %>%
    filter(!is.na(value.y)) %>% select(A, B) %>%
    left_join(pd_1 %>% mutate(A = cell_id) %>% select(A, celltype_A = celltype), by = "A") %>%
    left_join(pd_2 %>% mutate(B = cell_id) %>% select(B, celltype_B = celltype), by = "B")

### save this result because it will be used for step5_Astrocytes as well.
saveRDS(dat, paste0(work_path, example_i, ".MNN_pairs.rds"))

dat_num = dat %>% group_by(celltype_A, celltype_B) %>% tally() %>%
    dcast(celltype_A~celltype_B)
rownames(dat_num) = as.vector(dat_num[,1])
dat_num = dat_num[,-1]
dat_num[is.na(dat_num)] = 0

dat_num_rowSum = apply(dat_num, 1, sum)

include_rows = names(dat_num_rowSum[dat_num_rowSum >= 50])

dat_num_sub = dat_num[rownames(dat_num) %in% include_rows,]

row_names_order = c("Cajal-Retzius cells",
                    "Intermediate neuronal progenitors",
                    "GABAergic cortical interneurons",
                    "Thalamic neuronal precursors",
                    "Striatal projection neurons",
                    
                    "Di/mesencephalon GABAergic neurons",
                    "Di/mesencephalon glutamatergic neurons",
                    "Hypothalamic Sim1 neurons",
                    "Midbrain dopaminergic neurons",
                    "Choroid plexus",
                    
                    "Spinal dI1 interneurons",
                    "Neurons (Slc17a8+)",
                    "Cranial motor neurons",
                    "Spinal V3 interneurons",
                    "Spinal dI2 interneurons",
                    
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons",
                    "Spinal V0 interneurons",
                    "Spinal V1 interneurons",
                    
                    "Spinal V2a interneurons",
                    "Spinal V2b interneurons",
                    "Spinal cord motor neurons",
                    "Neural progenitor cells (Neurod1+)",
                    "Astrocytes")

col_names_order = c("Telencephalon",
                    "Dorsal telencephalon",
                    "Hypothalamus",
                    "Diencephalon",
                    "Midbrain",
                    "Hypothalamus (Sim1+)",
                    "Anterior floor plate",
                    "Midbrain-hindbrain boundary",
                    "Anterior roof plate",
                    "Hindbrain",
                    "Posterior floor plate",
                    "Spinal cord/r7/r8",
                    "Floorplate and p3 domain")

### Fig. 4g
### The number of mutual nearest neighbors (MNN) pairs between pairwise neuroectodermal territories (column) and their derivative cell types (row).

heatmap.2(as.matrix(dat_num_sub[row_names_order, col_names_order]), 
          col=magma,
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(8,8))


### mapping to the backbone UMAP

pd_back = readRDS(paste0(work_path, "Neuroectoderm_backbone_adata_scale.obs.rds"))
rownames(pd_back) = as.vector(pd_back$cell_id)

for(i in row_names_order){
    print(i)
    dat_sub = dat %>% filter(celltype_A == i) %>% group_by(B) %>% tally() %>% rename(cell_id = B, freq = n)
    df = pd_back %>% select(UMAP_1 = UMAP_2d_1, UMAP_2 = UMAP_2d_2, cell_id, day) %>% left_join(dat_sub, by = "cell_id") 
    df$freq[is.na(df$freq)] = 0
    
    name = gsub("/", "_", i)
    name = gsub(" ", "_", name)
    
    try(ggplot() +
            geom_point(data = df[sample(1:nrow(df),100000),], aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey80") +
            geom_point(data = df[df$freq != 0,], aes(x = UMAP_1, y = UMAP_2, color = freq), size=0.2) +
            theme_void() +
            scale_color_viridis() +
            theme(legend.position="none") + 
            ggsave(paste0(work_path, name, ".png"), width = 8, height = 6, dpi = 300), silent = T)
}

### assigning each cell in the backbone with one potential derivatives based on MNNs

dat_m = dat %>% filter(celltype_A %in% row_names_order) %>% group_by(B, celltype_A) %>% tally() %>%
    group_by(B) %>% slice_max(order_by = n, n = 1, with_ties = F) %>%
    rename(cell_id = B, celltype_derive = celltype_A, freq = n)

df = pd_back %>% select(UMAP_1 = UMAP_2d_1, UMAP_2 = UMAP_2d_2, cell_id, day) %>% left_join(dat_m, by = "cell_id") 
df$celltype_derive[is.na(df$celltype_derive)] = "other"
df$freq[is.na(df$freq)] = 0

### Fig. 4h
### The same UMAP as the patterned neuroectoderm, but with inferred progenitor cells colored by derivative cell type with the most frequent MNN pairs.

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey90") +
        geom_point(data = df[df$celltype_derive != "other",][sample(1:sum(df$celltype_derive != "other")),], aes(x = UMAP_1, y = UMAP_2, color = celltype_derive), size=0.3) +
        theme_void() +
        scale_color_manual(values=neuroectoderm_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0(work_path, "backbone_mapping.png"), width = 8, height = 6, dpi = 300), silent = T)

### For each neuronal subtype in Fig. 4g-h,  we selected the annotation in the patterned neuroectoderm to which the most inferred progenitors had been assigned, and plotted the distribution of timepoints for that subset of inferred progenitors.

dat = datB %>% left_join(datA, by = c("A" = "A", "B" = "B")) %>%
    filter(!is.na(value.y)) %>% select(A, B) %>%
    left_join(pd_1 %>% mutate(A = cell_id) %>% select(A, celltype_A = celltype, day_A = day_value), by = "A") %>%
    left_join(pd_2 %>% mutate(B = cell_id) %>% select(B, celltype_B = celltype, day_B = day_value), by = "B") %>%
    filter(celltype_A %in% row_names_order, celltype_B %in% col_names_order)

dat_x = dat %>% group_by(celltype_A, celltype_B) %>% tally() %>% 
    group_by(celltype_A) %>% slice_max(order_by = n, n = 1)
dat_x$celltype_A = factor(dat_x$celltype_A, levels=  row_names_order)
dat_x$celltype_B = factor(dat_x$celltype_B, levels=  col_names_order)
dat_x = dat_x %>% arrange(celltype_A, celltype_B) %>% as.data.frame()

label_list = c(paste0(c(0, 2:12, 14:18, 20:34), " somites"),
               "E10.25", "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
               "E12.25", "E12.5", "E12.75")

df = NULL
for(i in 1:nrow(dat_x)){
    print(i)
    
    dat_i = dat %>% filter(celltype_A %in% dat_x$celltype_A[i], celltype_B %in% dat_x$celltype_B[i]) %>%
        select(celltype_A, day_B)
    
    df = rbind(df, dat_i)
}

df$celltype_A = factor(df$celltype_A, levels = rev(row_names_order))

### Extended Data Fig. 10n

p = ggplot(df, aes(x = day_B, y = celltype_A, fill = celltype_A)) +
    geom_density_ridges() +
    scale_fill_manual(values=neuroectoderm_color_plate) +
    theme_ridges() + 
    labs(x = "", y = "") +
    scale_x_continuous(limits = c(1,43), breaks=c(1:43), labels = label_list) +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1, family = "Helvetica"), axis.text.y = element_text(color="black", family = "Helvetica")) +
    theme(legend.position = "none")


