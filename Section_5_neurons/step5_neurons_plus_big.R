
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons_plus_big"

pd = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)
pd$celltype_update_backup = pd$celltype_update

pd_x = readRDS(paste0(work_path, "/Neurons/Neurons_adata_scale.obs.rds"))
rownames(pd_x) = as.vector(pd_x$cell_id)

pd_1 = pd[pd$cell_id %in% pd_x$cell_id,]
pd_1 = pd_1[rownames(pd_x),]
pd_1$celltype_update = as.vector(pd_x$celltype_sub_clustering)

pd_2 = pd[pd$major_trajectory == "Intermediate_progenitors",]
pd_2$celltype_update = "Intermediate progenitors"

pd_3 = pd[pd$major_trajectory == "Ependymal_cells",]
pd_3$celltype_update = "Choroid plexus"

pd_4 = pd[pd$major_trajectory == "Oligodendrocytes",]
pd_4$celltype_update = "Oligodendrocytes"

pd_5 = pd[!pd$cell_id %in% pd_x$cell_id & !pd$major_trajectory %in% c("Intermediate_progenitors", "Ependymal_cells", "Oligodendrocytes"),]

pd_y = rbind(pd_1, pd_2, pd_3, pd_4, pd_5)
pd_y = pd_y[rownames(pd),]
pd = pd_y

saveRDS(pd, paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

emb = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.PCs.csv"), header=F, as.is=T)
rownames(emb) = rownames(pd)
colnames(emb) = paste0("PC_", 1:30)
emb = as.matrix(emb)

saveRDS(emb, paste0(work_path, "/Neurons/", example_i, "_adata_scale.PCs.rds"))

pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])
day_list_2 = names(table(pd$day))

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_2))
names(day_color_plate_2) = day_list_2


fig = plot_ly(pd[sample(1:nrow(pd), 250000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd[sample(1:nrow(pd), 250000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day, colors = day_color_plate_2) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_day.html"), selfcontained = FALSE, libdir = "tmp")


fig = plot_ly(pd[sample(1:nrow(pd), 250000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~major_trajectory, colors = major_trajectory_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_major_trajectory.html"), selfcontained = FALSE, libdir = "tmp")


############
### Based on Knn, find the nearest neighboring cells for individual neurons at E8.5

example_i = "Neurons_plus_big"

pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))
emb = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.PCs.rds"))

### merge some cell types
celltype = as.vector(pd$celltype_update)

pd$celltype = as.vector(celltype)

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
                  "Mesencephalon",
                  "Posterior floor plate",
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
                  "Intermediate progenitors",
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

### only extracting early stage cells for individual neurons

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

### mutually nearest neighbors

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

saveRDS(list(pd_1 = pd_1, pd_2 = pd_2, 
             nn_matrix = nn_matrix, nn_matrix_2 = nn_matrix_2), paste0(work_path, "/Neurons/", example_i, ".MNN.rds"))

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

saveRDS(dat, paste0(work_path, "/Neurons/", example_i, ".MNN_pairs.rds"))

library(reshape2)
dat_num = dat %>% group_by(celltype_A, celltype_B) %>% tally() %>%
    dcast(celltype_A~celltype_B)
rownames(dat_num) = as.vector(dat_num[,1])
dat_num = dat_num[,-1]
dat_num[is.na(dat_num)] = 0

dat_num_rowSum = apply(dat_num, 1, sum)

include_rows = names(dat_num_rowSum[dat_num_rowSum >= 50])

dat_num_sub = dat_num[rownames(dat_num) %in% include_rows,]

row_names_order = c("Cajal-Retzius cells",
                    "Intermediate progenitors",
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
                    "Mesencephalon",
                    "Hypothalamus (Sim1+)",
                    "Anterior floor plate",
                    "Midbrain-hindbrain boundary",
                    "Anterior roof plate",
                    "Hindbrain",
                    "Posterior floor plate",
                    "Spinal cord/r7/r8",
                    "Posterior roof plate")

library("gplots")
library(viridis) 
pdf("~/share/Neuron_Brain_MNN_big.pdf", 8, 5)
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
dev.off()

### mapping to the backbone UMAP
pd_back = readRDS(paste0(work_path, "/Neurons/Neurons_plus_big_backbone_adata_scale.obs.rds"))
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
            ggsave(paste0("~/share/", name, ".png"), width = 8, height = 6, dpi = 300), silent = T)
    
}

### assigning each cell in the backbone with one potential derivatives based on MNNs

dat_m = dat %>% filter(celltype_A %in% row_names_order) %>% group_by(B, celltype_A) %>% tally() %>%
    group_by(B) %>% slice_max(order_by = n, n = 1, with_ties = F) %>%
    rename(cell_id = B, celltype_derive = celltype_A, freq = n)

df = pd_back %>% select(UMAP_1 = UMAP_2d_1, UMAP_2 = UMAP_2d_2, cell_id, day) %>% left_join(dat_m, by = "cell_id") 
df$celltype_derive[is.na(df$celltype_derive)] = "other"
df$freq[is.na(df$freq)] = 0

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey90") +
        geom_point(data = df[df$celltype_derive != "other",][sample(1:sum(df$celltype_derive != "other")),], aes(x = UMAP_1, y = UMAP_2, color = celltype_derive), size=0.3) +
        theme_void() +
        scale_color_manual(values=neuron_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0("~/share/backbone_mapping.png"), width = 8, height = 6, dpi = 300), silent = T)


### making the scatter plot used for timing

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

for(i in 1:nrow(dat_x)){
    print(i)
    
    name = paste0(dat_x$celltype_A[i], "_", dat_x$celltype_B[i])
    name = gsub("/", "_", name)
    name = gsub(" ", "_", name)
    
    try(dat %>% filter(celltype_A %in% dat_x$celltype_A[i], celltype_B %in% dat_x$celltype_B[i]) %>% 
            group_by(day_A, day_B) %>% tally() %>% rename(Frequency = n) %>%
            ggplot() +
            geom_point(aes(x = day_B, y = day_A, color = Frequency, size = Frequency)) +
            scale_color_viridis(option = "A") +
            labs(x = dat_x$celltype_B[i], y = dat_x$celltype_A[i]) +
            theme_classic(base_size = 10) +
            scale_x_continuous(limits = c(1,43), breaks=c(1:43), labels = label_list) +
            scale_y_continuous(limits = c(1,43), breaks=c(1:43), labels = label_list) +
            theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1, family = "Helvetica"), axis.text.y = element_text(color="black", family = "Helvetica")) +
            ggsave(paste0("~/share/", name, ".png"), width = 6, height = 5, dpi = 300), silent = T)
}

df = NULL
for(i in 1:nrow(dat_x)){
    print(i)

    dat_i = dat %>% filter(celltype_A %in% dat_x$celltype_A[i], celltype_B %in% dat_x$celltype_B[i]) %>%
        select(celltype_A, day_B)
    
    df = rbind(df, dat_i)
}

df$celltype_A = factor(df$celltype_A, levels = rev(row_names_order))

library(ggridges)
p=ggplot(df, aes(x = day_B, y = celltype_A, fill = celltype_A)) +
    geom_density_ridges() +
    scale_fill_manual(values=neuron_color_plate) +
    theme_ridges() + 
    labs(x = "", y = "") +
    scale_x_continuous(limits = c(1,43), breaks=c(1:43), labels = label_list) +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1, family = "Helvetica"), axis.text.y = element_text(color="black", family = "Helvetica")) +
    theme(legend.position = "none")
pdf("~/share/Mapping_timing.pdf",10,10)
p
dev.off()

df_mean = df %>% group_by(celltype_A) %>% summarise(row_mean = mean(day_B))
df_sd = df %>% group_by(celltype_A) %>% summarise(row_sd = sd(day_B))
df = dat_x %>% left_join(df_mean, by = "celltype_A") %>% left_join(df_sd, by = "celltype_A") %>%
    mutate(low_cutoff = row_mean - 2*row_sd, up_cutoff = row_mean + 2*row_sd)
df$low_day = label_list[round(df$low_cutoff)]
df$high_day = label_list[round(df$up_cutoff)]
df$low_day[is.na(df$low_day)] = label_list[1]
df$high_day[is.na(df$high_day)] = label_list[43]
df$mean_day = label_list[round(df$row_mean)]

df_sub = df[,c("celltype_A", "low_day", "mean_day", "high_day")]

### Three cell types are removed in the heatmap, due to <50 MNN pairs are detected.
### Cerebellar Purkinje cells, Precerebellar neurons, Spinal dI6 interneurons

target_celltype = c("Cerebellar Purkinje cells", "Precerebellar neurons", "Spinal dI6 interneurons")

pd_1 = pd_integration[pd_integration$celltype %in% target_celltype,]
pd_2 = pd_integration[pd_integration$celltype %in% row_names_order,]

pd_1_x = pd_1 %>% group_by(celltype) %>% slice_min(order_by = day_value, n = 500, with_ties = F)
pd_1_y = pd_1_x %>% group_by(celltype) %>% tally() %>% filter(n >= 500) %>% pull(celltype)
pd_1_x = pd_1_x %>% filter(celltype %in% pd_1_y)

pd_1 = pd_integration[pd_integration$cell_id %in% pd_1_x$cell_id,]

pd_2_x = pd_2
pd_2 = pd_integration[pd_integration$cell_id %in% pd_2_x$cell_id,]

emb_1 = PC_integration[pd_integration$cell_id %in% pd_1_x$cell_id,]
emb_2 = PC_integration[pd_integration$cell_id %in% pd_2_x$cell_id,]

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

### mutually nearest neighbors

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

saveRDS(list(pd_1 = pd_1, pd_2 = pd_2, 
             nn_matrix = nn_matrix, nn_matrix_2 = nn_matrix_2), paste0(work_path, "/Neurons/", example_i, ".MNN_iter.rds"))

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

library(reshape2)
dat_num = dat %>% group_by(celltype_A, celltype_B) %>% tally() %>%
    dcast(celltype_A~celltype_B)
rownames(dat_num) = as.vector(dat_num[,1])
dat_num = dat_num[,-1]
dat_num[is.na(dat_num)] = 0

dat_num_rowSum = apply(dat_num, 1, sum)

include_rows = names(dat_num_rowSum[dat_num_rowSum >= 50])

dat_num_sub = dat_num[rownames(dat_num) %in% include_rows,]

dat_num_sub_x = matrix(0, nrow(dat_num_sub), 2)
rownames(dat_num_sub_x) = rownames(dat_num_sub)
colnames(dat_num_sub_x) = c("Choroid plexus", "Astrocytes")
dat_num_sub = cbind(dat_num_sub, dat_num_sub_x)

col_names_order = c("Cajal-Retzius cells",
                    "Intermediate progenitors",
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

row_names_order = c("Cerebellar Purkinje cells", "Precerebellar neurons", "Spinal dI6 interneurons")

library("gplots")
library(viridis) 
pdf("~/share/Neuron_Brain_MNN_big_iter.pdf", 8, 5)
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
dev.off()




########################################
### Making 2D UMAP on backbone cells ###
########################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons_plus_big_backbone"

pd = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)


p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.3) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Neurons/", example_i, ".2D_UMAP.day.png"), width = 6, height = 6, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.2) +
    theme_void() +
    scale_color_manual(values=neuron_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Neurons/", example_i, ".2D_UMAP.celltype_update.png"), width = 8, height = 6, dpi = 300)


fig = plot_ly(pd[sample(1:nrow(pd), 250000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update, colors = neuron_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

saveRDS(pd, paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

### create the monocle3 cds file
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons_plus_big_backbone"

pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

set.seed(1234)
pd_sub = pd[sample(1:nrow(pd), 100000),]

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = names(table(as.vector(pd_sub$embryo_id)))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE])
}

gene_count = gene_count[,rownames(pd_sub)]
obj = CreateSeuratObject(gene_count, meta.data = pd_sub)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Neurons/", example_i, ".cds.rds"))

### New supplementary figure with Territory marker UMAPs?
cds_sub = cds[,sample(1:ncol(cds), 50000)]

my_plot_cells(cds_sub, genes = c("Emx1","Emx2","Otx2","Pax6"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_1.png", dpi = 300, height = 2, width = 11)

my_plot_cells(cds_sub, genes = c("Wnt3a","Bmp4","Nkx2-1","Rax"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_2.png", dpi = 300, height = 2, width = 11)

my_plot_cells(cds_sub, genes = c("Dmbx1","Sim1","Shh","Foxa2"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_3.png", dpi = 300, height = 2, width = 11)

my_plot_cells(cds_sub, genes = c("En1","Pax2","Pax5","Fgf8"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_4.png", dpi = 300, height = 2, width = 11)

my_plot_cells(cds_sub, genes = c("Lmx1a","Msx1","Mafb","Egr2"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_5.png", dpi = 300, height = 2, width = 11)

my_plot_cells(cds_sub, genes = c("Hoxa2","Hoxb4","Hoxd4","Hoxc6"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_6.png", dpi = 300, height = 2, width = 11)


### backup - checking expression of TFs of 11 spinal interneurons in the backbone UMAP
my_plot_cells(cds, genes = c("Lhx9", "Lhx2", "Barhl1", "Mecom", "Pax3", "Prdm16"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_1.png", dpi = 300, height = 2, width = 14.4)

my_plot_cells(cds, genes = c("Isl1", "Otp", "Bnc2", "Tfap2b", "Pax2", "Pknox2"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_2.png", dpi = 300, height = 2, width = 14.4)

my_plot_cells(cds, genes = c("Lmx1b", "Ebf2", "Pax8", "Esrrg", "Onecut2", "Evx1"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_3.png", dpi = 300, height = 2, width = 14.4)

my_plot_cells(cds, genes = c("Foxp2", "Tfdp2", "Lhx4", "Vsx2", "Shox2"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_4.png", dpi = 300, height = 2, width = 12)

my_plot_cells(cds, genes = c("Zfpm2", "Gata3", "Tal1", "Sim1", "Nkx2-2", "Arid5b"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_5.png", dpi = 300, height = 2, width = 14.4)





###########################################################################
### subclustering Cerebellar Purkinje cells and Spinal dI2 interneurons ###
###########################################################################

example_i = "Neurons_plus_big"

pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

pd_sub = pd[pd$celltype_update %in% c("Cerebellar Purkinje cells", "Spinal dI2 interneurons"),]

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = names(table(as.vector(pd_sub$embryo_id)))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE])
}

gene_count = gene_count[,rownames(pd_sub)]

pd_sub = pd_sub[,!colnames(pd_sub) %in% c("leiden","leiden_res_1","leiden_res_5","leiden_res_20",
                                          "UMAP_1","UMAP_2","UMAP_3","UMAP_2d_1","UMAP_2d_2","celltype_update_backup")]

dat = readRDS(paste0(work_path, "/Neurons/", example_i, ".MNN_iter.rds"))
pd_1 = dat[["pd_1"]]
pd_2 = dat[["pd_2"]]
nn_matrix = dat[["nn_matrix"]]
nn_matrix_2 = dat[["nn_matrix_2"]]

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

cell_group = as.vector(pd_sub$celltype_update)
cell_group[pd_sub$cell_id %in% as.vector(pd_1$cell_id) & 
               pd_sub$celltype_update == "Cerebellar Purkinje cells"] = "Cerebellar Purkinje cells early 500"
cell_group[pd_sub$cell_id %in% as.vector(dat$A)] = "Cerebellar Purkinje cells MNNs"
cell_group[pd_sub$cell_id %in% as.vector(dat$B)] = "Spinal dI2 interneurons MNNs"

pd_sub$cell_group = as.vector(cell_group)

name = "Cerebellar_Purkinje"
Matrix::writeMM(t(gene_count), paste0(work_path, "/Neurons/", name, ".mtx"))
write.csv(pd_sub, paste0(work_path, "/Neurons/", name, ".df_cell.csv"))

pd_sub = read.csv(paste0(work_path, "/Neurons/", name, "_adata_scale.obs.csv"), header=T, row.names=1)

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.6) +
    theme_void() +
    #scale_color_manual(values=neuron_color_plate) +
    #theme(legend.position="none") + 
    scale_color_brewer(palette = "Set1") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    ggsave(paste0(work_path, "/Neurons/", name, ".celltype_update.png"), width = 7, height = 6, dpi = 300)


p = ggplot() +
    geom_point(data = pd_sub, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(data = pd_sub, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey90", size=0.8) +
    geom_point(data = pd_sub %>% filter(cell_group %in% c("Cerebellar Purkinje cells early 500")), aes(x = UMAP_2d_1, y = UMAP_2d_2, color = cell_group), size=1) +
    geom_point(data = pd_sub %>% filter(cell_group %in% c("Cerebellar Purkinje cells MNNs","Spinal dI2 interneurons MNNs")), aes(x = UMAP_2d_1, y = UMAP_2d_2, color = cell_group), size=1) +
    theme_void() +
    #scale_color_manual(values=neuron_color_plate) +
    #theme(legend.position="none") + 
    scale_color_brewer(palette = "Dark2") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    ggsave(paste0(work_path, "/Neurons/", name, ".cell_group.png"), width = 7, height = 6, dpi = 300)

day_list_sub = day_list[day_list %in% pd_sub$day]
pd_sub$day = factor(pd_sub$day, levels = day_list_sub)

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_sub))
names(day_color_plate_2) = day_list_sub

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.6) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Neurons/", name, ".day.png"), width = 7, height = 6, dpi = 300)


obj = CreateSeuratObject(gene_count, meta.data = pd_sub)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Neurons/", name, ".cds.rds"))

my_plot_cells(cds, genes = c("Foxd3","Lhx1","Lhx5","Pou4f1"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_6.png", dpi = 300, height = 2, width = 8)

my_plot_cells(cds, genes = c("Gad1","Gad2","Skor2"), how_many_rows = 1, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_7.png", dpi = 300, height = 2, width = 6)

