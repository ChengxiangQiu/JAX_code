
################################################
### checking the timeline for mesoderm cells ###
################################################


source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

### normalization by the total number of cells from each stage ####
pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd = pd_all[pd_all$day %in% c("E8.5","E8.75","E9.0","E9.25","E9.5","E9.75","E10.0"),]

pd_LPM = readRDS(paste0(work_path, "/LPM/LPM_adata_scale.obs.rds"))
pd_FM = readRDS(paste0(work_path, "/mtx/example/Facial_mesenchyme_adata_scale.obs.rds"))

pd_1 = pd[pd$cell_id %in% as.vector(pd_LPM$cell_id),]
pd_2 = pd[pd$cell_id %in% as.vector(pd_FM$cell_id),]
pd_3 = pd[!pd$cell_id %in% as.vector(pd_LPM$cell_id) & !pd$cell_id %in% as.vector(pd_FM$cell_id),]

pd_1_x = pd_1 %>% left_join(pd_LPM[,c("cell_id","celltype_sub_clustering")], by = "cell_id")
pd_2_x = pd_2 %>% left_join(pd_FM[,c("cell_id","celltype_sub_clustering")], by = "cell_id")

pd_1$celltype_update = paste0("LPM:",as.vector(pd_1_x$celltype_sub_clustering))
pd_2$celltype_update = paste0("FM:",as.vector(pd_2_x$celltype_sub_clustering))
pd_x = rbind(pd_1, pd_2, pd_3)
pd_x = pd_x[rownames(pd),]
pd = pd_x

celltype_include = c("LPM:Allantois",
                     "LPM:Amniotic mesoderm",
                     "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                     "LPM:Extraembryonic mesoderm",
                     "LPM:Gut mesenchyme",
                     "LPM:Renal capsule",
                     "LPM:Somatic mesoderm",
                     "LPM:Splanchnic mesoderm",
                     "Mesodermal progenitors (Tbx6+)",
                     "First heart field",
                     "Second heart field",
                     "Neuromesodermal progenitors",
                     
                     "Anterior intermediate mesoderm",
                     "LPM:Foregut mesenchyme",
                     "LPM:Hepatic mesenchyme",
                     "LPM:Proepicardium (Tbx18+)",
                     "Posterior intermediate mesoderm",
                     "Limb mesenchyme progenitors",
                     
                     "LPM:Meninges",
                     "LPM:Vascular smooth muscle cells",
                     "LPM:Mesothelial cells",
                     "LPM:Gonad progenitor cells",
                     "LPM:Lung mesenchyme")

cell_num = readRDS(paste0(work_path, "/embryo/cell_num_prediction.rds"))

pd$somite_count = factor(pd$somite_count, levels = somite_list)
pd$celltype_udpate = as.vector(pd$celltype_update)

x = pd %>% group_by(somite_count, day, celltype_update) %>% tally() %>% 
    left_join(pd_all %>% group_by(somite_count) %>% tally() %>% rename(total_n = n), by = "somite_count") %>%
    left_join(cell_num %>% select(day, cell_num_pred), by = "day") %>%
    mutate(estimated_num_log2 = log2(ceiling(cell_num_pred*n/total_n))) %>%
    filter(celltype_update %in% celltype_include)

x$somite_count = factor(x$somite_count, levels = somite_list)
x$celltype_update = factor(x$celltype_update, levels = celltype_include)

x$somite_count_value = as.vector(x$somite_count)
x$somite_count_value = as.numeric(gsub(" somites","",x$somite_count_value))

x_first_somite_count = x %>% group_by(celltype_update) %>% filter(n >= 10) %>% 
    slice_min(order_by = somite_count_value, n = 1) %>% mutate(first_somite_count_value = somite_count_value)

x_sub = x %>% left_join(x_first_somite_count %>% select(celltype_update, first_somite_count_value), by = "celltype_update") %>%
    filter(somite_count_value >= first_somite_count_value)


p = x_sub %>% 
    ggplot(aes(x=somite_count, y=estimated_num_log2, fill = celltype_update)) + 
    geom_bar(stat='identity') + facet_grid(celltype_update ~ .) + 
#    geom_bar(stat='identity') + facet_grid(celltype_update ~ ., switch = "y") + 
#    theme(strip.text.y.left = element_text(angle = 0)) +
    labs(x='',y='Log2(Estimated # of cells)') +
    scale_fill_manual(values=LPM_E85_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0("~/share/", "LPM_pre_E10_density.pdf"), 7, 12)
p
dev.off()






################################
### Focusing on somite 5-20 ####
################################




source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM_somite_5_20"

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd = pd_all[pd_all$somite_count %in% paste0(c(5:20), " somites"),]

pd_LPM = readRDS(paste0(work_path, "/LPM/LPM_adata_scale.obs.rds"))
pd_FM = readRDS(paste0(work_path, "/mtx/example/Facial_mesenchyme_adata_scale.obs.rds"))

pd_1 = pd[pd$cell_id %in% as.vector(pd_LPM$cell_id),]
pd_2 = pd[pd$cell_id %in% as.vector(pd_FM$cell_id),]
pd_3 = pd[!pd$cell_id %in% as.vector(pd_LPM$cell_id) & !pd$cell_id %in% as.vector(pd_FM$cell_id),]

pd_1_x = pd_1 %>% left_join(pd_LPM[,c("cell_id","celltype_sub_clustering")], by = "cell_id")
pd_2_x = pd_2 %>% left_join(pd_FM[,c("cell_id","celltype_sub_clustering")], by = "cell_id")

pd_1$celltype_update = paste0("LPM:",as.vector(pd_1_x$celltype_sub_clustering))
pd_2$celltype_update = paste0("FM:",as.vector(pd_2_x$celltype_sub_clustering))
pd_x = rbind(pd_1, pd_2, pd_3)
pd_x = pd_x[rownames(pd),]

### first backbone
pd_sub = pd_x[pd_x$celltype_update %in% c("LPM:Allantois",
                                          "LPM:Amniotic mesoderm",
                                          "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                                          "LPM:Extraembryonic mesoderm",
                                          "LPM:Gut mesenchyme",
                                          "LPM:Renal capsule",
                                          "LPM:Somatic mesoderm",
                                          "LPM:Splanchnic mesoderm",
                                          "Mesodermal progenitors (Tbx6+)",
                                          "First heart field",
                                          "Second heart field",
                                          "Neuromesodermal progenitors"),]
rownames(pd_sub) = as.vector(pd_sub$cell_id)

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

writeMM(t(gene_count), paste0(work_path, "/LPM/", example_i, "_backbone.gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/LPM/", example_i, "_backbone.df_cell.csv"))

### then, including both backbone and derivatives

pd_sub = pd_x[pd_x$celltype_update %in% c("LPM:Allantois",
                                          "LPM:Amniotic mesoderm",
                                          "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                                          "LPM:Extraembryonic mesoderm",
                                          "LPM:Gut mesenchyme",
                                          "LPM:Renal capsule",
                                          "LPM:Somatic mesoderm",
                                          "LPM:Splanchnic mesoderm",
                                          "Mesodermal progenitors (Tbx6+)",
                                          "First heart field",
                                          "Second heart field",
                                          "Neuromesodermal progenitors",
                                          
                                          "Anterior intermediate mesoderm",
                                          "LPM:Foregut mesenchyme",
                                          "LPM:Hepatic mesenchyme",
                                          "LPM:Proepicardium (Tbx18+)",
                                          "Posterior intermediate mesoderm",
                                          "Limb mesenchyme progenitors"),]
rownames(pd_sub) = as.vector(pd_sub$cell_id)

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

writeMM(t(gene_count), paste0(work_path, "/LPM/", example_i, "_all.gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/LPM/", example_i, "_all.df_cell.csv"))


### After running Scanpy to get embedding, plotting the backbone UMAP

example_i = "LPM_somite_5_20"

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_backbone_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.6) +
    theme_void() +
    scale_color_manual(values=LPM_E85_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_update.png"), width = 8, height = 6, dpi = 300)


### mapping derivatives to the original backbone based on MNN

example_i = "LPM_somite_5_20"

pd_backbone = read.csv(paste0(work_path, "/LPM/", example_i, "_backbone_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd_backbone) = as.vector(pd_backbone$cell_id)

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_all_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)
emb = read.csv(paste0(work_path, "/LPM/", example_i, "_all_adata_scale.PCs.csv"), header=F)
rownames(emb) = as.vector(pd$cell_id)
colnames(emb) = paste0("PC_", 1:30)
emb = as.matrix(emb)

pd$celltype = as.vector(pd$celltype_update)
pd_integration = pd
PC_integration = emb

pd_integration$day_value = as.numeric(gsub(" somites", "", as.vector(pd_integration$somite_count)))

### cell type of each column (backbone)
col_cell_type = c("LPM:Allantois",
                  "LPM:Amniotic mesoderm",
                  "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                  "LPM:Extraembryonic mesoderm",
                  "LPM:Gut mesenchyme",
                  "LPM:Renal capsule",
                  "LPM:Somatic mesoderm",
                  "LPM:Splanchnic mesoderm",
                  "Mesodermal progenitors (Tbx6+)",
                  "First heart field",
                  "Second heart field",
                  "Neuromesodermal progenitors")

### cell type of each row (derivatives)
row_cell_type = c("Anterior intermediate mesoderm",
                  "LPM:Foregut mesenchyme",
                  "LPM:Hepatic mesenchyme",
                  "LPM:Proepicardium (Tbx18+)",
                  "Posterior intermediate mesoderm",
                  "Limb mesenchyme progenitors")

### only extracting early stage cells for individual neurons

pd_1 = pd_integration[pd_integration$celltype %in% row_cell_type,]
pd_2 = pd_integration[pd_integration$celltype %in% col_cell_type,]

pd_1_x = pd_1 %>% group_by(celltype) %>% slice_min(order_by = day_value, n = 500, with_ties = F)
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

saveRDS(dat, paste0(work_path, "/LPM/", example_i, ".MNN_pairs.rds"))

library(reshape2)
dat_num = dat %>% group_by(celltype_A, celltype_B) %>% tally() %>%
    dcast(celltype_A~celltype_B)
rownames(dat_num) = as.vector(dat_num[,1])
dat_num = dat_num[,-1]
dat_num[is.na(dat_num)] = 0

dat_num_rowSum = apply(dat_num, 1, sum)
include_rows = names(dat_num_rowSum[dat_num_rowSum >= 50])
dat_num_sub = dat_num[rownames(dat_num) %in% include_rows,]

celltype_not_include = col_cell_type[!col_cell_type %in% colnames(dat_num_sub)]
for(i in celltype_not_include){
    new_name = c(colnames(dat_num_sub), i)
    dat_num_sub = cbind(dat_num_sub, rep(0, nrow(dat_num_sub)))
    colnames(dat_num_sub) = new_name
}

row_names_order = row_cell_type[row_cell_type %in% rownames(dat_num_sub)]
col_names_order = col_cell_type[col_cell_type %in% colnames(dat_num_sub)]

library("gplots")
library(viridis) 
pdf(paste0("~/share/", example_i, ".MNN_pairs.pdf"), 8, 5)
heatmap.2(as.matrix(dat_num_sub[row_names_order,col_names_order]), 
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

### assigning each cell in the backbone with one potential derivatives based on MNNs

dat_m = dat %>% filter(celltype_A %in% row_names_order) %>% group_by(B, celltype_A) %>% tally() %>%
    group_by(B) %>% slice_max(order_by = n, n = 1, with_ties = F) %>%
    rename(cell_id = B, celltype_derive = celltype_A, freq = n)

df = pd_backbone %>% select(UMAP_1 = UMAP_2d_1, UMAP_2 = UMAP_2d_2, cell_id, day) %>% left_join(dat_m, by = "cell_id") 
df$celltype_derive[is.na(df$celltype_derive)] = "other"
df$freq[is.na(df$freq)] = 0

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey90") +
        geom_point(data = df[df$celltype_derive != "other",][sample(1:sum(df$celltype_derive != "other")),], aes(x = UMAP_1, y = UMAP_2, color = celltype_derive), size=0.3) +
        theme_void() +
        scale_color_manual(values=LPM_E85_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0("~/share/", example_i, ".backbone_mapping.png"), width = 8, height = 6, dpi = 300), silent = T)





#################
### Re Jay's comment, Can we consider providing one anecdote here? For example, 
### what is the DE for the precursors of hepatic vs. foregut mesenchyme (and how 
### does that differ from the DE across the MNN jumps (e.g. precursors -> hepatic mesenchyme)? 
### Maybe we won't do anything here but can you take quick look? I worry reader 
### will want a little more insight than just the mappings..


dat = readRDS(paste0(work_path, "/LPM/", example_i, ".MNN_pairs.rds"))

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_all_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = names(table(as.vector(pd$embryo_id)))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd$cell_id, drop=FALSE])
}
gene_count = gene_count[,as.vector(pd$cell_id)]

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

### identifying prog. for Hepatic and Foregut based on MNN pairs (dat)

x1 = as.vector(unique(dat$B[dat$celltype_A == "LPM:Hepatic mesenchyme"]))
x2 = as.vector(unique(dat$A[dat$celltype_A == "LPM:Hepatic mesenchyme"]))
y1 = as.vector(unique(dat$B[dat$celltype_A == "LPM:Foregut mesenchyme"]))
y2 = as.vector(unique(dat$A[dat$celltype_A == "LPM:Foregut mesenchyme"]))

x1_y1 = x1[x1 %in% y1]
x1 = x1[!x1 %in% x1_y1]
y1 = y1[!y1 %in% x1_y1]

group = rep("other", ncol(obj))
group[obj$cell_id %in% x1] = "Hepatic_prog"
group[obj$cell_id %in% x2] = "Hepatic_itself"
group[obj$cell_id %in% y1] = "Foregut_prog"
group[obj$cell_id %in% y2] = "Foregut_itself"

obj$group = as.vector(group)


### first, comparing precursors between Hepatic_prog vs. Foregut_prog

obj_sub = subset(obj, subset = group %in% c("Hepatic_prog", "Foregut_prog"))
Idents(obj_sub) = as.vector(obj_sub$group)

res = FindMarkers(obj_sub, ident.1 = "Hepatic_prog", ident.2 = "Foregut_prog")
res = res %>% mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    filter(p_val_adj < 0.05)
res$high_in_which = if_else(res$avg_logFC > 0, "Hepatic_prog", "Foregut_prog")

res_1 = res

### second, comparing Hepatic_prog vs. Hepatic_itself

obj_sub = subset(obj, subset = group %in% c("Hepatic_prog", "Hepatic_itself"))
Idents(obj_sub) = as.vector(obj_sub$group)

res = FindMarkers(obj_sub, ident.1 = "Hepatic_prog", ident.2 = "Hepatic_itself")
res = res %>% mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    filter(p_val_adj < 0.05)
res$high_in_which = if_else(res$avg_logFC > 0, "Hepatic_prog", "Hepatic_itself")

res_2 = res


### third, comparing Foregut_prog vs. Foregut_itself

obj_sub = subset(obj, subset = group %in% c("Foregut_prog", "Foregut_itself"))
Idents(obj_sub) = as.vector(obj_sub$group)

res = FindMarkers(obj_sub, ident.1 = "Foregut_prog", ident.2 = "Foregut_itself")
res = res %>% mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    filter(p_val_adj < 0.05)
res$high_in_which = if_else(res$avg_logFC > 0, "Foregut_prog", "Foregut_itself")

res_3 = res

### if any overlap

x1 = as.vector(res_1$gene_short_name[res_1$high_in_which == "Hepatic_prog"])
x2 = as.vector(res_2$gene_short_name[res_2$high_in_which == "Hepatic_prog"])

y1 = as.vector(res_1$gene_short_name[res_1$high_in_which == "Foregut_prog"])
y2 = as.vector(res_3$gene_short_name[res_3$high_in_which == "Foregut_prog"])

write.csv(res_1, "~/share/compare_1_Hepatic_prog_vs_Foregut_prog.csv")
write.csv(res_2, "~/share/compare_2_Hepatic_prog_vs_Hepatic_itself.csv")
write.csv(res_3, "~/share/compare_3_Foregut_prog_vs_Foregut_itself.csv")














#################################
### Focusing on somite 26-34 ####
#################################




source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM_somite_26_34"

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd = pd_all[pd_all$somite_count %in% paste0(c(26:34), " somites"),]

pd_LPM = readRDS(paste0(work_path, "/LPM/LPM_adata_scale.obs.rds"))
pd_FM = readRDS(paste0(work_path, "/mtx/example/Facial_mesenchyme_adata_scale.obs.rds"))

pd_1 = pd[pd$cell_id %in% as.vector(pd_LPM$cell_id),]
pd_2 = pd[pd$cell_id %in% as.vector(pd_FM$cell_id),]
pd_3 = pd[!pd$cell_id %in% as.vector(pd_LPM$cell_id) & !pd$cell_id %in% as.vector(pd_FM$cell_id),]

pd_1_x = pd_1 %>% left_join(pd_LPM[,c("cell_id","celltype_sub_clustering")], by = "cell_id")
pd_2_x = pd_2 %>% left_join(pd_FM[,c("cell_id","celltype_sub_clustering")], by = "cell_id")

pd_1$celltype_update = paste0("LPM:",as.vector(pd_1_x$celltype_sub_clustering))
pd_2$celltype_update = paste0("FM:",as.vector(pd_2_x$celltype_sub_clustering))
pd_x = rbind(pd_1, pd_2, pd_3)
pd_x = pd_x[rownames(pd),]

### first backbone
pd_sub = pd_x[pd_x$celltype_update %in% c("LPM:Allantois",
                                          "LPM:Amniotic mesoderm",
                                          "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                                          "LPM:Extraembryonic mesoderm",
                                          "LPM:Gut mesenchyme",
                                          "LPM:Renal capsule",
                                          "LPM:Somatic mesoderm",
                                          "LPM:Splanchnic mesoderm",
                                          "Mesodermal progenitors (Tbx6+)",
                                          "First heart field",
                                          "Second heart field",
                                          "Neuromesodermal progenitors",
                                          
                                          "Anterior intermediate mesoderm",
                                          "LPM:Foregut mesenchyme",
                                          "LPM:Hepatic mesenchyme",
                                          "LPM:Proepicardium (Tbx18+)",
                                          "Posterior intermediate mesoderm",
                                          "Limb mesenchyme progenitors"),]
rownames(pd_sub) = as.vector(pd_sub$cell_id)

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

writeMM(t(gene_count), paste0(work_path, "/LPM/", example_i, "_backbone.gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/LPM/", example_i, "_backbone.df_cell.csv"))

### then, including both backbone and derivatives

pd_sub = pd_x[pd_x$celltype_update %in% c("LPM:Allantois",
                                          "LPM:Amniotic mesoderm",
                                          "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                                          "LPM:Extraembryonic mesoderm",
                                          "LPM:Gut mesenchyme",
                                          "LPM:Renal capsule",
                                          "LPM:Somatic mesoderm",
                                          "LPM:Splanchnic mesoderm",
                                          "Mesodermal progenitors (Tbx6+)",
                                          "First heart field",
                                          "Second heart field",
                                          "Neuromesodermal progenitors",
                                          
                                          "Anterior intermediate mesoderm",
                                          "LPM:Foregut mesenchyme",
                                          "LPM:Hepatic mesenchyme",
                                          "LPM:Proepicardium (Tbx18+)",
                                          "Posterior intermediate mesoderm",
                                          "Limb mesenchyme progenitors",
                                          
                                          "LPM:Meninges",
                                          "LPM:Vascular smooth muscle cells",
                                          "LPM:Mesothelial cells",
                                          "LPM:Gonad progenitor cells",
                                          "LPM:Lung mesenchyme"),]
rownames(pd_sub) = as.vector(pd_sub$cell_id)

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

writeMM(t(gene_count), paste0(work_path, "/LPM/", example_i, "_all.gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/LPM/", example_i, "_all.df_cell.csv"))



### After running Scanpy to get embedding, plotting the backbone UMAP

example_i = "LPM_somite_26_34"

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_backbone_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.6) +
    theme_void() +
    scale_color_manual(values=LPM_E85_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_update.png"), width = 8, height = 6, dpi = 300)



### mapping derivatives to the original backbone based on MNN

example_i = "LPM_somite_26_34"

pd_backbone = read.csv(paste0(work_path, "/LPM/", example_i, "_backbone_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd_backbone) = as.vector(pd_backbone$cell_id)

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_all_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)
emb = read.csv(paste0(work_path, "/LPM/", example_i, "_all_adata_scale.PCs.csv"), header=F)
rownames(emb) = as.vector(pd$cell_id)
colnames(emb) = paste0("PC_", 1:30)
emb = as.matrix(emb)

pd$celltype = as.vector(pd$celltype_update)
pd_integration = pd
PC_integration = emb

pd_integration$day_value = as.numeric(gsub(" somites", "", as.vector(pd_integration$somite_count)))

### cell type of each column (backbone)
col_cell_type = c("LPM:Allantois",
                  "LPM:Amniotic mesoderm",
                  "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                  "LPM:Extraembryonic mesoderm",
                  "LPM:Gut mesenchyme",
                  "LPM:Renal capsule",
                  "LPM:Somatic mesoderm",
                  "LPM:Splanchnic mesoderm",
                  "Mesodermal progenitors (Tbx6+)",
                  "First heart field",
                  "Second heart field",
                  "Neuromesodermal progenitors",
                  
                  "Anterior intermediate mesoderm",
                  "LPM:Foregut mesenchyme",
                  "LPM:Hepatic mesenchyme",
                  "LPM:Proepicardium (Tbx18+)",
                  "Posterior intermediate mesoderm",
                  "Limb mesenchyme progenitors")

### cell type of each row (derivatives)
row_cell_type = c("LPM:Meninges",
                  "LPM:Vascular smooth muscle cells",
                  "LPM:Mesothelial cells",
                  "LPM:Gonad progenitor cells",
                  "LPM:Lung mesenchyme")

### only extracting early stage cells for individual neurons

pd_1 = pd_integration[pd_integration$celltype %in% row_cell_type,]
pd_2 = pd_integration[pd_integration$celltype %in% col_cell_type,]

pd_1_x = pd_1 %>% group_by(celltype) %>% slice_min(order_by = day_value, n = 500, with_ties = F)
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

saveRDS(dat, paste0(work_path, "/LPM/", example_i, ".MNN_pairs.rds"))

library(reshape2)
dat_num = dat %>% group_by(celltype_A, celltype_B) %>% tally() %>%
    dcast(celltype_A~celltype_B)
rownames(dat_num) = as.vector(dat_num[,1])
dat_num = dat_num[,-1]
dat_num[is.na(dat_num)] = 0

dat_num_rowSum = apply(dat_num, 1, sum)
include_rows = names(dat_num_rowSum[dat_num_rowSum >= 50])
dat_num_sub = dat_num[rownames(dat_num) %in% include_rows,]

celltype_not_include = col_cell_type[!col_cell_type %in% colnames(dat_num_sub)]
for(i in celltype_not_include){
    new_name = c(colnames(dat_num_sub), i)
    dat_num_sub = cbind(dat_num_sub, rep(0, nrow(dat_num_sub)))
    colnames(dat_num_sub) = new_name
}

row_names_order = row_cell_type[row_cell_type %in% rownames(dat_num_sub)]
col_names_order = col_cell_type[col_cell_type %in% colnames(dat_num_sub)]

library("gplots")
library(viridis) 
pdf(paste0("~/share/", example_i, ".MNN_pairs.pdf"), 8, 5)
heatmap.2(as.matrix(dat_num_sub[row_names_order,col_names_order]), 
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

### assigning each cell in the backbone with one potential derivatives based on MNNs

dat_m = dat %>% filter(celltype_A %in% row_names_order) %>% group_by(B, celltype_A) %>% tally() %>%
    group_by(B) %>% slice_max(order_by = n, n = 1, with_ties = F) %>%
    rename(cell_id = B, celltype_derive = celltype_A, freq = n)

df = pd_backbone %>% select(UMAP_1 = UMAP_2d_1, UMAP_2 = UMAP_2d_2, cell_id, day) %>% left_join(dat_m, by = "cell_id") 
df$celltype_derive[is.na(df$celltype_derive)] = "other"
df$freq[is.na(df$freq)] = 0

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey90") +
        geom_point(data = df[df$celltype_derive != "other",][sample(1:sum(df$celltype_derive != "other")),], aes(x = UMAP_1, y = UMAP_2, color = celltype_derive), size=0.6) +
        theme_void() +
        scale_color_manual(values=LPM_E85_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0("~/share/", example_i, ".backbone_mapping.png"), width = 8, height = 6, dpi = 300), silent = T)



