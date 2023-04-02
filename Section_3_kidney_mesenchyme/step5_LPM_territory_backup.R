

#############################################################################################
### Re-embedding on cells from stages E85, after including some other mesoderm derivatives ##
#############################################################################################


source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM_E85"

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

pd_x = readRDS("/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/NovaSeq/run_4_E85/pd_3D.rds")
pd_y = pd %>% left_join(pd_x[,c("cell_id","celltype_sub_clustering")], by = "cell_id")
pd$celltype_sub_clustering = as.vector(pd_y$celltype_sub_clustering)

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_sub_clustering, colors = LPM_E85_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_sub_clustering.html"), selfcontained = FALSE, libdir = "tmp")

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(12)
names(day_color_plate_2) = paste0(c(0, 2:12), " somites")

pd$somite_count = factor(pd$somite_count, levels = somite_list[somite_list %in% pd$somite_count])
fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~somite_count, colors = day_color_plate_2) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_somite_count.html"), selfcontained = FALSE, libdir = "tmp")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.6) +
    theme_void() +
    scale_color_manual(values=LPM_E85_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_sub_clustering.png"), width = 8, height = 6, dpi = 300)

pd$somite_count = factor(pd$somite_count, levels = somite_list[somite_list %in% pd$somite_count])

x_table = table(pd$somite_count)
pd_1 = pd[pd$somite_count %in% names(x_table)[x_table <= 5000],]
pd_2 = pd[pd$somite_count %in% names(x_table)[x_table > 5000],]
pd_2_x = pd_2 %>% group_by(somite_count) %>% sample_n(5000)
pd_sub = rbind(pd_1, pd_2[pd_2$cell_id %in% pd_2_x$cell_id,])

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.4) +
    scale_color_manual(values=somite_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite0-12.day.2D_UMAP.png"), width = 8, height = 6, dpi = 300)


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

gene_count = gene_count[,rownames(pd)]

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_1", "UMAP_2", "UMAP_3")])
saveRDS(cds, paste0(work_path, "/LPM/", example_i, ".cds.rds"))

Idents(obj) = as.vector(obj$leiden_res_1)
res = FindMarkers(obj, ident.1 = "21", only.pos = T)
res_gene = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_ID")

write.csv(res_gene, "~/share/LPM_E85_cluster_21_marker_genes.csv")


######################
### Re Jay's suggestion - Sorry lastly can we try a 2D that excludes derm/scler (and some of the other exclusions)? I'm wondering if we can get a prettier map of the patterned mesoderm if we don't include derivatives (maybe even exclude first/second heart field?)
### Like what if we only included paraxial mesoderm, extraembryonic mesoderm, splanchnic/somatic and cardiopharyngeal. Maybe that would provide a cleaner 2D map that we could point to regions that various derivatives launched from?

### This will be a backbone

source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM_E85"

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

pd_x = readRDS("/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/NovaSeq/run_4_E85/pd_3D.rds")
pd_y = pd %>% left_join(pd_x[,c("cell_id","celltype_sub_clustering")], by = "cell_id")
pd$celltype_sub_clustering = as.vector(pd_y$celltype_sub_clustering)


pd_sub = pd_x[pd_x$celltype_sub_clustering %in% c("Mesodermal progenitors (Tbx6+)",
                                                  "LPM:Extraembryonic mesoderm",
                                                  "LPM:Somatic mesoderm",
                                                  "LPM:Splanchnic mesoderm",
                                                  "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                                                  "Neuromesodermal progenitors"),]
rownames(pd_sub) = as.vector(pd_sub$cell_id)

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = gsub("beth", "embryo", names(table(as.vector(pd_sub$embryo_id))))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE])
}

gene_count = gene_count[,rownames(pd_sub)]

writeMM(t(gene_count), paste0(work_path, "/LPM/", example_i, "_sub.gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/LPM/", example_i, "_sub.df_cell.csv"))



source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM_E85_sub"

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.6) +
    theme_void() +
    scale_color_manual(values=LPM_E85_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_sub_clustering.png"), width = 8, height = 6, dpi = 300)

pd$somite_stage = factor(pd$somite_stage, levels = paste0("stage_", c(0, 2:12)))

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_stage), size=0.6) +
    theme_void() +
    scale_color_viridis(discrete=TRUE) 
theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite_count.png"), width = 8, height = 6, dpi = 300)


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

celltype_include = c(names(table(pd$celltype_update[pd$major_trajectory %in% c("Mesoderm")])), 
                     "Neuromesodermal progenitors", "Anterior intermediate mesoderm", "Posterior intermediate mesoderm",
                     "First heart field", "Second heart field")

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
    geom_bar(stat='identity') + facet_grid(celltype_update ~ ., switch = "y") + 
    theme(strip.text.y.left = element_text(angle = 0)) +
    labs(x='',y='Log2(Estimated # of cells)') +
 #   scale_fill_manual(values=LPM_E85_color_plate) +
 #   theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0("~/share/", "LPM_pre_E10_density.pdf"), 10, 12)
p
dev.off()





################################
### Re Jay's suggestion, only embedding somites 7-16
################################


source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM_somite_7_16"

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd = pd_all[pd_all$somite_count %in% paste0(c(7:16), " somites"),]

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

pd_sub = pd_x[pd_x$celltype_update %in% c("Anterior intermediate mesoderm",
                                          "First heart field",
                                          "LPM:Allantois",
                                          "LPM:Amniotic mesoderm",
                                          "LPM:Cardiopharyngeal mesoderm (Tbx1+)",
                                          "LPM:Extraembryonic mesoderm",
                                          "LPM:Foregut mesenchyme",
                                          "LPM:Gut mesenchyme",
                                          "LPM:Hepatic mesenchyme",
                                          "LPM:Proepicardium (Tbx18+)",
                                          "LPM:Renal capsule",
                                          "LPM:Somatic mesoderm",
                                          "LPM:Splanchnic mesoderm",
                                          "Mesodermal progenitors (Tbx6+)",
                                          "Posterior intermediate mesoderm",
                                          "Second heart field",
                                          "Neuromesodermal progenitors"),]
rownames(pd_sub) = as.vector(pd_sub$cell_id)

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = gsub("beth", "embryo", names(table(as.vector(pd_sub$embryo_id))))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE])
}

gene_count = gene_count[,rownames(pd_sub)]

writeMM(t(gene_count), paste0(work_path, "/LPM/", example_i, ".gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/LPM/", example_i, ".df_cell.csv"))

### After running Scanpy to get embedding

pd = read.csv(paste0(work_path, "/LPM/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.6) +
    theme_void() +
    scale_color_manual(values=LPM_E85_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_update.png"), width = 8, height = 6, dpi = 300)

pd$somite_count = factor(pd$somite_count, levels = paste0(paste0(c(7:16), " somites")))

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.6) +
    theme_void() +
    scale_color_viridis(discrete=TRUE) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite_count.png"), width = 8, height = 6, dpi = 300)

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update, colors = LPM_E85_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

pd$my_cluster = paste0("cluster_", pd$leiden_res_1)
fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~my_cluster) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_leiden_res_1.html"), selfcontained = FALSE, libdir = "tmp")

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) = as.vector(obj$leiden_res_1)
res = FindMarkers(obj, ident.1 = "25",  only.pos = T) 
res = res %>% mutate(gene_ID = rownames(res)) %>% 
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>% filter(p_val_adj < 0.05)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/LPM/", example_i, ".cds.rds"))















