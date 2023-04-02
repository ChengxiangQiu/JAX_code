#########################################
### subclustering on renal trajectory ###
#########################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Renal"; print(example_i)

obj = readRDS(paste0(work_path, "/Renal/", example_i, "_obj.rds"))

Idents(obj) = as.vector(obj$leiden_res_1)
obj_sub = subset(obj, downsample = 1000)
res = FindMarkers(obj_sub, ident.1 = c(5), ident.2 = c(1), only.pos = T)
res = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID")
head(res, 20)

pd = read.csv(paste0(work_path, "/Renal/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)

celltype = rep(NA, nrow(pd))
celltype[pd$leiden_res_1 %in% c(5)] = "Posterior intermediate mesoderm"
celltype[pd$leiden_res_1 %in% c(1)] = "Anterior intermediate mesoderm"
celltype[pd$leiden_res_1 %in% c(7,18)] = "Ureteric bud"
celltype[pd$leiden_res_1 %in% c(21)] = "Collecting duct intercalated cells"
celltype[pd$leiden_res_1 %in% c(17)] = "Collecting duct principal cells"
celltype[pd$leiden_res_1 %in% c(0,3,6,9)] = "Metanephric mesenchyme"
celltype[pd$leiden_res_1 %in% c(2,22,12,10,16,15)] = "Nephron progenitors"
celltype[pd$leiden_res_1 %in% c(14)] = "Podocytes"
celltype[pd$leiden_res_1 %in% c(4,20,8,13)] = "Proximal tubule cells"
celltype[pd$leiden_res_1 %in% c(11)] = "Ascending loop of Henle"
celltype[pd$leiden_res_1 %in% c(19)] = "Distal convoluted tubule"

celltype[pd$leiden_res_5 %in% c(25)] = "Anterior intermediate mesoderm"
celltype[pd$leiden_res_5 %in% c(37,53,44,23)] = "Podocytes"

celltype[pd$leiden_res_20 %in% c(236,51,85)] = "Anterior intermediate mesoderm"
celltype[pd$leiden_res_20 %in% c(38,260)] = "Connecting tubule"
celltype[pd$leiden_res_20 %in% c(117,69,210,219)] = "Distal convoluted tubule"

pd$celltype_sub_clustering = as.vector(celltype)


saveRDS(pd, paste0(work_path, "/Renal/", example_i, "_adata_scale.obs.rds"))

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_sub_clustering)
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_sub_clustering.html"), selfcontained = FALSE, libdir = "tmp")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/mtx/plot/", example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.5) +
    scale_color_manual(values=day_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Renal/", example_i, ".day.2D_UMAP.png"), width = 6, height = 6, dpi = 300)



example_i = "Renal"
pd = readRDS(paste0(work_path, "/Renal/", example_i, "_adata_scale.obs.rds"))

x = as.vector(pd$day)
x[pd$day == "E9"] = "E9.0"
x[pd$day == "E10"] = "E10.0"
x[pd$day == "E11"] = "E11.0"
x[pd$day == "E12"] = "E12.0"
x[pd$day == "E13"] = "E13.0"
x[pd$day == "E14"] = "E14.0"
x[pd$day == "E15"] = "E15.0"
x[pd$day == "E16"] = "E16.0"
x[pd$day == "E17"] = "E17.0"
x[pd$day == "E18"] = "E18.0"
pd$day = as.vector(x)

x = data.frame(day = names(day_color_plate), day_id = c(0:42))
pd = pd %>% left_join(x, by = "day")
pd$day_id = as.numeric((pd$day_id))

pd$celltype_sub_clustering = factor(pd$celltype_sub_clustering, levels = rev(names(renal_color_plate)))

### normalization by the total number of cells from each stage
pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
x = as.vector(pd_all$day)
x[pd_all$day == "E9"] = "E9.0"
x[pd_all$day == "E10"] = "E10.0"
x[pd_all$day == "E11"] = "E11.0"
x[pd_all$day == "E12"] = "E12.0"
x[pd_all$day == "E13"] = "E13.0"
x[pd_all$day == "E14"] = "E14.0"
x[pd_all$day == "E15"] = "E15.0"
x[pd_all$day == "E16"] = "E16.0"
x[pd_all$day == "E17"] = "E17.0"
x[pd_all$day == "E18"] = "E18.0"
pd_all$day = as.vector(x)
x = pd %>% group_by(day, celltype_sub_clustering, day_id) %>% tally() %>% left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day")
x$normalize_num = ceiling(100000 * x$n/x$total_n)
pd_sim = data.frame(celltype_sub_clustering = rep(x$celltype_sub_clustering, x$normalize_num),
                    day_id = rep(x$day_id, x$normalize_num))

library(ggridges)
p=ggplot(pd_sim, aes(x = day_id, y = celltype_sub_clustering, fill = celltype_sub_clustering)) +
    geom_density_ridges() +
    scale_fill_manual(values=renal_color_plate) +
    theme_ridges() + 
    labs(x = "", y = "") +
    theme(legend.position = "none")
pdf(paste0(work_path, "/Renal/Renal_density.pdf"),10,10)
p
dev.off()


cds = readRDS("~/work/jax/rna_seq/Renal/Renal_cds2.rds")

my_plot_cells(cds, genes = c("Pax2", "Pax8", "Sim1", "Lhx1", 
                             "Ret", "Wnt11", "Gdnf", "Wt1", 
                             "Osr1", "Hoxc6","Six2","Eya1"), how_many_rows = 3) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("~/work/jax/rna_seq/Renal/renal_1.png", dpi = 300, height = 6, width = 8)



my_plot_cells(cds, genes = c("Nphs1","Nphs2","Slc27a2","Lrp2",
                             "Umod","Slc12a1","Slc12a3","Pvalb",
                             "Atp6v1g3","Atp6v0d2","Aqp2","Hsd11b2"), how_many_rows = 3) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("~/work/jax/rna_seq/Renal/renal_2.png", dpi = 300, height = 6, width = 8)


my_plot_cells(cds, genes = c("Wnt11","Etv5","Wnt7b","Tacstd2"), how_many_rows = 1) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("~/work/jax/rna_seq/Renal/renal_3.png", dpi = 300, height = 2, width = 8)


my_plot_cells(cds, genes = c("Aqp2","Aqp3","Aqp4","Aqp5"), how_many_rows = 2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("~/work/jax/rna_seq/Renal/renal_4.png", dpi = 300, height = 4, width = 4)




my_plot_cells(cds, genes = c("Wnt11","Etv5","Wnt7b","Tacstd2"), how_many_rows = 2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("~/work/jax/rna_seq/Renal/renal_4.png", dpi = 300, height = 4, width = 4)


p = my_plot_cells(cds[,sample(1:ncol(cds),5000)], genes = c("Wnt11","Ret","Etv5","Wnt7b","Tacstd2"), how_many_rows = 1)
pdf("~/work/jax/rna_seq/Renal/test.pdf")
p
dev.off()


###
### reannotte those cells in the epithelium trajectory (see step4_annotation.R)
###


### output data matrix in 10X format, for Jay and Beth ###

obj = readRDS(paste0(work_path, "/Renal/Renal_obj.rds"))
pd = readRDS(paste0(work_path, "/Renal/Renal_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)
pd = pd[colnames(obj),]

count = GetAssayData(obj, slot = "counts")

gene = mouse_gene[rownames(count),]
gene$mean_exp = Matrix::rowMeans(count)
gene_x = gene %>% group_by(gene_short_name) %>% slice_max(order_by = mean_exp, n = 1) 
gene_1 = gene[rownames(gene) %in% gene_x$gene_ID,]
gene_2 = gene[!rownames(gene) %in% gene_x$gene_ID,]
gene_2$gene_short_name = paste0(gene_2$gene_short_name, ".dup")
gene = rbind(gene_1, gene_2)
gene = gene[rownames(count),]

rownames(count) = as.vector(gene$gene_short_name)

save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data/share"

Matrix::writeMM(count, paste0(save_path, "/matrix.mtx"))
write.table(rownames(count), paste0(save_path, "/genes.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(colnames(count), paste0(save_path, "/barcodes.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(as.vector(pd$celltype_sub_clustering), paste0(save_path, "/celltype.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(as.vector(pd$day), paste0(save_path, "/day.tsv"), row.names=F, col.names=F, sep="\t", quote=F)


pd_sub = pd[,!colnames(pd) %in% c("celltype","keep", "leiden", "celltype_update", "celltype_sub_clustering", "UMAP_1", 
                                  "UMAP_2", "UMAP_3","UMAP_2d_1","UMAP_2d_2", "group", "batch", "cell_id")]
pd_sub$celltye = pd$celltype_sub_clustering
write.table(pd_sub, paste0(save_path, "/meta.tsv"), row.names=T, col.names=T, sep="\t", quote=F)



#### plot day group ####


example_i = "Renal"
pd = readRDS(paste0(work_path, "/Renal/", example_i, "_adata_scale.obs.rds"))

x = as.vector(pd$day)
x[pd$day == "E9"] = "E9.0"
x[pd$day == "E10"] = "E10.0"
x[pd$day == "E11"] = "E11.0"
x[pd$day == "E12"] = "E12.0"
x[pd$day == "E13"] = "E13.0"
x[pd$day == "E14"] = "E14.0"
x[pd$day == "E15"] = "E15.0"
x[pd$day == "E16"] = "E16.0"
x[pd$day == "E17"] = "E17.0"
x[pd$day == "E18"] = "E18.0"
pd$day = as.vector(x)

day_group = rep("Before E10", nrow(pd))

day_group[pd$day %in% c("E10.0", "E10.25", "E10.5", "E10.75")] = "E10-E11"
day_group[pd$day %in% c("E11.0", "E11.25", "E11.5", "E11.75")] = "E11-E12"
day_group[pd$day %in% c("E12.0", "E12.25", "E12.5", "E12.75")] = "E12-E13"
day_group[pd$day %in% c("E13.0", "E13.25", "E13.5", "E13.75")] = "E13-E14"
day_group[pd$day %in% c("E14.0", "E14.25", "E14.333", "E14.75")] = "E14-E15"
day_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"
day_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"
day_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"
day_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"
day_group[pd$day %in% c("P0")] = "P0"

day_group_list = c("Before E10", "E10-E11", "E11-E12", "E12-E13", "E13-E14", "E14-E15", "E15-E16", 
                   "E16-E17", "E17-E18", "E18-P0", "P0")
pd$day_group = factor(day_group, levels = day_group_list)

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_group_list))
names(day_color_plate_2) = day_group_list

day_group_table = table(pd$day_group)
pd_1 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table > 5000]) %>% group_by(day_group) %>% sample_n(5000) %>% as.data.frame()
pd_2 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table <= 5000]) %>% as.data.frame()
pd_sub = rbind(pd_1, pd_2)

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.5) +
    scale_color_manual(values=day_color_plate_2) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".day_group.2D_UMAP.png"), width = 6, height = 6, dpi = 300)






### normalization by the total number of cells from each stage ####
pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))

cell_num = readRDS(paste0(work_path, "/embryo/cell_num_prediction.rds"))

pd$celltype_update = as.vector(pd$celltype_sub_clustering)
pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])

x = pd %>% group_by(day, celltype_update) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    left_join(cell_num %>% select(day, cell_num_pred), by = "day") %>%
    mutate(estimated_num_log2 = log2(ceiling(cell_num_pred*n/total_n))) %>%
    mutate(normalized_num_log2 = log2(ceiling(100000*n/total_n))) 

x$day = factor(x$day, levels = day_list[day_list %in% x$day])
x$celltype_update = factor(x$celltype_update, levels = names(renal_color_plate))

x$day_value = as.vector(x$day)
x$day_value[x$day == "P0"] = "E19"
x$day_value = as.numeric(gsub("E","",x$day_value))

x_first_day = x %>% group_by(celltype_update) %>% filter(n >= 10) %>% 
    slice_min(order_by = day_value, n = 1) %>% mutate(first_day_value = day_value)

x_sub = x %>% left_join(x_first_day %>% select(celltype_update, first_day_value), by = "celltype_update") %>%
    filter(day_value >= first_day_value)


p = x_sub %>% 
    ggplot(aes(x=day, y=estimated_num_log2, fill = celltype_update)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(celltype_update)) + 
    labs(x='',y='Log2(Estimated # of cells)') +
    scale_fill_manual(values=renal_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0("~/share/", "Renal_density_2.pdf"), 6, 8)
p
dev.off()


################
### Re-embedding 4, 5, 12 to find subpopulations

source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol8/projects/cxqiu/work/jax/rna_seq"

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
pd_all$somite_count[is.na(pd_all$somite_count)] = "missing"

pd_sub = pd_all[pd_all$celltype_update %in% c("Collecting duct intercalated cells",
                                              "Connecting tubule",
                                              "Collecting duct principal cells"),]
pd_sub = pd_sub[,c("cell_id", "day", "somite_count", "embryo_id", "embryo_sex",
                   "sequencing_batch", "log2_umi", "G2M.Score", "S.Score",
                   "major_trajectory", "celltype_update")]

celltype_list = read.table(paste0(work_path, "/mtx/celltype_list.txt"), header=F, as.is=T, sep="\t")
names(celltype_list) = c("celltype_update", "celltype_num", "celltype_name")
celltype_list = celltype_list[celltype_list$celltype_update %in% c("Collecting duct intercalated cells",
                                                                   "Connecting tubule",
                                                                   "Collecting duct principal cells"),]

df_gene = read.csv("/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mtx/df_gene.csv", header=T, row.names=1, as.is=T)

gene_count = NULL
cell_name_list = NULL
for(i in 1:nrow(celltype_list)){
    print(celltype_list$celltype_name[i])
    gene_count_i = readMM(paste0(save_path, "/celltype/", celltype_list$celltype_name[i], ".gene_count.mtx"))
    pd_i = read.csv(paste0(save_path, "/celltype/", celltype_list$celltype_name[i], ".df_cell.csv"), row.names=1, header=T)
    gene_count_i = t(gene_count_i)
    rownames(gene_count_i) = rownames(df_gene)
    cell_name_list = c(cell_name_list, rownames(pd_i))
    gene_count = cbind(gene_count, gene_count_i)
}
colnames(gene_count) = as.vector(cell_name_list)
gene_count = gene_count[,rownames(pd_sub)]

example_i = "Renal_CDI"
Matrix::writeMM(t(gene_count), paste0(work_path, "/Renal/", example_i, ".gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, "/Renal/", example_i, ".df_cell.csv"))


pd = read.csv(paste0(work_path, "/Renal/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Renal/", example_i, ".cds.rds"))



day_group = rep("Before E15", nrow(pd))

day_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"
day_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"
day_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"
day_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"
day_group[pd$day %in% c("P0")] = "P0"

day_group_list = c("Before E15", "E15-E16", 
                   "E16-E17", "E17-E18", "E18-P0", "P0")
pd$day_group = factor(day_group, levels = day_group_list)

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_group_list))
names(day_color_plate_2) = day_group_list

day_group_table = table(pd$day_group)
pd_1 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table > 1000]) %>% group_by(day_group) %>% sample_n(500) %>% as.data.frame()
pd_2 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table <= 1000]) %>% as.data.frame()
pd_sub = rbind(pd_1, pd_2)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=1) +
    scale_color_manual(values=day_color_plate_2) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".day_group.2D_UMAP.png"), width = 5, height = 5, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=1) +
    scale_color_manual(values=renal_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_update.2D_UMAP.png"), width = 5, height = 5, dpi = 300)


p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.3) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Renal/", example_i, ".2D_UMAP.day.png"), width = 6, height = 6, dpi = 300)

saveRDS(pd, paste0(work_path, "/Renal/", example_i, "_adata_scale.obs.rds"))


gene_list = c("Atp6v1b1","Kit","Slc4a1","Slc26a4","Aqp2","Aqp4")

my_plot_cells(cds, genes = gene_list, how_many_rows = 2, cell_size = 0.6) 

my_plot_cells(cds, genes = gene_list, how_many_rows = 3, cell_size = 0.6) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("renal_5.png", dpi = 300, height = 6, width = 4)




