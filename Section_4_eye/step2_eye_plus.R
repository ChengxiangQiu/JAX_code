day_include = c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
"E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
"E12.25", "E12.5", "E12.75", "E13.0")

celltype_include = c("Amacrine cells",
                     "Amacrine/Horizontal precursor cells",
                     "Bipolar precursor cells",
                     "Cholinergic amacrine cells",
                     "Ciliary margin cells",
                     "Cone precursor cells",
                     "Horizontal cells",
                     "Naive retinal progenitor cells",
                     "Photoreceptor precursor cells",
                     "PV-containing retinal ganglion cells",
                     "Retinal ganglion cells",
                     "Retinal progenitor cells",
                     "Rod precursor cells",
                     "Eye field",
                     "Retinal pigment cells")

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Eye_plus2"
pd = read.csv(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])

pd$leiden_res_5 == 10 ### Optic stalk cells
pd$leiden_res_20 == 228 ### Iris pigment epithelium


day_group = rep("Before E13", nrow(pd))

day_group[pd$day %in% c("E13.0", "E13.25", "E13.5", "E13.75")] = "E13-E14"
day_group[pd$day %in% c("E14.0", "E14.25", "E14.333", "E14.75")] = "E14-E15"
day_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"
day_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"
day_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"
day_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"
day_group[pd$day %in% c("P0")] = "P0"

pd$day_group = factor(day_group, levels = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                            "E16-E17", "E17-E18", "E18-P0", "P0"))

library(RColorBrewer)
day_group_color_plate=rev(brewer.pal(11,"Spectral"))
day_group_color_plate=colorRampPalette(day_group_color_plate)(8)
names(day_group_color_plate) = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                 "E16-E17", "E17-E18", "E18-P0", "P0")

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update, colors = eye_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day, colors = day_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_day.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd %>% filter(celltype_update == "Eye field"), x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day, colors = day_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_Eye_field_day.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day_group, colors = day_group_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_day_group.html"), selfcontained = FALSE, libdir = "tmp")

pd$my_cluster = paste0("cluster_", pd$leiden_res_5)
fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~my_cluster) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_cluster.html"), selfcontained = FALSE, libdir = "tmp")


x_table = table(pd$day_group)
pd_1 = pd[pd$day_group %in% names(x_table)[x_table <= 5000],]
pd_2 = pd[pd$day_group %in% names(x_table)[x_table > 5000],]
pd_2_x = pd_2 %>% group_by(day_group) %>% sample_n(5000)
pd_sub = rbind(pd_1, pd_2[pd_2$cell_id %in% pd_2_x$cell_id,])

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day_group, colors = day_group_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_day_group.html"), selfcontained = FALSE, libdir = "tmp")


### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
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
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(obj, paste0(work_path, "/Eye/", example_i, ".obj.rds"))

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Eye/", example_i, ".cds.rds"))


p = pd_sub[sample(1:nrow(pd_sub)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.3) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.day_group_downsample.png"), width = 6, height = 6, dpi = 300)

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.celltype_update_downsample.png"), width = 6, height = 6, dpi = 300)





### normalization by the total number of cells from each stage ####
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


cell_num = readRDS(paste0(work_path, "/embryo/cell_num_prediction.rds"))

x = pd %>% group_by(day, celltype_update) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    left_join(cell_num %>% select(day, cell_num_pred), by = "day") %>%
    mutate(estimated_num_log2 = log2(ceiling(cell_num_pred*n/total_n))) %>%
    mutate(normalized_num_log2 = log2(ceiling(100000*n/total_n))) 

x$day = factor(x$day, levels = day_list[day_list %in% x$day])
x$celltype_update = factor(x$celltype_update, levels = names(eye_color_plate))

x$day_value = as.vector(x$day)
x$day_value[x$day == "P0"] = "E19"
x$day_value = as.numeric(gsub("E","",x$day_value))

x_first_day = x %>% group_by(celltype_update) %>% filter(n >= 10) %>% 
    slice_min(order_by = day_value, n = 1) %>% mutate(first_day_value = day_value)

x_sub = x %>% left_join(x_first_day %>% select(celltype_update, first_day_value), by = "celltype_update") %>%
    filter(day_value >= first_day_value)

p = x %>%
    ggplot(aes(day_value, normalized_num_log2, color = celltype_update)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 1), axis.text.y = element_text(color="black")) +
    scale_x_continuous(breaks=seq(8.5, 19, 0.5)) +
    scale_color_manual(values=eye_color_plate) +
    NoLegend()

pdf(paste0(work_path, "/Eye/Eye_density_2.pdf"), 8, 5)
p
dev.off()



p = x_sub %>% 
    ggplot(aes(x=day, y=estimated_num_log2, fill = celltype_update)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(celltype_update)) + 
    labs(x='',y='Log2(Estimated # of cells)') +
    scale_fill_manual(values=eye_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0(work_path, "/Eye/Eye_density_3.pdf"), 6.5, 12)
p
dev.off()



cds_sub = cds[,sample(1:ncol(cds), 50000)]

my_plot_cells(cds_sub, genes = c("Pax2","Rax","Tyr","Oca2","Ccnd1","Sox2"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_1.png", dpi = 300, height = 2, width = 12)

my_plot_cells(cds_sub, genes = c("Neurog2","Otx2","Vsx2","Neurod1","Otx1","Gja1"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_2.png", dpi = 300, height = 2, width = 12)

my_plot_cells(cds_sub, genes = c("Neurod4","Crx","Cnga3","Thrb","Nrl","Nr2e3"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_3.png", dpi = 300, height = 2, width = 12)

my_plot_cells(cds_sub, genes = c("Pou4f2","Isl1","Pou6f2","Pvalb","Ptf1a","Tfap2b"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_4.png", dpi = 300, height = 2, width = 12)

my_plot_cells(cds_sub, genes = c("Gad1","Nr4a2","Sall3","Lhx1","Ebf1","Chat"), how_many_rows = 1, cell_size = 0.2) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_5.png", dpi = 300, height = 2, width = 12)




#########
### Focusing on days <= E12.5


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Eye_plus2_early"
pd = read.csv(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])

obj = readRDS(paste0(work_path, "/Eye/Eye_plus2.obj.rds"))
gene_count = GetAssayData(obj, slot = "counts")
gene_count = gene_count[,rownames(pd)]
obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Eye/", example_i, ".cds.rds"))

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.35) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.celltype_update.png"), width = 6, height = 6, dpi = 300)



pd$day_x = as.numeric(as.vector(gsub("E","", as.vector(pd$day))))

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_x), size=0.35) +
    theme_void() +
    scale_color_viridis(option = "B", breaks = c(8.5, 9.5, 10.5, 11.5, 12.5), labels = c("E8.5","E9.5","E10.5","E11.5","E12.5")) + 
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.day.png"), width = 6, height = 6, dpi = 300)

pdf(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.day.pdf"))
pd[sample(1:nrow(pd), 1000),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_x), size=0.35) +
    theme_void() +
    scale_color_viridis(option = "B", breaks = c(8.5, 9.5, 10.5, 11.5, 12.5), labels = c("E8.5","E9.5","E10.5","E11.5","E12.5")) 
dev.off()

my_plot_cells(cds_1, genes = c("Pax2","Pax6","Rax",
                               "Vax1","Tyr","Fgf15"), how_many_rows = 2, cell_size = 0.35) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_6.png", dpi = 300, height = 4, width = 6)









#########
### Focusing on iris


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Eye_plus2_iris"
pd = read.csv(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)
rownames(pd) = as.vector(pd$cell_id)

pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])

day_group = rep("Before E13", nrow(pd))

day_group[pd$day %in% c("E13.0", "E13.25", "E13.5", "E13.75")] = "E13-E14"
day_group[pd$day %in% c("E14.0", "E14.25", "E14.333", "E14.75")] = "E14-E15"
day_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"
day_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"
day_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"
day_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"
day_group[pd$day %in% c("P0")] = "P0"

pd$day_group = factor(day_group, levels = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                            "E16-E17", "E17-E18", "E18-P0", "P0"))


obj = readRDS(paste0(work_path, "/Eye/Eye_plus2.obj.rds"))
gene_count = GetAssayData(obj, slot = "counts")
gene_count = gene_count[,rownames(pd)]
obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Eye/", example_i, ".cds.rds"))



p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.3) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.day_group.png"), width = 6, height = 6, dpi = 300)

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.celltype_update.png"), width = 6, height = 6, dpi = 300)


my_plot_cells(cds, genes = c("Cdh18", "Zic1", "Tbx20", "Wnt16","Tyr","Oca2"), how_many_rows = 2, cell_size = 0.3) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_7.png", dpi = 300, height = 4, width = 6)





