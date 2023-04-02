#########################################
### subclustering on Eye trajectory ###
#########################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Eye"; print(example_i)

obj = readRDS(paste0(work_path, "/Eye/", example_i, "_obj.rds"))

Idents(obj) = as.vector(obj$leiden_res_1)
obj_sub = subset(obj, downsample = 1000)
res = FindMarkers(obj_sub, ident.1 = c(18), only.pos = T)
res = FindMarkers(obj_sub, ident.1 = c(18), ident.2 = c(11,12,20), only.pos = T)
res = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID")
head(res, 20)

pd = read.csv(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)

celltype = rep(NA, nrow(pd))
celltype[pd$leiden_res_1 %in% c(0,1,2,3,4,6)] = "Naive retinal progenitor cells"
celltype[pd$leiden_res_1 %in% c(5,15)] = "Retinal progenitor cells"
celltype[pd$leiden_res_1 %in% c(14)] = "Bipolar precursor cells"
celltype[pd$leiden_res_1 %in% c(7)] = "Amacrine/Horizontal precursor cells"
celltype[pd$leiden_res_1 %in% c(11,12,18)] = "Amacrine cells"
celltype[pd$leiden_res_1 %in% c(20)] = "Cholinergic amacrine cells"
celltype[pd$leiden_res_1 %in% c(22)] = "Horizontal cells"
celltype[pd$leiden_res_1 %in% c(8)] = "Photoreceptor precursor cells"
celltype[pd$leiden_res_1 %in% c(19,10)] = "Cone precursor cells"
celltype[pd$leiden_res_1 %in% c(17)] = "Rod precursor cells"
celltype[pd$leiden_res_1 %in% c(13,9,16)] = "Retinal ganglion cells"
celltype[pd$leiden_res_1 %in% c(21)] = "PV-containing retinal ganglion cells"

celltype[pd$celltype_update == "Ciliary margin cells"] = "Ciliary margin cells"

pd_1 = readRDS(paste0(work_path, "/mtx/sub_clustering/CNS_neurons_pd.rds"))
pd_2 = readRDS(paste0(work_path, "/mtx/sub_clustering/Eye_and_other_pd.rds"))
pd_x = rbind(pd_1[,c("leiden_res_1","celltype_update")], pd_2[,c("leiden_res_1","celltype_update")])
pd_x = pd_x[rownames(pd),]
pd$orig_cluster = as.vector(pd_x$leiden_res_1)

x = celltype

x[pd$major_trajectory == "CNS_neurons" & 
      celltype %in% c("Naive retinal progenitor cells","Retinal progenitor cells","Bipolar precursor cells","Ciliary margin cells","Photoreceptor precursor cells","Cone precursor cells","Rod precursor cells") &
      pd$orig %in% c(35)] = "Amacrine/Horizontal precursor cells"
x[pd$major_trajectory == "CNS_neurons" & 
      celltype %in% c("Naive retinal progenitor cells","Retinal progenitor cells","Bipolar precursor cells","Ciliary margin cells","Photoreceptor precursor cells","Cone precursor cells","Rod precursor cells") &
      pd$orig %in% c(37)] = "Retinal ganglion cells"

x[pd$major_trajectory == "Eye_and_other" & 
      celltype %in% c("Retinal ganglion cells","PV-containing retinal ganglion cells","Amacrine/Horizontal precursor cells","Amacrine cells","Horizontal cells","Cholinergic amacrine cells") &
      pd$orig %in% c(2)] = "Naive retinal progenitor cells"
x[pd$major_trajectory == "Eye_and_other" & 
      celltype %in% c("Retinal ganglion cells","PV-containing retinal ganglion cells","Amacrine/Horizontal precursor cells","Amacrine cells","Horizontal cells","Cholinergic amacrine cells") &
      pd$orig %in% c(4,7,14)] = "Retinal progenitor cells"
x[pd$major_trajectory == "Eye_and_other" & 
      celltype %in% c("Retinal ganglion cells","PV-containing retinal ganglion cells","Amacrine/Horizontal precursor cells","Amacrine cells","Horizontal cells","Cholinergic amacrine cells") &
      pd$orig %in% c(5)] = "Cone precursor cells"
x[pd$major_trajectory == "Eye_and_other" & 
      celltype %in% c("Retinal ganglion cells","PV-containing retinal ganglion cells","Amacrine/Horizontal precursor cells","Amacrine cells","Horizontal cells","Cholinergic amacrine cells") &
      pd$orig %in% c(9)] = "Photoreceptor precursor cells"

celltype = x

pd$celltype_sub_clustering = as.vector(celltype)


saveRDS(pd, paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.rds"))


fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_sub_clustering)
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_sub_clustering.html"), selfcontained = FALSE, libdir = "tmp")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.3) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)



example_i = "Eye"
pd = readRDS(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.rds"))

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

pd$celltype_sub_clustering = factor(pd$celltype_sub_clustering, levels = rev(names(eye_color_plate)))

###### remaking the 2D UMAP colored by day

color_group = rep(NA, nrow(pd))

color_group[pd$day %in% c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
                  "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", "E12.25", "E12.5", "E12.75")] = "Before E13"

color_group[pd$day %in% c("E13.0", "E13.25", "E13.5", "E13.75")] = "E13-E14"

color_group[pd$day %in% c("E14.0", "E14.25", "E14.333", "E14.75")] = "E14-E15"

color_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"

color_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"

color_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"

color_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"

color_group[pd$day %in% c("P0")] = "P0"

pd$color_group = factor(color_group, levels = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                                "E16-E17", "E17-E18", "E18-P0", "P0"))
x_table = table(pd$color_group)
pd_1 = pd[pd$color_group %in% names(x_table)[x_table <= 5000],]
pd_2 = pd[pd$color_group %in% names(x_table)[x_table > 5000],]
pd_2_x = pd_2 %>% group_by(color_group) %>% sample_n(5000)
pd_sub = rbind(pd_1, pd_2[pd_2$cell_id %in% pd_2_x$cell_id,])

color_group_plate = c("Before E13" = "#5E4FA2", 
                      "E13-E14" = "#48A0B2", 
                      "E14-E15" = "#A1D9A4", 
                      "E15-E16" = "#EDF7A3", 
                      "E16-E17" = "#FEE899", 
                      "E17-E18" = "#FBA45C", 
                      "E18-P0" = "#E25249", 
                      "P0" = "#9E0142")

p = pd_sub[sample(1:nrow(pd_sub)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = color_group), size=0.5) +
    theme_void() +
    scale_color_manual(values=color_group_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Eye/", example_i, ".2D_UMAP.day_downsample.png"), width = 6, height = 6, dpi = 300)


df = data.frame(color_group = names(color_group_plate),
                UMAP_2d_1 = sample(1:20, length(color_group_plate)),
                UMAP_2d_2 = sample(1:20, length(color_group_plate)))
pdf("~/work/jax/rna_seq/Eye/Eye.2D_UMAP.day_downsample.pdf", 5, 5)
df %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = color_group), size=0.5) +
    theme_void() +
    scale_color_manual(values=color_group_plate) +
    guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()


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

x = pd %>% group_by(day, celltype_sub_clustering, day_id) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day")
x$normalize_num = ceiling(100000 * x$n/x$total_n)
pd_sim = data.frame(celltype_sub_clustering = rep(x$celltype_sub_clustering, x$normalize_num),
                    day_id = rep(x$day_id, x$normalize_num))

pd_sim$celltype_sub_clustering = factor(pd_sim$celltype_sub_clustering, levels = rev(names(eye_color_plate)))

library(ggridges)
p=ggplot(pd_sim, aes(x = day_id, y = celltype_sub_clustering, fill = celltype_sub_clustering)) +
    geom_density_ridges() +
    scale_fill_manual(values=eye_color_plate) +
    theme_ridges() + 
    labs(x = "", y = "") +
    theme(legend.position = "none")
pdf(paste0(work_path, "/Eye/Eye_density.pdf"),10,10)
p
dev.off()


### Re jay's suggestion
### Can you regenerate the figure that was originally there w/ the timing of various 
### annotated subsets as a set of line graphs? However, can you make it ONE figure with 
### a different color for each line, corresponding to cell type, and the y-axis should 
### be on a log scale. Does that make sense? We want to do absolute numbers of cells 
### (not proportions) but the estimated absolute number should be scaled in such a way 
### that the estimate is for the whole embryo. Does that make sense? So we actually want 
### to estimate literally how many rod precursor cells there are in the entire embryo 
### at each timepoint and plot that.


cell_num = readRDS(paste0(work_path, "/embryo/cell_num_prediction.rds"))

x = pd %>% group_by(day, celltype_sub_clustering) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    left_join(cell_num %>% select(day, cell_num_pred), by = "day") %>%
    mutate(estimate_num_log2 = log2(ceiling(cell_num_pred*n/total_n)))

x$day = factor(x$day, levels = day_list[day_list %in% x$day])
x$celltype_sub_clustering = factor(x$celltype_sub_clustering, levels = names(eye_color_plate))

p = x %>% 
    ggplot(aes(x=day, y=estimate_num_log2, fill = celltype_sub_clustering)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(celltype_sub_clustering)) + 
    labs(x='',y='Log2(Estimated # of cells)') +
    scale_fill_manual(values=eye_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

    
pdf(paste0(work_path, "/Eye/Eye_density_2.pdf"),8,12)
p
dev.off()




###
### reannotte those cells in the CNS_neurons and Eye_and_other trajectory (see step4_annotation.R)
###


my_plot_cells(cds, genes = c("Gad1","Bhlhe22","Nr4a2","Isl1","Chat","Tnc","Runx2","Lhx9"), how_many_rows = 2, cell_size = 0.4) +
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("~/work/jax/figures/tmp/eye_1.png", dpi = 300, height = 3, width = 6)

p = my_plot_cells(cds[,sample(1:ncol(cds),5000)], genes = c("Gad1","Bhlhe22","Nr4a2","Isl1","Chat","Tnc","Runx2","Lhx9"), how_many_rows = 2)
pdf("~/work/jax/figures/tmp/eye_1.pdf")
p
dev.off()




#####################
### Re Jay's suggestion, checking the TF corresponding to each subpopulations of RGCs
#####################


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Eye"
pd = readRDS(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)
obj = readRDS(paste0(work_path, "/Eye/", example_i, "_obj.rds"))
print(table(pd$celltype_sub_clustering))
gene_count = GetAssayData(obj, slot = "counts")

pd_sub = pd[pd$celltype_sub_cluster %in% c("Retinal ganglion cells"),]
gene_count_sub = gene_count[,rownames(pd_sub)]

#writeMM(t(gene_count_sub), paste0(work_path, "/Eye/", example_i, "_RGC_gene_count.mtx"))
#write.csv(pd_sub, paste0(work_path, "/Eye/", example_i, "_RGC_df_cell.csv"))


obj_sub = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
obj_sub$group = "RGC"
obj_sub = doClusterSeurat(obj_sub, doClustering = F)

pd_x = read.csv(paste0(work_path, "/Eye/", example_i, "_RGC_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
print(sum(rownames(pd_x) == colnames(obj_sub)))
emb = pd_x[,c("UMAP_2d_1", "UMAP_2d_2")]
colnames(emb) = c("UMAP_1", "UMAP_2")
rownames(emb) = rownames(emb)
obj_sub@reductions$umap@cell.embeddings = as.matrix(emb)

leiden_cluster = rep(0, nrow(pd_x))
leiden_cluster[pd_x$leiden == 13] = 1
leiden_cluster[pd_x$leiden == 5] = 2
leiden_cluster[pd_x$leiden == 0] = 3
leiden_cluster[pd_x$leiden == 8] = 4
leiden_cluster[pd_x$leiden == 21] = 5
leiden_cluster[pd_x$leiden == 15] = 6
leiden_cluster[pd_x$leiden == 19] = 7
leiden_cluster[pd_x$leiden == 11] = 8
leiden_cluster[pd_x$leiden == 17] = 9
leiden_cluster[pd_x$leiden == 20] = 10
leiden_cluster[pd_x$leiden == 16] = 11
leiden_cluster[pd_x$leiden == 12] = 12
leiden_cluster[pd_x$leiden == 18] = 13
leiden_cluster[pd_x$leiden == 6] = 14
leiden_cluster[pd_x$leiden == 4] = 15
pd_x$leiden_cluster = as.vector(leiden_cluster)

obj_sub$leiden_cluster = as.vector(pd_x$leiden_cluster)
pdf("~/share/RGC.pdf", 5, 5)
DimPlot(obj_sub, group.by = "leiden_cluster", label = T) + NoLegend()
dev.off()


p = pd_x %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = color_group), size=0.45) +
    theme_void() +
    scale_color_manual(values=color_group_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/RGC_day.png"), width = 5, height = 5, dpi = 300)

pd_x$leiden_cluster = paste0("c", as.vector(pd_x$leiden_cluster))
p = pd_x %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = leiden_cluster), size=0.45) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/RGC_leiden_cluster.png"), width = 5, height = 5, dpi = 300)


Idents(obj_sub) = as.vector(obj_sub$leiden_cluster)

mouse_TF = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/TF_list/AnimalTFDB/Mus_musculus_TF.txt", sep="\t", header=T, as.is=T)
mouse_TF_overlap = intersect(as.vector(mouse_TF$Ensembl), rownames(obj_sub))

result = FindAllMarkers(subset(obj_sub, leiden_cluster != 0), features = mouse_TF_overlap)
rownames(result) = NULL
result$gene_ID = result$gene
result = result[result$p_val_adj < 0.05,]
x = result %>% left_join(mouse_gene, by = "gene_ID")
result$gene_short_name = as.vector(x$gene_short_name)
result$gene = NULL
write.csv(result, paste0(work_path, "/Eye/Eye_RGC_top_TFs.csv"), row.names=F)

result_sub = result %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 3)

data = GetAssayData(obj_sub, slot = "data")

result_sub$cluster = factor(result_sub$cluster, levels = c(1:15))
result_sub = result_sub %>% arrange(cluster)

exp = NULL
for(i in c(1:15)){
    print(i)
    exp = cbind(exp, Matrix::rowMeans(data[,obj_sub$leiden_cluster == i]))
}

colnames(exp) = paste0("cluster_", c(1:15))

gene_list = unique(as.vector(result_sub$gene_short_name))
mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% gene_list,]
rownames(mouse_gene_sub) = mouse_gene_sub$gene_short_name
mouse_gene_sub = mouse_gene_sub[gene_list,]
exp_sub = exp[as.vector(mouse_gene_sub$gene_ID),]
rownames(exp_sub) = rownames(mouse_gene_sub)

saveRDS(exp_sub, paste0(work_path, "/Eye/Eye_RGC_heatmap_dat.rds"))

dat = readRDS("~/work/jax/rna_seq/Eye/Eye_RGC_heatmap_dat.rds")

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("~/work/jax/rna_seq/Eye/Eye_RGC_heatmap.pdf", 8, 5)
heatmap.2(as.matrix(t(dat)), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()






##############################################
### doing PCA on Rod vs. Cone populations  ###
##############################################


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Eye"
pd = readRDS(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)
obj = readRDS(paste0(work_path, "/Eye/", example_i, "_obj.rds"))
print(table(pd$celltype_sub_clustering))
gene_count = GetAssayData(obj, slot = "counts")


### extracing the sub-populations 
i = "Photoreceptor"

pd_x = pd[pd$celltype_sub_clustering %in% c("Photoreceptor precursor cells", "Cone precursor cells", "Rod precursor cells"),]
gene_count_x = gene_count[,rownames(pd_x)]
obj_x = CreateSeuratObject(gene_count_x, meta.data = pd_x)

npcs = 30
reduction.key = "PC_"
seed.use = 42

obj_x = NormalizeData(obj_x, normalization.method = "LogNormalize", scale.factor = 10000)
obj_x = FindVariableFeatures(obj_x, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj_x)
obj_x = ScaleData(obj_x, verbose = FALSE, features = rownames(obj_x))
scale_dat = GetAssayData(obj_x, slot = "scale.data")
print(dim(scale_dat))

set.seed(seed = seed.use)
pca.results <- irlba::irlba(A = t(x = scale_dat[genes_include,]), nv = npcs)
feature.loadings <- pca.results$v
set.seed(seed = seed.use)
cell.embeddings <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- genes_include
colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
rownames(x = cell.embeddings) <- colnames(obj_x)
colnames(x = cell.embeddings) <- colnames(x = feature.loadings)

stdev <- pca.results$d/sqrt(max(1, ncol(scale_dat) - 1))
eigValues = (stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)

res = list(cell.embeddings = cell.embeddings,
           feature.loadings = feature.loadings,
           varExplained = varExplained)
saveRDS(res, paste0(work_path, "/Eye/", example_i, "_adata_scale.", i, ".PCA.rds"))

res = readRDS(paste0(work_path, "/Eye/", example_i, "_adata_scale.", i, ".PCA.rds"))
emb = res[["cell.embeddings"]]
emb = emb[rownames(pd_x),]
pd_x = cbind(pd_x, emb[,c(1:3)])
print(res[["varExplained"]])

fig = plot_ly(pd_x, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~celltype_sub_clustering) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (26.1%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (16.8%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (9.9%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_Photoreceptor_PCA_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

pd_x$day = factor(pd_x$day, levels = day_list[day_list %in% pd_x$day])
fig = plot_ly(pd_x, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~day, colors = day_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (26.1%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (16.8%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (9.9%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_Photoreceptor_PCA_day.html"), selfcontained = FALSE, libdir = "tmp")

fig = plot_ly(pd_x, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~G2M.Score) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (26.1%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (16.8%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (9.9%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_Photoreceptor_PCA_G2M.Score.html"), selfcontained = FALSE, libdir = "tmp")

### boxplot of PC1 across timepoints for different subpopulations

pd_x$day = factor(pd_x$day, day_list[day_list %in% pd_x$day])
pd_x$celltype_sub_clustering = factor(pd_x$celltype_sub_clustering, levels = c("Photoreceptor precursor cells",
                                                                               "Cone precursor cells", "Rod precursor cells"))


pdf("Rod_Cone_PC1.pdf", 6, 4.5)
pd_x %>% 
    ggplot(aes(day, PC_1, fill = celltype_sub_clustering)) + geom_boxplot(outlier.shape = NA) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 10) +
    scale_fill_brewer(palette="Set1") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(-20, 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 1), axis.text.y = element_text(color="black"))
dev.off()


#pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
#x = as.vector(pd_all$day)
#x[pd_all$day == "E9"] = "E9.0"
#x[pd_all$day == "E10"] = "E10.0"
#x[pd_all$day == "E11"] = "E11.0"
#x[pd_all$day == "E12"] = "E12.0"
#x[pd_all$day == "E13"] = "E13.0"
#x[pd_all$day == "E14"] = "E14.0"
#x[pd_all$day == "E15"] = "E15.0"
#x[pd_all$day == "E16"] = "E16.0"
#x[pd_all$day == "E17"] = "E17.0"
#x[pd_all$day == "E18"] = "E18.0"
#pd_all$day = as.vector(x)


### making the 2D density plot between real day (x-axis) and PC3 (y-axis)

pseudotime_bin_break = c(min(pd_x$PC_1), seq(-20, max(pd_x$PC_1), length.out=10))
pd_x$pseudotime_bin = cut(pd_x$PC_1,
                         breaks = pseudotime_bin_break,
                         include.lowest = T,
                         right = F) 

pd_x$day = factor(pd_x$day, levels = day_list[day_list %in% pd_x$day])
pd_x$pseudotime_bin = factor(pd_x$pseudotime_bin, levels = levels(pd_x$pseudotime_bin))

day_bin_list = levels(pd_x$day)
pseudotime_bin_list = levels(pd_x$pseudotime_bin)

dat_x = NULL
for(i in day_bin_list){
    for(j in pseudotime_bin_list){
        dat_x = rbind(dat_x, data.frame(day = i, pseudotime_bin = j, n = 0, stringsAsFactors = F))
    }
}


dat = pd_x[pd_x$celltype_sub_clustering %in% c("Cone precursor cells", "Photoreceptor precursor cells"),]

x = dat %>% group_by(day, pseudotime_bin) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day")
x$normalize_num = ceiling(100000 * x$n/x$total_n)

x1 = dat_x %>% left_join(x %>% select(day, pseudotime_bin, normalize_num), by = c("day", "pseudotime_bin"))
x1$normalize_num[is.na(x1$normalize_num)] = 0
x1$pseudotime_bin = factor(x1$pseudotime_bin, levels = rev(pseudotime_bin_list))

library(reshape2)
dat = x1 %>% select(day, pseudotime_bin, normalize_num) %>%
    dcast(pseudotime_bin~day)
dat[is.na(dat)] = 0
rownames(dat) = dat[,1]
dat = dat[,-1]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("~/share/Cone_PC1.pdf", 8, 5)
heatmap.2(as.matrix(dat), 
          col=Colors, 
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()






















###########################################################
### projecting Retinal cells onto the PCA space of NMP  ###
###########################################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "posterior_embryo"

gene_count = readRDS(paste0(work_path, "/somitogenesis/", example_i, "_adata_scale.gene_count.rds"))

i = "NMP_Mesoderm"

pd_x = read.csv(paste0(work_path, "/somitogenesis/", example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)

gene_count_x = gene_count[,rownames(pd_x)]
obj_x = CreateSeuratObject(gene_count_x, meta.data = pd_x)

npcs = 30
reduction.key = "PC_"
seed.use = 42

obj_x = NormalizeData(obj_x, normalization.method = "LogNormalize", scale.factor = 10000)
obj_x = FindVariableFeatures(obj_x, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj_x)
obj_x = ScaleData(obj_x, verbose = FALSE, features = rownames(obj_x))
scale_dat = GetAssayData(obj_x, slot = "scale.data")
print(dim(scale_dat))

set.seed(seed = seed.use)
pca.results <- irlba::irlba(A = t(x = scale_dat[genes_include,]), nv = npcs)
feature.loadings <- pca.results$v
set.seed(seed = seed.use)
cell.embeddings <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- genes_include
colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
rownames(x = cell.embeddings) <- colnames(obj_x)
colnames(x = cell.embeddings) <- colnames(x = feature.loadings)

stdev <- pca.results$d/sqrt(max(1, ncol(scale_dat) - 1))
eigValues = (stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)

pd_x$PC_1 = as.vector(cell.embeddings[,1])
pd_x$PC_2 = as.vector(cell.embeddings[,2])
pd_x$PC_3 = as.vector(cell.embeddings[,3])

pd_1 = pd_x

### calculate the PCA manually
set.seed(seed = seed.use)
irlba_res <- my_sparse_prcomp_irlba(Matrix::t(scale_dat[genes_include,]), 
                                    n = npcs, 
                                    center = F, 
                                    scale. = F)
preproc_res <- irlba_res$x

norm_dat = GetAssayData(obj_x, slot = "data")
norm_dat_sub = norm_dat[genes_include,]
scale_dat_sub = scale_dat[genes_include,]

norm_dat_sub_mean = apply(norm_dat_sub, 1, mean)
norm_dat_sub_sd = apply(norm_dat_sub, 1, sd)

scale_dat_sub_reconstruct = t(scale(t(norm_dat_sub), norm_dat_sub_mean, norm_dat_sub_sd))
cell.embeddings_reconstruct = t(scale_dat_sub_reconstruct) %*% irlba_res$rotation


### now reading EYE data, and projecting it onto the above PCA space
example_i = "Eye"
pd = readRDS(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)
obj = readRDS(paste0(work_path, "/Eye/", example_i, "_obj.rds"))
print(table(pd$celltype_sub_clustering))
gene_count = GetAssayData(obj, slot = "counts")

### extracing the sub-populations 
i = "Photoreceptor"

pd_x = pd[pd$celltype_sub_clustering %in% c("Photoreceptor precursor cells", "Cone precursor cells", "Rod precursor cells"),]
gene_count_x = gene_count[,rownames(pd_x)]
obj_x = CreateSeuratObject(gene_count_x, meta.data = pd_x)
obj_x = NormalizeData(obj_x, normalization.method = "LogNormalize", scale.factor = 10000)

norm_dat = GetAssayData(obj_x, slot = "data")
norm_dat_sub = norm_dat[genes_include,]

scale_dat_sub_reconstruct = t(scale(t(norm_dat_sub), norm_dat_sub_mean, norm_dat_sub_sd))
cell.embeddings_reconstruct = t(scale_dat_sub_reconstruct) %*% irlba_res$rotation

pd_x$PC_1 = as.vector(cell.embeddings_reconstruct[,1])
pd_x$PC_2 = as.vector(cell.embeddings_reconstruct[,2])
pd_x$PC_3 = as.vector(cell.embeddings_reconstruct[,3])

pd_2 = pd_x
pd_2$celltype_update = pd_2$celltype_sub_clustering




df = rbind(pd_1[,c("cell_id", "somite_count", "day", "PC_1", "PC_2", "PC_3", "celltype_update")],
           pd_2[,c("cell_id", "somite_count", "day", "PC_1", "PC_2", "PC_3", "celltype_update")])

color_group = as.vector(df$day)

color_group[df$day %in% c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
                          "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", "E12.25", "E12.5", "E12.75")] = "Before E13"

color_group[df$day %in% c("E13.0", "E13.25", "E13.5", "E13.75")] = "E13-E14"

color_group[df$day %in% c("E14.0", "E14.25", "E14.333", "E14.75")] = "E14-E15"

color_group[df$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"

color_group[df$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"

color_group[df$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"

color_group[df$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"

color_group[df$day %in% c("P0")] = "P0"

color_group[!df$celltype_update %in% c("Photoreceptor precursor cells", "Cone precursor cells", "Rod precursor cells")] = "NMP"

df$color_group = factor(color_group, levels = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                                "E16-E17", "E17-E18", "E18-P0", "P0", "NMP"))

color_group_plate = c("Before E13" = "#5E4FA2", 
                      "E13-E14" = "#48A0B2", 
                      "E14-E15" = "#A1D9A4", 
                      "E15-E16" = "#EDF7A3", 
                      "E16-E17" = "#FEE899", 
                      "E17-E18" = "#FBA45C", 
                      "E18-P0" = "#E25249", 
                      "P0" = "#9E0142",
                      "NMP" = "grey80")


fig = plot_ly(df, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~color_group, colors = color_group_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (21.0%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (13.9%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (11.1%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(save_path, "/example/", example_i, "_Photoreceptor_PCA_test.html"), selfcontained = FALSE, libdir = "tmp")

saveRDS(df, paste0(work_path, "/Eye/", example_i, "_project_PCA_NMPs.rds"))


### making the 2D density plot between real day (x-axis) and PC3 (y-axis)

dat = df[df$celltype_update %in% c("Rod precursor cells", "Photoreceptor precursor cells"),]

pseudotime_bin_num = 50

pseudotime_bin_break = seq(min(dat$PC_3), max(dat$PC_3), length.out = pseudotime_bin_num)
dat$pseudotime_bin = cut(dat$PC_3,
                            breaks = pseudotime_bin_break,
                            include.lowest = T,
                            right = F) 

dat$day = factor(dat$day, levels = day_list[day_list %in% dat$day])
dat$pseudotime_bin = factor(dat$pseudotime_bin, levels = rev(levels(dat$pseudotime_bin)))

#pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
#x = as.vector(pd_all$day)
#x[pd_all$day == "E9"] = "E9.0"
#x[pd_all$day == "E10"] = "E10.0"
#x[pd_all$day == "E11"] = "E11.0"
#x[pd_all$day == "E12"] = "E12.0"
#x[pd_all$day == "E13"] = "E13.0"
#x[pd_all$day == "E14"] = "E14.0"
#x[pd_all$day == "E15"] = "E15.0"
#x[pd_all$day == "E16"] = "E16.0"
#x[pd_all$day == "E17"] = "E17.0"
#x[pd_all$day == "E18"] = "E18.0"
#pd_all$day = as.vector(x)
x = dat %>% group_by(day, pseudotime_bin) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day")
x$normalize_num = ceiling(100000 * x$n/x$total_n)
x$day = factor(x$day, levels = day_list[day_list %in% x$day])

library(reshape2)
dat = x %>% select(day, pseudotime_bin, normalize_num) %>%
    dcast(pseudotime_bin~day)
dat[is.na(dat)] = 0
rownames(dat) = dat[,1]
dat = dat[,-1]


library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("~/share/test.pdf", 8, 5)
heatmap.2(as.matrix(dat), 
          col=Colors, 
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

















############################# BACKUP #############################################




#################
### Re Jay's suggestion, 
### OK so right now you have a 1D timeline of distribution of cell types as a function of time. But what if there was a 2D heatmap where one axis was time and one was pseudotime, along one particular differentiation trajecotry (i.e. the pseudotime leading from retinal progenitors to retinal ganglion. Then the plot is basically a heat map of density (normalized in some way). there are very few P0 cells on the "bridge" of cells that leads to retinal ganglion cells, for example. This would be evident in the heatmap as the earlier pseuotime region only have earlier real-time cells.
### In contrast, photoreceptor precursors have tons of P0 cells, so the fact that that process of differentiation was still going at P0 would also be evident in its 2D heatmap
###  i don't know if we need to show all to get the concept across. it seems to me like if we could do a good job of the path to 8 (maybe would only need pseudotime for cluster 8 ) and 12 (maybe only need 12) and then to 5 to 6 and 5 to 7 (where i think the timing is a little different based on literature and panel c) that would get the key points across. trying to keep somewhat simple so we don't get bogged down, key is to illustrate the concept
### We could test out concept with 2D maps for 8, 5+6, 5+7 and 12


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
example_i = "Eye"
pd = readRDS(paste0(work_path, "/Eye/", example_i, "_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)
obj = readRDS(paste0(work_path, "/Eye/", example_i, "_obj.rds"))
print(table(pd$celltype_sub_clustering))
gene_count = GetAssayData(obj, slot = "counts")


######################################################
### First, focusing on Rod
######################################################

pd_sub = pd[pd$celltype_sub_cluster %in% c("Photoreceptor precursor cells", "Cone precursor cells"),]
gene_count_sub = gene_count[,rownames(pd_sub)]
obj_sub = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)
obj_sub = FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = 2500)
HVG = VariableFeatures(obj_sub)

cds_sub = doObjectTransform(obj_sub, transform_to = "monocle3")
cds_sub = preprocess_cds(cds_sub, num_dim = 30, use_genes = HVG)
cds_sub = reduce_dimension(cds_sub, umap.min_dist = 0.3, umap.n_neighbors = 30)
cds_sub = cluster_cells(cds_sub)
cds_sub = learn_graph(cds_sub)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin=c("E11.75", "E12.25", "E12.5", "E12.75")){
    cell_ids <- which(colData(cds)[, "day"] %in% time_bin)
    
    closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
}

cds_sub = order_cells(cds_sub, root_pr_nodes=get_earliest_principal_node(cds_sub))

saveRDS(cds_sub, paste0(work_path, "/Eye/", example_i, "_adata_scale.cds_Cone.rds"))

print(sum(colnames(cds_sub) == rownames(pd_sub)))
pd_sub$pseudotime = cds_sub@principal_graph_aux[["UMAP"]]$pseudotime
pd_sub_x = pd_sub %>% left_join(pd_num, by = "embryo_id") 
pd_sub$embryo_cell_num = as.vector(pd_sub_x$n)
pd_sub$embryo_cell_num_rev = 1/pd_sub$embryo_cell_num
day_x = as.vector(pd_sub$day)
day_x[pd_sub$day == "P0"] = "E19"
pd_sub$day_value = as.numeric(gsub("E","",day_x))

pd_sub = pd_sub[pd_sub$pseudotime != Inf,]

### manually calculate the heatmap
day_bin_num = 10

day_bin_break = seq(min(pd_sub$day_value), max(pd_sub$day_value), length.out = day_bin_num)
pd_sub$day_bin = cut(pd_sub$day_value,
                     breaks = day_bin_break,
                     include.lowest = T,
                     right = F) 

pseudotime_bin_num = 50

pseudotime_bin_break = seq(min(pd_sub$pseudotime), max(pd_sub$pseudotime), length.out = pseudotime_bin_num)
pd_sub$pseudotime_bin = cut(pd_sub$pseudotime,
                            breaks = pseudotime_bin_break,
                            include.lowest = T,
                            right = F) 

df = pd_sub %>% group_by(day_bin, pseudotime_bin) %>% 
    summarise(count_norm = sum(embryo_cell_num_rev))
df$density_norm = df$count_norm/sum(df$count_norm) 

saveRDS(df,  paste0(work_path, "/Eye/", example_i, "_adata_scale.df_Cone.rds"))


### locally checking if it looks good
setwd("~/work/jax/rna_seq/Eye/")
cds_sub = readRDS("Eye_adata_scale.cds_Cone.rds")
df = readRDS("Eye_adata_scale.df_Cone.rds")
plot_cells(cds_sub, color_cells_by = "day")
plot_cells(cds_sub, color_cells_by = "partition")
plot_cells(cds_sub, color_cells_by = "pseudotime")

# Bin size control + color palette
ggplot(df, aes(x=day_bin, y=pseudotime_bin, color = count_norm) ) +
    geom_point(size = 5) +
    scale_color_continuous(type = "viridis") +
    theme_bw()

