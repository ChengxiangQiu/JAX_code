#######################################
### subclustering on LPM trajectory ###
#######################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "LPM"; print(example_i)

obj = readRDS(paste0(work_path, "/mtx/example/", example_i, "_obj.rds"))

Idents(obj) = as.vector(obj$leiden_res_1)
obj_sub = subset(obj, downsample = 500)
res = FindMarkers(obj_sub, ident.1 = c(20), only.pos = T)
res = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID")
head(res, 20)

count = GetAssayData(obj, slot = "counts")
count_x = count["ENSMUSG00000021765",as.vector(pd$cell_id)]

pd = read.csv(paste0(work_path, "/mtx/example/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
pd_all = pd_all[rownames(pd),]
pd$embryo_sex = as.vector(pd_all$embryo_sex)
pd$Fst = as.vector(count["ENSMUSG00000021765",as.vector(pd$cell_id)])

fig = plot_ly(pd[sample(1:nrow(pd), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~embryo_sex)
saveWidget(fig, paste0(save_path, "/example/", example_i, "_embryo_sex.html"), selfcontained = FALSE, libdir = "tmp")


celltype = rep(NA, nrow(pd))
celltype[pd$leiden_res_1 %in% c(3,13,26)] = "Lung mesenchyme"
celltype[pd$leiden_res_1 %in% c(21)] = "Hepatic mesenchyme"
celltype[pd$leiden_res_1 %in% c(16)] = "Renal stromal cells"
celltype[pd$leiden_res_1 %in% c(17)] = "Cardiopharyngeal mesoderm (Tbx1+)"
celltype[pd$leiden_res_1 %in% c(5)] = "Vascular smooth muscle cells"
celltype[pd$leiden_res_1 %in% c(30)] = "Meninges"
celltype[pd$leiden_res_1 %in% c(28)] = "Airway smooth muscle cells"
celltype[pd$leiden_res_1 %in% c(12)] = "Amniotic mesoderm"
celltype[pd$leiden_res_1 %in% c(7)] = "Allantois"
celltype[pd$leiden_res_1 %in% c(1)] = "Extraembryonic mesoderm"
celltype[pd$leiden_res_1 %in% c(2,6,11,14,23)] = "Splanchnic mesoderm"
celltype[pd$leiden_res_1 %in% c(0,24,19,4)] = "Somatic mesoderm"
celltype[pd$leiden_res_1 %in% c(20)] = "Renal capsule"
celltype[pd$leiden_res_1 %in% c(8,9,22)] = "Gut mesenchyme"
celltype[pd$leiden_res_1 %in% c(18)] = "Foregut mesenchyme"
celltype[pd$leiden_res_1 %in% c(29)] = "Gonad progenitor cells"
celltype[pd$leiden_res_1 %in% c(25)] = "Proepicardium (Tbx18+)"
celltype[pd$leiden_res_1 %in% c(27,10,15)] = "Mesothelial cells"

celltype[pd$leiden_res_5 %in% c(18)] = "Vascular smooth muscle cells (Pparg+)"
celltype[pd$leiden_res_5 %in% c(14,12,81)] = "Hepatic mesenchyme"
celltype[pd$leiden_res_5 %in% c(13)] = "Gastrointestinal smooth muscle cells"
celltype[pd$leiden_res_5 %in% c(93)] = "Sertoli cells"

celltype[pd$leiden_res_1 == 29 &
         pd$leiden_res_5 != 93 & 
         pd$embryo_sex == "F" &
         pd$Fst > 0 &
         pd$celltype == "Granulosa cells"] = "Granulosa cells"

pd$celltype_sub_clustering = as.vector(celltype)

saveRDS(pd, paste0(work_path, "/mtx/example/", example_i, "_adata_scale.obs.rds"))

fig = plot_ly(pd[sample(1:nrow(pd), 200000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_sub_clustering, colors = LPM_color_plate)
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_sub_clustering.html"), selfcontained = FALSE, libdir = "tmp")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.2) +
    theme_void() +
    scale_color_manual(values=LPM_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/mtx/example/", example_i, ".2D_UMAP.png"), width = 8, height = 6, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.2) +
    scale_color_manual(values=day_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/mtx/example/", example_i, ".day.2D_UMAP.png"), width = 8, height = 6, dpi = 300)





############ BACKUP ###############

example_i = "LPM"
pd = readRDS(paste0(work_path, "/mtx/example/", example_i, "_adata_scale.obs.rds"))

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

pd$celltype_sub_clustering = factor(pd$celltype_sub_clustering, levels = rev(names(LPM_color_plate)))

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

pd_sim$celltype_sub_clustering = factor(pd_sim$celltype_sub_clustering, levels = rev(names(LPM_color_plate)))

library(ggridges)
p=ggplot(pd_sim, aes(x = day_id, y = celltype_sub_clustering, fill = celltype_sub_clustering)) +
    geom_density_ridges() +
    scale_fill_manual(values=LPM_color_plate) +
    theme_ridges() + 
    labs(x = "", y = "") +
    theme(legend.position = "none")
pdf(paste0(work_path, "/mtx/example/LPM_density.pdf"),10,10)
p
dev.off()

df = pd %>% group_by(celltype_sub_clustering) %>% sample_n(15) %>% as.data.frame()
p = df %>% ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=3) +
    theme_void() +
    scale_color_manual(values=LPM_color_plate)
pdf(paste0(work_path, "/mtx/example/LPM_color_plate.pdf"),10,10)
p
dev.off()

####################################
### making gene expression plots ###
####################################

setwd("~/work/jax/rna_seq/LPM")
cds = readRDS("LPM_cds2.rds")

# lung mesenchyme (Tbx5+, Tbx4+) (Arora, Metzger, and Papaioannou 2012), 
# hepatic mesenchyme (Reln+) (Dobie et al. 2019), 
# gut mesenchyme (Nkx2-3+) (Pabst et al. 1997), 
# foregut mesenchyme (Barx1+) (Jayewickreme and Shivdasani 2015), 
# Amniotic mesoderm (Postn)

# renal stromal cells (Foxd1+, Tcf21+) (Finer et al. 2022). 
# renal capsule (Pax2+, Pax8+) (Bouchard et al. 2002)

# Second, some smooth muscle cells (Acta2+, Myh11+) show organ-specific pattern, including 
# airway smooth muscle cells (Trpc6+, Tbx5+) (Godin and Rousseau 2007), 
# gastrointestinal smooth muscle cells (Nkx2-3+) (Pabst et al. 1997). 

# Next, some cell types were highly leveraged spatial information to better annotate, including 
# proepicardium (Tbx18+) (Quijada et al. 2021), 
# mesothelial cells (Msln+) (Hassan and Ho 2008), 
# Meninges (Vtn+) (Seiffert et al. 1995). 



my_plot_cells(cds, genes = c("Tbx5","Tbx4","Reln","Nkx2-3", "Barx1", "Postn"), how_many_rows = 1) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("LPM_gene_exp_1.png", dpi = 300, height = 2, width = 16)

my_plot_cells(cds, genes = c("Foxd1","Tcf21","Pax2","Pax8","Vtn","Mmp9"), how_many_rows = 1) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("LPM_gene_exp_2.png", dpi = 300, height = 2, width = 16)

my_plot_cells(cds, genes = c("Acta2","Myh11","Trpc6","Tbx18","Upk3b","Msln"), how_many_rows = 1) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("LPM_gene_exp_3.png", dpi = 300, height = 2, width = 16)


p = my_plot_cells(cds[,sample(1:ncol(cds),5000)], genes = c("Wnt11","Ret","Etv5","Wnt7b","Tacstd2"), how_many_rows = 1)
pdf("~/work/jax/rna_seq/Renal/test.pdf")
p
dev.off()





### output data matrix in 10X format, for Jay and Beth ###

obj = readRDS(paste0(work_path, "/mtx/example/LPM_obj.rds"))
pd = readRDS(paste0(work_path, "/mtx/example/LPM_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)
pd = pd[colnames(obj),]
count = GetAssayData(obj, slot = "counts")

if(ncol(count) > 100000){
    index = sample(1:ncol(count), 100000)
    pd = pd[index,]
    count = count[,index]
}

gene = mouse_gene[rownames(count),]
gene$mean_exp = Matrix::rowMeans(count)
gene_x = gene %>% group_by(gene_short_name) %>% slice_max(order_by = mean_exp, n = 1) 
gene_1 = gene[rownames(gene) %in% gene_x$gene_ID,]
gene_2 = gene[!rownames(gene) %in% gene_x$gene_ID,]
gene_2$gene_short_name = paste0(gene_2$gene_short_name, ".dup")
gene = rbind(gene_1, gene_2)
gene = gene[rownames(count),]

rownames(count) = as.vector(gene$gene_short_name)

save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data/share/LPM"

Matrix::writeMM(count, paste0(save_path, "/matrix.mtx"))
write.table(rownames(count), paste0(save_path, "/genes.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(colnames(count), paste0(save_path, "/barcodes.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(as.vector(pd$celltype_sub_clustering), paste0(save_path, "/celltype.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(as.vector(pd$day), paste0(save_path, "/day.tsv"), row.names=F, col.names=F, sep="\t", quote=F)


pd_sub = pd[,!colnames(pd) %in% c("celltype","keep", "leiden", "celltype_update", "celltype_sub_clustering", "UMAP_1", 
                                  "UMAP_2", "UMAP_3","UMAP_2d_1","UMAP_2d_2", "group", "batch", "cell_id")]
pd_sub$celltye = pd$celltype_sub_clustering
write.table(pd_sub, paste0(save_path, "/meta.tsv"), row.names=T, col.names=T, sep="\t", quote=F)






