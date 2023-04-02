#########################################
### subclustering on renal trajectory ###
#########################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons"; print(example_i)

obj = readRDS(paste0(work_path, "/Neurons/", example_i, "_obj.rds"))

Idents(obj) = as.vector(obj$leiden_res_1)
obj_sub = subset(obj, downsample = 1000)
res = FindMarkers(obj_sub, ident.1 = c(21,6), ident.2 = c(9), only.pos = T)
res = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID")
head(res, 20)

pd = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)

celltype = rep(NA, nrow(pd))
celltype[pd$leiden_res_1 %in% c(2,18,23,25)] = "Spinal dI1 interneurons"
celltype[pd$leiden_res_1 %in% c(10,12)] = "Spinal dI2 interneurons"
celltype[pd$leiden_res_1 %in% c(27)] = "Spinal dI3 interneurons"
celltype[pd$leiden_res_1 %in% c(3,7,22,28)] = "Spinal dI4 interneurons"
celltype[pd$leiden_res_1 %in% c(0,8,14)] = "Spinal dI5 interneurons"
#celltype[pd$leiden_res_1 %in% c(4)] = "Spinal dI6 interneurons"  
#celltype[pd$leiden_res_1 %in% c(20)] = "Spinal V0 interneurons"     
#celltype[pd$leiden_res_1 %in% c(4,6,9,21)] = "Spinal V1 interneurons" 
celltype[pd$leiden_res_1 %in% c(11)] = "Spinal V2a interneurons"
celltype[pd$leiden_res_1 %in% c(19)] = "Spinal V2b interneurons"
celltype[pd$leiden_res_1 %in% c(29)] = "Spinal V3 interneurons"
celltype[pd$leiden_res_1 %in% c(5,15,17,24)] = "Di/mesencephalon glutamatergic neurons"
celltype[pd$leiden_res_1 %in% c(1)] = "Di/mesencephalon GABAergic neurons"
celltype[pd$leiden_res_1 %in% c(13,16)] = "Hypothalamic Sim1 neurons"
celltype[pd$leiden_res_1 %in% c(26)] = "Midbrain dopaminergic neurons"
celltype[pd$leiden_res_1 %in% c(31)] = "Spinal cord motor neuron progenitors"
celltype[pd$leiden_res_1 %in% c(30)] = "Precerebellar neurons"

celltype[pd$leiden_res_5 %in% c(50,38)] = "Striatal projection neurons"

pd$celltype_sub_clustering = as.vector(celltype)

pdx = readRDS(paste0(work_path, "/Neurons/Neurons_sub_1.rds"))
pd1 = pd[!is.na(celltype),]
pd2 = pd[is.na(celltype),]
pdx = pdx[as.vector(pd2$cell_id),]
pd2$celltype_sub_clustering = as.vector(pdx$celltype_sub_clustering)

pdx = rbind(pd1, pd2)
pdx = pdx[rownames(pd),]
pd$celltype_sub_clustering = as.vector(pdx$celltype_sub_clustering)

saveRDS(pd, paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_sub_clustering, colors = neuron_color_plate)
saveWidget(fig, paste0(save_path, "/example/", example_i, "_celltype_sub_clustering.html"), selfcontained = FALSE, libdir = "tmp")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.2) +
    theme_void() +
    scale_color_manual(values=neuron_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Neurons/", example_i, ".2D_UMAP.png"), width = 8, height = 6, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.2) +
    scale_color_manual(values=neuron_day_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Neurons/", example_i, ".day.2D_UMAP.png"), width = 8, height = 6, dpi = 300)



example_i = "Neurons"
pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

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

x = data.frame(day = names(neuron_day_color_plate), day_id = c(0:16))
pd = pd %>% left_join(x, by = "day")
pd$day_id = as.numeric((pd$day_id))

pd$celltype_sub_clustering = factor(pd$celltype_sub_clustering, levels = rev(names(neuron_color_plate)))

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

pd_sim$celltype_sub_clustering = factor(pd_sim$celltype_sub_clustering, levels = rev(names(neuron_color_plate)))

library(ggridges)
p=ggplot(pd_sim, aes(x = day_id, y = celltype_sub_clustering, fill = celltype_sub_clustering)) +
    geom_density_ridges() +
    scale_fill_manual(values=neuron_color_plate) +
    theme_ridges() + 
    labs(x = "", y = "") +
    theme(legend.position = "none")
pdf(paste0(work_path, "/Neurons/Neurons_density.pdf"),10,10)
p
dev.off()

df = pd %>% group_by(day) %>% sample_n(15) %>% as.data.frame()
df$day_x = as.numeric(as.vector(gsub("E","", as.vector(df$day))))
p = ggplot(df, aes(UMAP_2d_1, UMAP_2d_2, color = day_x)) + 
    geom_point(size = 3) + 
    theme_classic(base_size = 12) + 
    scale_color_viridis(option = "B", breaks = c(8.75, 9.75, 10.75, 11.75, 12.75), labels = c("E8.75","E9.75","E10.75","E11.75","E12.75")) + 
    labs(x = "UMAP_1", y = "UMAP_2") + 
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"), legend.key.height= unit(0.8, 'cm'))
pdf("../example/backup/Neurons_day_color_plate.pdf")
p
dev.off()

df = pd %>% group_by(celltype_sub_clustering) %>% sample_n(15) %>% as.data.frame()
p = df %>% ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=3) +
    theme_void() +
    scale_color_manual(values=neuron_color_plate)
pdf("../example/backup/Neurons_color_plate.pdf")
p
dev.off()


### Neurons marker genes 
my_plot_cells(cds2, genes = c("Otx2","Hoxb4","Hoxd4","Slc32a1","Gad1","Slc17a6"), how_many_rows = 2, cell_size = 0.5) + theme_void() + NoLegend() + theme(strip.text.x = element_blank()) + ggsave("./example/Neurons_marker_gene.png", width = 8, height = 4, dpi = 300)

x = my_plot_cells(cds2[,sample(1:ncol(cds2), 20000)], genes = c("Otx2","Hoxb4","Hoxd4","Slc32a1","Gad1","Slc17a6"))
pdf("./example/Neurons_marker_gene.pdf")
x
dev.off()

### Neurons marker genes 
my_plot_cells(cds2, genes = c("Eomes","Pax6","Satb2","Pou3f2","Pou3f3","Tbr1","Bcl11b","Fezf2","Kcnab1","Chrna5","Syt6","Foxp2"), how_many_rows = 3) + theme_void() + NoLegend() + theme(strip.text.x = element_blank()) + ggsave("./example/IP_marker_gene.png", width = 8, height = 6, dpi = 300)

x = my_plot_cells(cds2[,sample(1:ncol(cds2), 5000)], genes = c("Eomes","Pax6","Satb2","Pou3f2","Pou3f3","Tbr1","Bcl11b","Fezf2","Kcnab1","Chrna5","Syt6","Foxp2"))
pdf("./example/IP_marker_gene.pdf")
x
dev.off()

#################################################
### PCA analysis on those spinal interneurons ###
#################################################


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons"; print(example_i)
pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))
obj = readRDS(paste0(work_path, "/Neurons/", example_i, "_obj.rds"))
rownames(pd) = as.vector(pd$cell_id)
pd = pd[colnames(obj),]
obj$celltype_sub_clustering = pd$celltype_sub_clustering

### PCA analysis
obj$keep = obj$celltype_sub_clustering %in% c("Spinal dI1 interneurons",
                                              "Spinal dI2 interneurons",
                                              "Spinal dI3 interneurons",
                                              "Spinal dI4 interneurons",
                                              "Spinal dI5 interneurons",
                                              "Spinal dI6 interneurons",
                                              "Spinal V0 interneurons",
                                              "Spinal V1 interneurons",
                                              "Spinal V2a interneurons",
                                              "Spinal V2b interneurons",
                                              "Spinal V3 interneurons")
obj = subset(obj,subset = keep)
Idents(obj) = as.vector(obj$celltype_sub_clustering)
obj = subset(obj, downsample = 10000)

npcs = 30
reduction.key = "PC_"
seed.use = 42

obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)
obj = ScaleData(object = obj, verbose = FALSE, features = rownames(obj))
scale_dat = GetAssayData(obj, slot = "scale.data")
print(dim(scale_dat))

saveRDS(obj, paste0(work_path, "/Neurons/Neurons_PCA_obj.rds"))

set.seed(seed = seed.use)
pca.results <- irlba::irlba(A = t(x = scale_dat[genes_include,]), nv = npcs)
feature.loadings <- pca.results$v
set.seed(seed = seed.use)
cell.embeddings <- pca.results$u %*% diag(pca.results$d)

rownames(x = feature.loadings) <- genes_include
colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
rownames(x = cell.embeddings) <- colnames(obj)
colnames(x = cell.embeddings) <- colnames(x = feature.loadings)

stdev <- pca.results$d/sqrt(max(1, ncol(scale_dat) - 1))
eigValues = (stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)

res = list(cell.embeddings = cell.embeddings,
           feature.loadings = feature.loadings,
           varExplained = varExplained,
           df_cell = data.frame(obj[[]]))
saveRDS(res, paste0(work_path, "/Neurons/Neurons_PCA_res.rds"))


###### locally ######

res = readRDS("~/work/jax/rna_seq/Neurons/Neurons_PCA_res.rds")
df = data.frame(res[["df_cell"]], res[["cell.embeddings"]])


x = as.vector(df$day)
x[df$day == "E9"] = "E9.0"
x[df$day == "E10"] = "E10.0"
x[df$day == "E11"] = "E11.0"
x[df$day == "E12"] = "E12.0"
x[df$day == "E13"] = "E13.0"
x[df$day == "E14"] = "E14.0"
x[df$day == "E15"] = "E15.0"
x[df$day == "E16"] = "E16.0"
x[df$day == "E17"] = "E17.0"
x[df$day == "E18"] = "E18.0"
df$day = as.vector(x)


x = gsub("Spinal ", "", df$celltype_sub_clustering)
x = gsub(" interneurons", "", x)
df$celltype = x

df$celltype = factor(df$celltype, levels = c("dI1","dI2","dI3","dI4","dI5","dI6",
                                             "V0","V1","V2a","V2b","V3"))
df$day = factor(df$day, levels = names(neuron_day_color_plate))

neuron_sub_color_plate = c("dI1" = '#636EFA',
                       "dI2" = '#EF553B',
                       "dI3" = '#00CC96',
                       "dI4" = '#AB63FA',
                       "dI5" = '#FFA15A',
                       "dI6" = '#19D3F3',
                       "V0" = '#FF6692',
                       "V1" = '#B6E880',
                       "V2a" = '#FECB52',
                       "V2b" = '#1F77B4',
                       "V3" = '#FF7F0E')

# making the 3D-PCA or the boxplot for individual PCs

t1 = list(family = 'Helvetica',
          size = 25,
          color = "black")
t2 = list(family = 'Helvetica',
          size = 15,
          color = "grey")
fig1 = plot_ly(df, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~celltype, colors = neuron_sub_color_plate, size = 1) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (31.1%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (8.6%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (5.7%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))),
           showlegend = FALSE)

x_table = table(df$day)
df_1 = df[df$day %in% names(x_table)[x_table <= 1000],]
df_2 = df[df$day %in% names(x_table)[x_table > 1000],] %>% group_by(day) %>% sample_n(1000)
df_sub = rbind(df_1, df_2)
fig2 = plot_ly(df_sub, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~day, colors = neuron_day_color_plate, size = 1) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (31.1%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (8.6%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (5.7%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))),
           showlegend = FALSE)

library(gridExtra) 

p1 = df %>%
    ggplot( aes(day, PC_4, fill = day)) + 
    geom_boxplot() + 
    labs(x="", y="PC_4 (5.0%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_day_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

p2 = df %>%
    ggplot( aes(celltype, PC_4, fill = celltype)) + 
    geom_boxplot() + 
    labs(x="", y="PC_4 (5.0%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_sub_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

pdf("~/work/jax/figures/tmp/Neurons_PCA_4.pdf", 10, 5)
grid.arrange(p1, p2, nrow=1, ncol=2) 
dev.off()

p = df %>% 
    sample_n(20000) %>%
    ggplot(aes(PC_1, log2_umi, color = day)) + 
    geom_point() + 
    labs(x="PC_1 (31.1%)", y="Log2_umi", title="") +
    theme_classic(base_size = 20) +
    scale_color_manual(values=neuron_day_color_plate) + 
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    NoLegend()
pdf("~/work/jax/figures/tmp/Neurons_PCA_1_umi.pdf", 8, 5)
p
dev.off()

p + theme_void() + NoLegend() +
    ggsave("~/work/jax/figures/tmp/Neurons_PCA_1_umi.png", width = 9.727149, height = 5, dpi = 300)


# calculate the correlation of UMI and day with each PCs
day_id = c(1:length(neuron_day_color_plate))
names(day_id) = names(neuron_day_color_plate)
df$day_id = as.numeric(day_id[as.vector(df$day)])

res_day = NULL
res_log2_umi = NULL
for(i in 1:30){
    fit = cor.test(df[[paste0("PC_", i)]], df$day_id)
    res_day = rbind(res_day, data.frame(PC = paste0("PC_", i), corr = fit$estimate, pval = fit$p.value))
    fit = cor.test(df[[paste0("PC_", i)]], df$log2_umi)
    res_log2_umi = rbind(res_log2_umi, data.frame(PC = paste0("PC_", i), corr = fit$estimate, pval = fit$p.value))
}

write.csv(res_day, "~/work/jax/rna_seq/Neurons/Neurons_PCA_correlation_with_day.csv", row.names=F)
write.csv(res_log2_umi, "~/work/jax/rna_seq/Neurons/Neurons_PCA_correlation_with_log2_umi.csv", row.names=F)



# calculate the correlation of genes with each PC after downsamling to 20K cells

source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data/example"

res = readRDS(paste0(work_path, "/Neurons/Neurons_PCA_res.rds"))
obj = readRDS(paste0(work_path, "/Neurons/Neurons_PCA_obj.rds"))
pd_x = res[["df_cell"]]
emb = res[["cell.embeddings"]]

gene_count = GetAssayData(obj, slot = "counts")
gene_count_x = gene_count[,rownames(emb)]
obj_x = CreateSeuratObject(gene_count_x, meta.data = pd_x)

obj_x = NormalizeData(obj_x, normalization.method = "LogNormalize", scale.factor = 10000)
obj_x = FindVariableFeatures(obj_x, selection.method = "vst", nfeatures = 2500)
obj_x = ScaleData(obj_x, verbose = FALSE)
scale_dat = GetAssayData(obj_x, slot = "scale.data")

if(nrow(emb) > 20000){
    keep = sample(1:nrow(emb), 20000)
} else {
    keep = 1:nrow(emb)
}

emb_sub = emb[keep,]
scale_dat_sub = scale_dat[,keep]

corr_matrix = matrix(NA, 2500, 4)
pval_matrix = matrix(NA, 2500, 4)
for(j in 1:4){
    for(k in 1:2500){
        print(k)
        fit = cor.test(emb_sub[,j], scale_dat_sub[k,])
        corr_matrix[k,j] = fit$estimate
        pval_matrix[k,j] = fit$p.value
    }
}

mouse_gene_sub = mouse_gene[rownames(scale_dat),]

df = data.frame(gene_ID = rep(rownames(scale_dat), 4),
                corr = c(corr_matrix[,1], corr_matrix[,2], corr_matrix[,3], corr_matrix[,4]),
                pval = c(pval_matrix[,1], pval_matrix[,2], pval_matrix[,3], pval_matrix[,4]),
                PC = rep(c("PC_1", "PC_2", "PC_3", "PC_4"), each = 2500),
                gene_short_name = rep(as.vector(mouse_gene_sub$gene_short_name), 4), stringsAsFactors = F)

saveRDS(df,  paste0(work_path, "/Neurons/Neurons_PCA_corr_genes.rds"))

df_sub = df %>% filter(!is.na(pval), PC == "PC_1") %>%
    mutate(qval = p.adjust(pval, method = "fdr")) 
x1 = mean(df_sub$corr) - sd(df_sub$corr)
x2 = mean(df_sub$corr) + sd(df_sub$corr)
df_1 = df_sub %>% filter(qval < 0.05) %>%
    filter(corr < x1 | corr > x2) %>%
    arrange(corr) %>%
    select(-pval)

df_sub = df %>% filter(!is.na(pval), PC == "PC_2") %>%
    mutate(qval = p.adjust(pval, method = "fdr")) 
x1 = mean(df_sub$corr) - sd(df_sub$corr)
x2 = mean(df_sub$corr) + sd(df_sub$corr)
df_2 = df_sub %>% filter(qval < 0.05) %>%
    filter(corr < x1 | corr > x2) %>%
    arrange(corr) %>%
    select(-pval)

df_sub = df %>% filter(!is.na(pval), PC == "PC_3") %>%
    mutate(qval = p.adjust(pval, method = "fdr")) 
x1 = mean(df_sub$corr) - sd(df_sub$corr)
x2 = mean(df_sub$corr) + sd(df_sub$corr)
df_3 = df_sub %>% filter(qval < 0.05) %>%
    filter(corr < x1 | corr > x2) %>%
    arrange(corr) %>%
    select(-pval)

df_sub = df %>% filter(!is.na(pval), PC == "PC_4") %>%
    mutate(qval = p.adjust(pval, method = "fdr")) 
x1 = mean(df_sub$corr) - sd(df_sub$corr)
x2 = mean(df_sub$corr) + sd(df_sub$corr)
df_4 = df_sub %>% filter(qval < 0.05) %>%
    filter(corr < x1 | corr > x2) %>%
    arrange(corr) %>%
    select(-pval)

df_out = rbind(df_1, df_2, df_3, df_4)

write.csv(df_out,  paste0(work_path, "/Neurons/Neurons_PCA_corr_genes_significant.csv"))










#############################################################
### for each interneuron, what are the top expressed TFs? ### 
#############################################################


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

obj = readRDS(paste0(work_path, "/Neurons/Neurons_PCA_obj.rds"))
obj = CreateSeuratObject(GetAssayData(obj, slot = "counts"), meta = data.frame(obj[[]]))
Idents(obj) = as.vector(obj$celltype_sub_clustering)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

mouse_TF = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/TF_list/AnimalTFDB/Mus_musculus_TF.txt", sep="\t", header=T, as.is=T)
mouse_TF_overlap = intersect(as.vector(mouse_TF$Ensembl), rownames(obj))

result = FindAllMarkers(obj, features = mouse_TF_overlap)
rownames(result) = NULL
result$gene_ID = result$gene
x = result %>% left_join(mouse_gene, by = "gene_ID")
result$gene_short_name = as.vector(x$gene_short_name)
result$gene = NULL
write.csv(result, paste0(work_path, "/Neurons/Neurons_interneurons_top_TFs.csv"), row.names=F)

result_sub = result %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 3)

data = GetAssayData(obj, slot = "data")

celltype_list = names(neuron_color_plate)[1:11]
result_sub$cluster = factor(result_sub$cluster, levels = celltype_list)
result_sub = result_sub %>% arrange(cluster)

exp = NULL
for(i in celltype_list){
    print(i)
    exp = cbind(exp, Matrix::rowMeans(data[,obj$celltype_sub_clustering == i]))
}

x = gsub("Spinal ", "", celltype_list)
x = gsub(" interneurons", "", x)
colnames(exp) = x

gene_list = unique(as.vector(result_sub$gene_short_name))
mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% gene_list,]
rownames(mouse_gene_sub) = mouse_gene_sub$gene_short_name
mouse_gene_sub = mouse_gene_sub[gene_list,]
exp_sub = exp[as.vector(mouse_gene_sub$gene_ID),]
rownames(exp_sub) = rownames(mouse_gene_sub)

saveRDS(exp_sub, paste0(work_path, "/Neurons/Neurons_heatmap_dat.rds"))

dat = readRDS("~/work/jax/rna_seq/Neurons/Neurons_heatmap_dat.rds")

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("~/work/jax/figures/tmp/Neurons_heatmap.pdf", 8, 5)
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


################################################################################################
### Visualizing CNS neurons trajectory and Intermediate progenitor trajectory in the 3D UMAP ###
################################################################################################

source("~/work/scripts/tome/utils.R")

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

pd = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))

pd$tmp = pd$major_trajectory
pd$tmp[!pd$major_trajectory %in% c("Intermediate_progenitors", "CNS_neurons", "Brain_and_spinal_cord" )] = "other"

tmp_color_plate = c("Brain_and_spinal_cord"     = "#f96100",
                    "Intermediate_progenitors"  = "#2e0ab7",
                    "CNS_neurons"               = "#e5c000",
                    "other"                     = "#B3B3B3")

pd_sub = pd[sample(1:nrow(pd),300000),]

sub1 = pd_sub[!pd_sub$tmp %in% c("other"),]
sub2 = pd_sub[pd_sub$tmp %in% c("other"),]

t1 = list(family = 'Helvetica',
          size = 25,
          color = "black")
t2 = list(family = 'Helvetica',
          size = 15,
          color = "grey")

fig <- plotly::plot_ly(sub1) %>%
    plotly::add_trace(x = ~UMAP_1, 
                      y = ~UMAP_2, 
                      z = ~UMAP_3,
                      type = 'scatter3d', 
                      size=I(30), 
                      alpha = I(1),
                      color = ~tmp,
                      colors = tmp_color_plate,
                      mode="markers",
                      showlegend=FALSE) %>%
    plotly::add_markers(x = sub2$UMAP_1, 
                        y = sub2$UMAP_2,
                        z = sub2$UMAP_3, 
                        color = I("lightgray"),
                        size=I(30),
                        marker=list(opacity = 0.1), 
                        showlegend=FALSE)  %>% layout(
                            scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                                         yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                                         zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                                         camera = list(eye = list(x = 1, y = 1, z = 1)))) %>% hide_colorbar()

saveWidget(fig, paste0(save_path, "/whole_dataset/whole_dataset_scale_major_trajectory_highlight_neurons.html"), selfcontained = FALSE, libdir = "tmp")


### making an easy plot to show the time series between "direct" and "indirect" trajectories

pd$day = factor(pd$day, levels = day_list)
pd_1 = pd[pd$major_trajectory == "CNS_neurons",]
pd_2 = pd[pd$major_trajectory == "Intermediate_progenitors",]
pd_3 = pd[pd$major_trajectory == "Intermediate_progenitors" & pd$celltype_update == "Intermediate progenitors",]
pd_4 = pd[pd$major_trajectory == "Intermediate_progenitors" & pd$celltype_update != "Intermediate progenitors",]

x1 = pd_1 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x2 = pd_2 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x = x1 %>% select(day, frac) %>% rename(direct_frac = frac) %>% left_join(x2 %>% select(day, frac) %>% rename(indirect_frac = frac), by = "day")
x$indirect_frac[is.na(x$indirect_frac)] = 0
x = data.frame(day = rep(x$day, 2),
               frac = c(x$direct_frac, x$indirect_frac),
               major_trajectory = rep(c("CNS_neurons","Intermediate_progenitors"), each = nrow(x)))

p = x %>% 
    ggplot(aes(x=day, y=frac, fill = day)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(major_trajectory)) + 
    labs(x='',y='% of cells') +
    scale_fill_manual(values=day_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0(work_path, "/Neurons/Neurons_density_1.pdf"), 6, 3)
p
dev.off()


####################################
### subclustering on dI1 neurons ###
####################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons"; print(example_i)
pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))
obj = readRDS(paste0(work_path, "/Neurons/", example_i, "_obj.rds"))
rownames(pd) = as.vector(pd$cell_id)
pd = pd[colnames(obj),]
obj$celltype_sub_clustering = pd$celltype_sub_clustering
obj = subset(obj,subset = celltype_sub_clustering == "Spinal dI1 interneurons")

obj$group = "Spinal dI1 interneurons"
obj_processed = doClusterSeurat(obj)
saveRDS(obj_processed, paste0(work_path, "/Neurons/dI1_clustering_obj.rds"))

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds)

reducedDims(cds)$UMAP = 
    as.matrix(Embeddings(obj_processed, reduction = "umap"))
saveRDS(cds, paste0(work_path, "/Neurons/dI1_clustering_cds.rds"))


obj = readRDS("~/work/jax/rna_seq/Neurons/dI1_clustering_obj.rds")
Idents(obj) = as.vector(obj$leiden_res_1)
obj_sub = subset(obj, downsample = 2000)
res = FindMarkers(obj_sub, ident.1 = c(18), ident.2 = c(23,25), only.pos = T)
res = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID")

anno = rep(NA, ncol(obj))
anno[obj$leiden_res_1 == 2] = "Progenitor cells"
anno[obj$leiden_res_1 == 25] = "dI1 interneurons (Ndst3+)"
anno[obj$leiden_res_1 == 18] = "dI1 interneurons (Vcan+)"
anno[obj$leiden_res_1 == 23] = "dI1 interneurons (Hoxc6+)"

obj$anno = anno
pd = data.frame(day = obj$day, 
                anno = obj$anno,
                UMAP_2d_1 = Embeddings(obj, reduction = "umap")[,1],
                UMAP_2d_2 = Embeddings(obj, reduction = "umap")[,2])

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

pd$day = factor(pd$day, levels = names(neuron_day_color_plate))


p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = anno), size=0.6) +
    theme_void() +
    scale_color_manual(values=dI1_neuron_color_plate) +
    theme(legend.position="none") + 
    ggsave("./example/dI1_neurons.2D_UMAP.png", width = 5, height = 4, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.6) +
    scale_color_manual(values=neuron_day_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("./example/dI1_neurons.day.2D_UMAP.png", width = 5, height = 4, dpi = 300)

plot_cells(cds, genes = c("Fras1","Ndst3","Hs3st4","Gpc5","Hoxc6","Vcan","Mki67","Pax3","Slc17a6"), cell_size = 0.5) +
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) + 
    ggsave("./example/dI1_neurons_marker_genes.png", width = 10, height = 8, dpi = 300)

plot_cells(cds[,sample(1:ncol(cds),10000)], genes = c("Fras1","Ndst3","Hs3st4","Gpc5","Hoxc6","Vcan","Mki67","Pax3","Slc17a6"), cell_size = 0.5) +
    ggsave("./example/dI1_neurons_marker_genes.pdf", width = 10, height = 8, dpi = 300)





############################################################################
### Remaking the UMAP of intermediate neurons with different time series ###
############################################################################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

trajectory_i = "Intermediate_progenitors"
pd = readRDS(paste0(work_path, "/mtx/sub_clustering/", trajectory_i, "_pd.rds"))

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


x_table = table(pd$day_group)
pd_1 = pd[pd$day_group %in% names(x_table)[x_table <= 10000],]
pd_2 = pd[pd$day_group %in% names(x_table)[x_table > 10000],]
pd_2_x = pd_2 %>% group_by(day_group) %>% sample_n(10000)
pd_sub = rbind(pd_1, pd_2[pd_2$cell_id %in% pd_2_x$cell_id,])

library(RColorBrewer)
day_group_color_plate=rev(brewer.pal(11,"Spectral"))
day_group_color_plate=colorRampPalette(day_group_color_plate)(8)
names(day_group_color_plate) = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                 "E16-E17", "E17-E18", "E18-P0", "P0")

p = pd_sub[sample(1:nrow(pd_sub)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.3) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "/Neurons/", trajectory_i, ".2D_UMAP.day_group_downsample.png"), width = 6, height = 6, dpi = 300)






############################# Backup ############################################

#############################################
### subclustering the spinal interneurons ###

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

example_i = "Neurons"; print(example_i)

pd = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.csv"), header=T, row.names=1)

pd_sub = pd[pd$leiden_res_1 %in% c(4,6,9,20,21),]
pdx = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
pdx = pdx[rownames(pd_sub),]
pd_sub$embryo_id = as.vector(pdx$embryo_id)


### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
fd = readRDS("/net/gs/vol1/home/cxqiu/work/jax/rna_seq/df_gene.rds")
fd = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & fd$chr %in% paste0("chr", c(1:19, "M")),]


embryo_list = as.vector(names(table(pd_sub$embryo_id)))

gene_count = NULL
for(i in 1:length(embryo_list)){
    embryo_i = embryo_list[i]; print(paste0(i, "/", embryo_i))
    gene_count_i = readRDS(paste0(work_path, "/embryo/", embryo_i, "_gene_count.rds"))
    gene_count_i = gene_count_i[rownames(gene_count_i) %in% rownames(fd),
                                colnames(gene_count_i) %in% as.vector(pd_sub$cell_id)]
    print(dim(gene_count_i))
    gene_count = cbind(gene_count, gene_count_i)
}

rownames(pd_sub) = as.vector(pd_sub$cell_id)
gene_count = gene_count[,rownames(pd_sub)]
pd_sub$leiden = pd_sub$UMAP_1 = pd_sub$UMAP_2 = pd_sub$UMAP_3 = pd_sub$UMAP_2d_1 = pd_sub$UMAP_2d_2 = NULL

obj = CreateSeuratObject(gene_count, meta.data = pd_sub)
saveRDS(obj, paste0(work_path, "/Neurons/Neurons_sub_obj.rds"))

obj$group = "Neurons_sub"
obj_processed = doClusterSeurat(obj)
saveRDS(obj_processed, paste0(work_path, "/Neurons/Neurons_sub_obj_processed.rds"))

obj$keep = 1:ncol(obj) %in% sample(1:ncol(obj),50000)
cds = doObjectTransform(subset(obj, subset = keep), transform_to = "monocle3")
cds = preprocess_cds(cds)

reducedDims(cds)$UMAP = 
    as.matrix(Embeddings(obj_processed, reduction = "umap")[obj$keep,])
saveRDS(cds, paste0(work_path, "/Neurons/Neurons_sub_cds.rds"))

obj_processed = readRDS(paste0(work_path, "/Neurons/Neurons_sub_obj_processed.rds"))
pd = data.frame(obj_processed[[]])

celltype = rep("Spinal V1 interneurons", nrow(pd))
celltype[pd$RNA_snn_res.1 %in% c(2,8,11,19)] = "Spinal V0 interneurons"    
celltype[pd$RNA_snn_res.1 %in% c(10,12,13,17,18,20,22)] = "Spinal dI6 interneurons"    
celltype[pd$RNA_snn_res.1 %in% c(21) & pd$leiden_res_5 == 1] = "Spinal V2b interneurons"  
celltype[pd$leiden_res_1 %in% c(20)] = "Spinal V0 interneurons"   
pd$celltype_sub_clustering = celltype
saveRDS(pd, paste0(work_path, "/Neurons/Neurons_sub_1.rds"))





