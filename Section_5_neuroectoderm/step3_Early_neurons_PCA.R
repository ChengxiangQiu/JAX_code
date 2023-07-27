
##################################
### Section - 5, Neuroectoderm ###
##################################

###############################################
### PCA analysis on the spinal interneurons ###
###############################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Neurons"; print(example_i)
pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))
rownames(pd) = as.vector(pd$cell_id)

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count = doExtractData(pd, mouse_gene_sub)
obj = CreateSeuratObject(gene_count, pd)

### PCA analysis

obj$keep = obj$neurons_sub_clustering %in% c("Spinal dI1 interneurons",
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
Idents(obj) = as.vector(obj$neurons_sub_clustering)
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

saveRDS(obj, paste0(work_path, "Neurons_PCA_obj.rds"))
saveRDS(res, paste0(work_path, "Neurons_PCA_res.rds"))


################################
### Plotting PC1 x PC2 x PC3 ###
################################

df = data.frame(res[["df_cell"]], res[["cell.embeddings"]])

x = gsub("Spinal ", "", df$neurons_sub_clustering)
x = gsub(" interneurons", "", x)
df$celltype = x

df$celltype = factor(df$celltype, levels = c("dI1","dI2","dI3","dI4","dI5","dI6",
                                             "V0","V1","V2a","V2b","V3"))
day_list = names(day_color_plate)
df$day = factor(df$day, levels = day_list[day_list %in% df$day])

neuron_sub_color_plate = c("dI1" = '#636EFA',
                           "dI2" = '#EF553B',
                           "dI3" = '#00CC96',
                           "dI4" = '#AB63FA',
                           "dI5" = '#FFA15A',
                           "dI6" = '#19D3F3',
                           "V0"  = '#FF6692',
                           "V1"  = '#B6E880',
                           "V2a" = '#FECB52',
                           "V2b" = '#1F77B4',
                           "V3"  = '#FF7F0E')

### Fig. 5g

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

### Fig. 5h

p1 = df %>%
    ggplot( aes(day, PC_1, fill = day)) + 
    geom_boxplot() + 
    labs(x="", y="PC_1 (31.1%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_day_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

p2 = df %>%
    ggplot( aes(celltype, PC_1, fill = celltype)) + 
    geom_boxplot() + 
    labs(x="", y="PC_1 (31.1%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_sub_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

pdf(paste0(work_path, "Neurons_PCA_1.pdf"), 10, 5)
grid.arrange(p1, p2, nrow=1, ncol=2) 
dev.off()


p1 = df %>%
    ggplot( aes(day, PC_2, fill = day)) + 
    geom_boxplot() + 
    labs(x="", y="PC_2 (8.6%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_day_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

p2 = df %>%
    ggplot( aes(celltype, PC_2, fill = celltype)) + 
    geom_boxplot() + 
    labs(x="", y="PC_2 (8.6%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_sub_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

pdf(paste0(work_path, "Neurons_PCA_2.pdf"), 10, 5)
grid.arrange(p1, p2, nrow=1, ncol=2) 
dev.off()


p1 = df %>%
    ggplot( aes(day, PC_3, fill = day)) + 
    geom_boxplot() + 
    labs(x="", y="PC_3 (5.7%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_day_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

p2 = df %>%
    ggplot( aes(celltype, PC_3, fill = celltype)) + 
    geom_boxplot() + 
    labs(x="", y="PC_3 (5.7%)", title="") +
    theme_classic(base_size = 20) +
    scale_fill_manual(values=neuron_sub_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()

pdf(paste0(work_path, "Neurons_PCA_3.pdf"), 10, 5)
grid.arrange(p1, p2, nrow=1, ncol=2) 
dev.off()


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

pdf(paste0(work_path, "Neurons_PCA_4.pdf"), 10, 5)
grid.arrange(p1, p2, nrow=1, ncol=2) 
dev.off()


#############################################################
### calculate the correlation of genes with each PC after ###
#############################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

res = readRDS(paste0(work_path, "Neurons_PCA_res.rds"))
obj = readRDS(paste0(work_path, "Neurons_PCA_obj.rds"))
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

saveRDS(df,  paste0(work_path, "Neurons_PCA_corr_genes.rds"))

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

write.csv(df_out,  paste0(work_path, "Neurons_PCA_corr_genes_significant.csv"))




