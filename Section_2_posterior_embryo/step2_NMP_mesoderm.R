
#####################################
### Section - 2, Posterior embryo ###
#####################################

########################
### Analysis on NMPs ###
########################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "posterior_embryo"

i = "NMP_Mesoderm"

pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)

### 2D UMAP of NMPs, with cells are colored by their initial cell types (Fig. 2c)
p = pd_x %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.35) +
    theme_void() +
    scale_color_manual(values=posterior_embryo_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "NMPs_celltype.png"), width = 4, height = 3, dpi = 300)

### 2D UMAP of NMPs, with cells are colored by their timepoints (Fig. 2c)
p = pd_x %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.35) +
    theme_void() +
    scale_color_manual(values=somite_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "NMPs_somite_count.png"), width = 4, height = 3, dpi = 300)

###############################
### Performing PCA analysis ###
###############################

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count_x = doExtractData(pd_x, mouse_gene_sub)
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

emb = res[["cell.embeddings"]]
emb = emb[rownames(pd_x),]
pd_x = cbind(pd_x, emb[,c(1:3)])
print(res[["varExplained"]])
pd_x$somite_count = factor(pd_x$somite_count, levels = names(somite_color_plate))

### making 3D PCA plot (Fig. 2e)
fig = plot_ly(pd_x, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~celltype_update, colors = posterior_embryo_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (21.0%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (13.9%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (11.1%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(work_path, example_i, "_NMP_Mesoderm_PCA_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

#################################################################################
### making scatter plot between each PCs and gene expression or somite counts ###
#################################################################################

gene_count_x = gene_count[,rownames(pd_x)]
emb = emb[rownames(pd_x),]

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[c("ENSMUSG00000074637","ENSMUSG00000030699","ENSMUSG00000020160","ENSMUSG00000062327","ENSMUSG00000009900","ENSMUSG00000024987"),]
rownames(gene_count_x) = c("Sox2","Tbx6","Mesi1","T","Wnt3a","Cyp26a1")
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(as.vector(gene_count_x[1,]), as.vector(gene_count_x[2,])),
                gene = c(rep("Sox2", nrow(emb)), rep("Tbx6", nrow(emb))),
                PC_1 = c(as.vector(emb[,1]), as.vector(emb[,1])))

### Fig. 2f
df %>%
    ggplot(aes(PC_1, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 


df = data.frame(exp = c(as.vector(gene_count_x[3,]), as.vector(gene_count_x[4,]), as.vector(gene_count_x[5,]), as.vector(gene_count_x[6,])),
                gene = c(rep("Mesi1", nrow(emb)), rep("T", nrow(emb)), rep("Wnt3a", nrow(emb)), rep("Cyp26a1", nrow(emb))),
                PC_3 = c(as.vector(emb[,3]), as.vector(emb[,3]), as.vector(emb[,3]), as.vector(emb[,3])))

### Fig. 2f
df %>%
    ggplot(aes(PC_3, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set2")

df  = data.frame(somite_count = pd_x$somite_count,
                 PC_2 = as.vector(emb[,2]))
df$somite_count = factor(df$somite_count, levels = names(somite_color_plate))
df$somite = as.vector(gsub(" somites", "", df$somite_count))
df$somite = factor(df$somite, levels = c(0, 2:12, 14:18, 20:34))

### Fig. 2f
df %>%
    ggplot( aes(somite, PC_2, fill = somite_count)) + 
    geom_boxplot(outlier.shape = NA) + 
    labs(x="", y="", title="") +
    theme_classic(base_size = 5) +
    scale_fill_manual(values=somite_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black")) +
    NoLegend()


