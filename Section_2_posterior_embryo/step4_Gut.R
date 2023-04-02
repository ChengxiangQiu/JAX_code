
#####################################
### Section - 2, Posterior embryo ###
#####################################

#######################
### Analysis on Gut ###
#######################

source("JAX_help_code.R")
source("JAX_color_code.R")

example_i = "posterior_embryo"

i = "Gut"

pd_x = read.csv(paste0(example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)

p = pd_x %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.8) +
    theme_void() +
    scale_color_manual(values=celltype_color_plate) +
    theme(legend.position="none") + 
    ggsave("Gut_celltype.png", width = 4, height = 3, dpi = 300)

p = pd_x %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.8) +
    theme_void() +
    scale_color_manual(values=somite_color_plate) +
    theme(legend.position="none") + 
    ggsave("Gut_day.png", width = 4, height = 3, dpi = 300)


###############################
### Performing PCA analysis ###
###############################

gene_count = readRDS("posterior_embryo_gene_count.rds")
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

emb = res[["cell.embeddings"]]
emb = emb[rownames(pd_x),]
pd_x = cbind(pd_x, emb[,c(1:3)])
print(res[["varExplained"]])
pd_x$somite_count = factor(pd_x$somite_count, levels = names(somite_color_plate))

fig = plot_ly(pd_x, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~somite_count, colors = somite_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (28.7%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (11.4%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (6.7%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(example_i, "_Gut_PCA_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

