
###################################
### Section - 1, Basic analysis ###
###################################

#########################################
### Pseudobulk analysis using Monocle ###
#########################################

source("JAX_help_code.R")
source("JAX_color_code.R")

cds = readRDS("embryo_cds.rds")

### identifying the highly variable genes
obj = doObjectTransform(cds, transform_to = "seurat")
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
gene_use = VariableFeatures(obj)

### performing PCA analysis
set.seed(2016)
FM = monocle3:::normalize_expr_data(cds, 
                                    norm_method = "log", 
                                    pseudo_count = 1)
FM = FM[gene_use,]

num_dim = 10
scaling = TRUE
set.seed(2016)
irlba_res = my_sparse_prcomp_irlba(Matrix::t(FM), 
                                   n = min(num_dim, min(dim(FM)) - 1), 
                                   center = scaling, 
                                   scale. = scaling)
preproc_res = irlba_res$x
row.names(preproc_res) = colnames(cds)

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)

df = data.frame(embryo_id = rownames(preproc_res),
                PC_1 = preproc_res[,1],
                PC_2 = preproc_res[,2],
                PC_3 = preproc_res[,3],
                day = as.vector(cds$day),
                embryo_sex = as.vector(cds$embryo_sex))
df$day = factor(df$day, levels = names(day_color_plate))

fig = plot_ly(df, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~day, colors = day_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (77.3%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (9.9%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (4.2%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))),
           showlegend = FALSE)
saveWidget(fig, "embryo_pca_day.html", selfcontained = FALSE, libdir = "tmp")

sex_color_plate = c("F" = "#ff0000",
                    "M" = "#0000FF")
fig = plot_ly(df, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~embryo_sex, colors = sex_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='PC_1 (77.3%)', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='PC_2 (9.9%)', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='PC_3 (4.2%)', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))),
           showlegend = FALSE)
saveWidget(fig, "embryo_pca_sex.html", selfcontained = FALSE, libdir = "tmp")













