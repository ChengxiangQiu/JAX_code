
#####################################
### Section - 2, Posterior embryo ###
#####################################

####################################
### Making 3D UMAP visualization ###
####################################

source("JAX_help_code.R")
source("JAX_color_code.R")

example_i = "posterior_embryo"

pd = read.csv(paste0(example_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
rownames(pd) = as.vector(pd$cell_id)
pd$somite_count = factor(pd$somite_count, levels = names(somite_color_plate))
pd$leiden_res_1 = paste0("cluster_", pd$leiden_res_1)

### temporal clustering three cell clusters 
cluster_tmp = rep("NMP_Mesoderm", nrow(pd))
cluster_tmp[pd$leiden_res_1 %in% c("cluster_14")] = "Notochord"
cluster_tmp[pd$leiden_res_1 %in% paste0("cluster_", c(8,9,11,18))] = "Gut"
pd$cluster_tmp = as.vector(cluster_tmp)

write.csv(paste0(example_i, "_adata_scale.obs.csv"))

### making 3D UMAP, with cells are colored by their initial cell type annotations (Fig. 2a)
fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update, colors = posterior_embryo_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(example_i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

### making 3D UMAP, with cells are colored by somite counts (Fig. 2b)
fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~somite_count, colors = somite_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(example_i, "_somite_count.html"), selfcontained = FALSE, libdir = "tmp")




