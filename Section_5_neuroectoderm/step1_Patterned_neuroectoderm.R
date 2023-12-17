
##################################
### Section - 5, Neuroectoderm ###
##################################

###############################################################
### Making 2D UMAP visualization of patterned neuroectoderm ###
###############################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Neuroectoderm_backbone"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. 4a

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=neuroectoderm_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)


#############################################################################
### Making 3D UMAP visualization of patterned neuroectoderm + derivatives ###
#############################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Neuroectoderm_derivative"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. 4b

fig = plot_ly(pd[sample(1:nrow(pd), 250000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~major_trajectory, colors = major_trajectory_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(work_path, example_i, "_major_trajectory.html"), selfcontained = FALSE, libdir = "tmp")

### Fig. 4c

fig = plot_ly(pd[sample(1:nrow(pd), 250000),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day, colors = day_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(work_path, example_i, "_day.html"), selfcontained = FALSE, libdir = "tmp")

### BACKUP ###

celltype_list = c("Telencephalon",
                  "Dorsal telencephalon",
                  "Hypothalamus",
                  "Diencephalon",
                  "Midbrain",
                  "Hypothalamus (Sim1+)",
                  "Anterior floor plate",
                  "Midbrain-hindbrain boundary",
                  "Anterior roof plate",
                  "Hindbrain",
                  "Floorplate and p3 domain",
                  "Spinal cord/r7/r8",
                  "Posterior roof plate")

day_list = c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
             "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
             "E12.25", "E12.5", "E12.75")

emb = as.matrix(pd[,c("UMAP_1","UMAP_2","UMAP_3")])

dist_1 = list()
for(day_i in day_list){
    print(day_i)
    emb_x = emb[pd$day == day_i & pd$celltype_update %in% celltype_list,]
    if(nrow(emb_x) > 10000){
        emb_x = emb_x[sample(1:nrow(emb_x), 10000),]
    }
    dist_1[[day_i]] = c(rdist(emb_x))
    
}

dist_2 = list()
for(day_i in day_list){
    print(day_i)
    emb_x = emb[pd$day == day_i & !pd$celltype_update %in% celltype_list,]
    if(nrow(emb_x) > 10000){
        emb_x = emb_x[sample(1:nrow(emb_x), 10000),]
    }
    dist_2[[day_i]] = c(rdist(emb_x))
}

df = NULL
for(day_i in day_list){
    df = rbind(df,
               data.frame(day = day_i,
                          dist = mean(dist_1[[day_i]]),
                          group = "patterned_neuroectoderm", stringsAsFactors = FALSE))
    df = rbind(df,
               data.frame(day = day_i,
                          dist = mean(dist_2[[day_i]]),
                          group = "derived_cell_types", stringsAsFactors = FALSE))
}
df$day = factor(df$day, levels = day_list)

df$day = factor(df$day, levels = rev(day_list))
p = df %>%
    ggplot(aes(x=day, y=dist, color=group, group=group)) +
    geom_line() +
    geom_point() +
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    coord_flip()





























