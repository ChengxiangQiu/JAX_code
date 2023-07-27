
######################################
### Section - 3, Kidney_mesenchyme ###
######################################

##############################################################################
### Making 2D UMAP visualization for lateral plate & intermediate mesoderm ###
##############################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "LPM"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. 3f
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = lateral_plate_mesoderm_sub_clustering), size=0.2) +
    theme_void() +
    scale_color_manual(values=LPM_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 8, height = 6, dpi = 300)


x_table = table(pd$day)
pd_1 = pd %>% filter(day %in% names(x_table)[x_table > 10000]) %>% group_by(day) %>% sample_n(10000) %>% as.data.frame()
pd_2 = pd %>% filter(day %in% names(x_table)[x_table <= 10000]) %>% as.data.frame()
pd_sub = rbind(pd_1, pd_2)
pd_sub$day = factor(pd_sub$day, levels = names(day_color_plate))

### Fig. 3f (sub panel on the top left)
p = ggplot() +
    geom_point(data = pd_sub, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.35) +
    geom_point(data = pd_sub[sample(1:nrow(pd_sub)),], aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.2) +
    scale_color_manual(values=day_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".day.2D_UMAP.png"), width = 8, height = 6, dpi = 300)


