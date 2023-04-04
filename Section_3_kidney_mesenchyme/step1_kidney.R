
######################################
### Section - 3, Kidney_mesenchyme ###
######################################

####################################
### Making 2D UMAP visualization ###
####################################

source("JAX_help_code.R")
source("JAX_color_code.R")

example_i = "Renal"; print(example_i)

pd = readRDS(paste0(example_i, "_adata_scale.obs.rds"))

### Fig. 3a
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)


### Fig. 3b
### Here, we pooled cells before E10
day_group = rep("Before E10", nrow(pd))

day_group[pd$day %in% c("E10.0", "E10.25", "E10.5", "E10.75")] = "E10-E11"
day_group[pd$day %in% c("E11.0", "E11.25", "E11.5", "E11.75")] = "E11-E12"
day_group[pd$day %in% c("E12.0", "E12.25", "E12.5", "E12.75")] = "E12-E13"
day_group[pd$day %in% c("E13.0", "E13.25", "E13.5", "E13.75")] = "E13-E14"
day_group[pd$day %in% c("E14.0", "E14.25", "E14.333", "E14.75")] = "E14-E15"
day_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"
day_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"
day_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"
day_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"
day_group[pd$day %in% c("P0")] = "P0"

day_group_list = c("Before E10", "E10-E11", "E11-E12", "E12-E13", "E13-E14", "E14-E15", "E15-E16", 
                   "E16-E17", "E17-E18", "E18-P0", "P0")
pd$day_group = factor(day_group, levels = day_group_list)

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_group_list))
names(day_color_plate_2) = day_group_list

day_group_table = table(pd$day_group)
pd_1 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table > 5000]) %>% group_by(day_group) %>% sample_n(5000) %>% as.data.frame()
pd_2 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table <= 5000]) %>% as.data.frame()
pd_sub = rbind(pd_1, pd_2)

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.5) +
    scale_color_manual(values=day_color_plate_2) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(example_i, ".day_group.2D_UMAP.png"), width = 6, height = 6, dpi = 300)



#########################################################
### Making 2D UMAP visualization on a subset of cells ###
#########################################################

### Here, we re-embedded cells from 
### a) Collecting duct intercalated cells
### b) Connecting tubule
### c) Collecting duct principal cells

source("JAX_help_code.R")
source("JAX_color_code.R")

example_i = "Renal_CDI"

pd = readRDS(paste0(example_i, "_adata_scale.obs.rds"))

### Fig. 3d
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)


day_group = rep("Before E15", nrow(pd))

day_group[pd$day %in% c("E15.0", "E15.25", "E15.5", "E15.75")] = "E15-E16"
day_group[pd$day %in% c("E16.0", "E16.25", "E16.5", "E16.75")] = "E16-E17"
day_group[pd$day %in% c("E17.0", "E17.25", "E17.5", "E17.75")] = "E17-E18"
day_group[pd$day %in% c("E18.0", "E18.25", "E18.5", "E18.75")] = "E18-P0"
day_group[pd$day %in% c("P0")] = "P0"

day_group_list = c("Before E15", "E15-E16", 
                   "E16-E17", "E17-E18", "E18-P0", "P0")
pd$day_group = factor(day_group, levels = day_group_list)

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_group_list))
names(day_color_plate_2) = day_group_list

day_group_table = table(pd$day_group)
pd_1 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table > 1000]) %>% group_by(day_group) %>% sample_n(500) %>% as.data.frame()
pd_2 = pd %>% filter(day_group %in% names(day_group_table)[day_group_table <= 1000]) %>% as.data.frame()
pd_sub = rbind(pd_1, pd_2)

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.5) +
    scale_color_manual(values=day_color_plate_2) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(example_i, ".day_group.2D_UMAP.png"), width = 6, height = 6, dpi = 300)



