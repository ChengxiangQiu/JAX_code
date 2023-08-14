
######################################
### Section - 3, Kidney_mesenchyme ###
######################################

####################################################
### Making 2D UMAP visualization of Kidney cells ###
####################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Renal"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. 3a
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)


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
    ggsave(paste0(work_path, example_i, ".day_group.2D_UMAP.png"), width = 6, height = 6, dpi = 300)



#########################################################
### Making 2D UMAP visualization on a subset of cells ###
#########################################################

### Here, we re-embedded cells from 
### a) Collecting duct intercalated cells
### b) Connecting tubule
### c) Collecting duct principal cells

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Renal_CDI"

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. 3d
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)


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
    ggsave(paste0(work_path, example_i, ".day_group.2D_UMAP.png"), width = 6, height = 6, dpi = 300)


#######################################################
### Ploting cells for before E18.75, E18.75, and P0 ###
#######################################################

### Fig. S10d

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Renal"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

day_group = rep("Early", nrow(pd))
day_group[pd$day == "E18.75"] = "E18.75"
day_group[pd$day == "P0"] = "P0"
pd$day_group = as.vector(day_group)

day_group_color_plate = c("Early" = "#a46cb7",
                          "E18.75" = "#7aa457",
                          "P0" = "#cb6a49",
                          "Other" = "grey90")

pd$tmp = if_else(pd$day_group == "Early", "Early", "Other")
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
    geom_point(data = subset(pd, tmp == 'Early'),
               aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, "_1.png"), width = 6, height = 6, dpi = 300)

pd$tmp = if_else(pd$day_group == "E18.75", "E18.75", "Other")
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
    geom_point(data = subset(pd, tmp == 'E18.75'),
               aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, "_2.png"), width = 6, height = 6, dpi = 300)

pd$tmp = if_else(pd$day_group == "P0", "P0", "Other")
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
    geom_point(data = subset(pd, tmp == 'P0'),
               aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, "_3.png"), width = 6, height = 6, dpi = 300)


########################################################
### Estimating cell number from individual timepoint ###
########################################################

### normalization by the total number of cells from each stage

example_i = "Renal"; print(example_i)

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))
cell_num = readRDS(paste0(work_path, "cell_num_prediction.rds"))

x = as.vector(pd_all$day)
x[pd_all$day == "E8.0-E8.5"] = "E8.5"
pd_all$day = as.vector(x)

day_list = names(day_color_plate)
pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])
pd_all$day = factor(pd_all$day, levels = day_list[day_list %in% pd_all$day])

x = pd %>% group_by(day, celltype_update) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    left_join(cell_num %>% select(day, cell_num_pred), by = "day") %>%
    mutate(estimated_num_log2 = log2(ceiling(cell_num_pred*n/total_n))) %>%
    mutate(normalized_num_log2 = log2(ceiling(100000*n/total_n))) 

x$day = factor(x$day, levels = day_list[day_list %in% x$day])
x$celltype_update = factor(x$celltype_update, levels = names(renal_color_plate))

x$day_value = as.vector(x$day)
x$day_value[x$day == "P0"] = "E19"
x$day_value = as.numeric(gsub("E","",x$day_value))

x_first_day = x %>% group_by(celltype_update) %>% filter(n >= 10) %>% 
    slice_min(order_by = day_value, n = 1) %>% mutate(first_day_value = day_value)

x_sub = x %>% left_join(x_first_day %>% select(celltype_update, first_day_value), by = "celltype_update") %>%
    filter(day_value >= first_day_value)

### Fig. S10b
p = x_sub %>% 
    ggplot(aes(x=day, y=estimated_num_log2, fill = celltype_update)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(celltype_update)) + 
    labs(x='',y='Log2(Estimated # of cells)') +
    scale_fill_manual(values=renal_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))


