
########################
### Section - 4, Eye ###
########################

########################################################
### Making 3D and 2D UMAP visualization of Eye cells ###
########################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Eye"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

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

day_group_color_plate=rev(brewer.pal(11,"Spectral"))
day_group_color_plate=colorRampPalette(day_group_color_plate)(8)
names(day_group_color_plate) = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                 "E16-E17", "E17-E18", "E18-P0", "P0")

### Fig. 4a

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_update, colors = eye_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(work_path, example_i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")

x_table = table(pd$day_group)
pd_1 = pd[pd$day_group %in% names(x_table)[x_table <= 5000],]
pd_2 = pd[pd$day_group %in% names(x_table)[x_table > 5000],]
pd_2_x = pd_2 %>% group_by(day_group) %>% sample_n(5000)
pd_sub = rbind(pd_1, pd_2[pd_2$cell_id %in% pd_2_x$cell_id,])

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~day_group, colors = day_group_color_plate) %>% 
    layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                        yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                        zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2),
                        camera = list(eye = list(x = -0.8, y = 2, z = 1.5))))
saveWidget(fig, paste0(work_path, example_i, "_day_group.html"), selfcontained = FALSE, libdir = "tmp")


### Fig. S10a

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    theme(legend.position="none") +
    ggsave(paste0(work_path, example_i, ".2D_UMAP.celltype_update.png"), width = 6, height = 6, dpi = 300)

p = pd_sub[sample(1:nrow(pd_sub)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.3) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.day_group.png"), width = 6, height = 6, dpi = 300)



########################################################
### Estimating cell number from individual cell type ###
########################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
x = as.vector(pd_all$day)
x[pd_all$day == "E8.0-E8.5"] = "E8.5"
pd_all$day = as.vector(x)

cell_num = readRDS(paste0(work_path, "cell_num_prediction.rds"))

x = pd %>% group_by(day, celltype_update) %>% tally() %>% 
    left_join(pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    left_join(cell_num %>% select(day, cell_num_pred), by = "day") %>%
    mutate(estimated_num_log2 = log2(ceiling(cell_num_pred*n/total_n))) %>%
    mutate(normalized_num_log2 = log2(ceiling(100000*n/total_n))) 

day_list = names(day_color_plate)
x$day = factor(x$day, levels = day_list[day_list %in% x$day])
x$celltype_update = factor(x$celltype_update, levels = names(eye_color_plate))

x$day_value = as.vector(x$day)
x$day_value[x$day == "P0"] = "E19"
x$day_value = as.numeric(gsub("E","",x$day_value))

x_first_day = x %>% group_by(celltype_update) %>% filter(n >= 10) %>% 
    slice_min(order_by = day_value, n = 1) %>% mutate(first_day_value = day_value)

x_sub = x %>% left_join(x_first_day %>% select(celltype_update, first_day_value), by = "celltype_update") %>%
    filter(day_value >= first_day_value)

### Fig. 4b

p = x %>%
    ggplot(aes(day_value, normalized_num_log2, color = celltype_update)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 1), axis.text.y = element_text(color="black")) +
    scale_x_continuous(breaks=seq(8.5, 19, 0.5)) +
    scale_color_manual(values=eye_color_plate) +
    NoLegend()

### Fig. S10e

p = x_sub %>% 
    ggplot(aes(x=day, y=estimated_num_log2, fill = celltype_update)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(celltype_update)) + 
    labs(x='',y='Log2(Estimated # of cells)') +
    scale_fill_manual(values=eye_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))



#############################
### Subclustering of RGCs ###
#############################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "Eye_RGC_adata_scale.obs.rds"))

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
day_group_color_plate=rev(brewer.pal(11,"Spectral"))
day_group_color_plate=colorRampPalette(day_group_color_plate)(8)
names(day_group_color_plate) = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                 "E16-E17", "E17-E18", "E18-P0", "P0")

### Fig. 4d

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.45) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "RGCs_day.png"), width = 5, height = 5, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = leiden_cluster), size=0.45) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "RGCs_leiden_cluster.png"), width = 5, height = 5, dpi = 300)


dat = readRDS(paste0(work_path, "Eye_RGC_heatmap_dat.rds"))

### Fig. 4e

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
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



#################################
### Focusing on days <= E12.5 ###
#################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "Eye_early_adata_scale.obs.rds"))

### Fig. S10c

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.35) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Eye_early_celltype.png"), width = 5, height = 5, dpi = 300)

pd$day_x = as.numeric(as.vector(gsub("E","", as.vector(pd$day))))

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_x), size=0.35) +
    theme_void() +
    scale_color_viridis(option = "B", breaks = c(8.5, 9.5, 10.5, 11.5, 12.5), labels = c("E8.5","E9.5","E10.5","E11.5","E12.5")) + 
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Eye_early_day.png"), width = 5, height = 5, dpi = 300)



########################
### Focusing on iris ###
########################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "Eye_iris_adata_scale.obs.rds"))

### Fig. S10f

p = pd[sample(1:nrow(pd)),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.35) +
    theme_void() +
    scale_color_manual(values=eye_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Eye_iris_celltype.png"), width = 5, height = 5, dpi = 300)

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
day_group_color_plate=rev(brewer.pal(11,"Spectral"))
day_group_color_plate=colorRampPalette(day_group_color_plate)(8)
names(day_group_color_plate) = c("Before E13", "E13-E14", "E14-E15", "E15-E16", 
                                 "E16-E17", "E17-E18", "E18-P0", "P0")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day_group), size=0.45) +
    theme_void() +
    scale_color_manual(values=day_group_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Eye_iris_day.png"), width = 5, height = 5, dpi = 300)

