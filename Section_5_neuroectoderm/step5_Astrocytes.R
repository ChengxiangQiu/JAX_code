
##################################
### Section - 5, Neuroectoderm ###
##################################

##############################################
### Analyzing astrocytes from stages < E13 ###
##############################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Astrocytes"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. S19a

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.3) +
    theme_void() +
    scale_color_manual(values=astrocytes_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)



##############################################################
### Compositions changing over time for different subtypes ###
##############################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Astrocytes"; print(example_i)

pd = readRDS(paste0(work_path, "df_cell.rds"))
pd_sub = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

day_list = names(day_color_plate)
pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])
pd_1 = pd[pd$cell_id %in% pd_sub$cell_id[pd_sub$celltype_sub_clustering == "VA1 astrocytes"],]
pd_2 = pd[pd$cell_id %in% pd_sub$cell_id[pd_sub$celltype_sub_clustering == "VA2 astrocytes"],]
pd_3 = pd[pd$cell_id %in% pd_sub$cell_id[pd_sub$celltype_sub_clustering == "VA3 astrocytes"],]
pd_4 = pd[pd$cell_id %in% pd_sub$cell_id[pd_sub$celltype_sub_clustering == "Anterior astrocytes"],]

x1 = pd_1 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x2 = pd_2 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x3 = pd_3 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x4 = pd_4 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x = x2 %>% select(day, frac) %>% rename(VA2 = frac) %>% 
    left_join(x1 %>% select(day, frac) %>% rename(VA1 = frac), by = "day") %>% 
    left_join(x3 %>% select(day, frac) %>% rename(VA3 = frac), by = "day") %>% 
    left_join(x4 %>% select(day, frac) %>% rename(AA = frac), by = "day")
x$VA1[is.na(x$VA1)] = 0
x$VA2[is.na(x$VA2)] = 0
x$VA3[is.na(x$VA3)] = 0
x$AA[is.na(x$AA)] = 0
x = data.frame(day = rep(x$day, 4),
               frac = c(x$VA1, x$VA2, x$VA3, x$AA),
               major_trajectory = rep(c("VA1","VA2","VA3","AA"), each = nrow(x)))

x$major_trajectory = factor(x$major_trajectory, levels = c("VA1","VA2","VA3","AA"))

### Fig. S19b

p = x %>% 
    ggplot(aes(x=day, y=frac, fill = day)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(major_trajectory)) + 
    labs(x='',y='% of cells') +
    scale_fill_manual(values=day_color_plate_2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))



###############################################################################
### Mapping different subtypes of astrocytes to their potential progenitors ###
###############################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Neuroectoderm_derivative"
name = "Astrocytes"

pd_sub = readRDS(paste0(work_path, name, "_adata_scale.obs.rds"))

### this result was calculated by step4_Mapping_neuroectoderm_derivatives.R
dat = readRDS(paste0(work_path, example_i, ".MNN_pairs.rds"))

pd_back = readRDS(paste0(work_path, "Neuroectoderm_backbone_adata_scale.obs.rds"))
rownames(pd_back) = as.vector(pd_back$cell_id)

celltype_sub_clustering_list = paste0("VA", c(1:3), " astrocytes")

### Fig. S19d

for(i in celltype_sub_clustering_list){
    print(i)
    pd_sub_i = pd_sub %>% filter(celltype_sub_clustering == i) %>% pull(cell_id)
    dat_sub = dat %>% filter(A %in% pd_sub_i) %>% group_by(B) %>% tally() %>% rename(cell_id = B, freq = n)
    df = pd_back %>% select(UMAP_1 = UMAP_2d_1, UMAP_2 = UMAP_2d_2, cell_id, day) %>% left_join(dat_sub, by = "cell_id") 
    df$freq[is.na(df$freq)] = 0
    
    name_i = gsub("/", "_", i)
    name_i = gsub(" ", "_", name_i)
    
    try(ggplot() +
            geom_point(data = df[sample(1:nrow(df),100000),], aes(x = UMAP_1, y = UMAP_2), size=0.5, color = "grey80") +
            geom_point(data = df[df$freq != 0,], aes(x = UMAP_1, y = UMAP_2, color = freq), size=0.5) +
            theme_void() +
            scale_color_viridis() +
            theme(legend.position="none") + 
            ggsave(paste0(work_path, name_i, ".png"), width = 8, height = 6, dpi = 300), silent = T)
    
}



