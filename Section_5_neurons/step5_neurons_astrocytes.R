### Focusing on Astrocytes from stages < E13


name = "Astrocytes"

pd_sub = read.csv(paste0(work_path, "/Neurons/", name, "_adata_scale.obs.csv"), header=T, row.names=1)

celltype_sub_clustering = rep("VA2 astrocytes", nrow(pd_sub))
celltype_sub_clustering[pd_sub$leiden_res_1 %in% c(4,6,7,11,12)] = "VA3 astrocytes"
celltype_sub_clustering[pd_sub$leiden_res_1 %in% c(1,9)] = "VA1 astrocytes"
celltype_sub_clustering[pd_sub$leiden_res_1 %in% c(14)] = "Anterior astrocytes"

pd_sub$celltype_sub_clustering = as.vector(celltype_sub_clustering)

obj = CreateSeuratObject(gene_count, meta.data = pd_sub)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Neurons/", name, ".cds.rds"))

### mapping to the backbone UMAP

dat = readRDS(paste0(work_path, "/Neurons/", example_i, ".MNN_pairs.rds"))

pd_back = readRDS(paste0(work_path, "/Neurons/Neurons_plus_big_backbone_adata_scale.obs.rds"))
rownames(pd_back) = as.vector(pd_back$cell_id)

celltype_sub_clustering_list = paste0("VA", c(1:3), " astrocytes")

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
            ggsave(paste0("~/share/", name_i, ".png"), width = 8, height = 6, dpi = 300), silent = T)
    
}

name = "Astrocytes"

astrocytes_color_plate = c("VA1 astrocytes" = "#9ac25e",
                           "VA2 astrocytes" = "#b35948",
                           "VA3 astrocytes" = "#b98b3a",
                           "Anterior astrocytes" = "#8d50a9")

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_sub_clustering), size=0.8) +
    theme_void() +
    scale_color_manual(values=astrocytes_color_plate) +
    theme(legend.position="none") + 
    #scale_color_brewer(palette = "Set2") +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    ggsave(paste0("~/share/", name, ".celltype_sub_clustering.png"), width = 6, height = 6, dpi = 300)


day_list_sub = day_list[day_list %in% pd_sub$day]
pd_sub$day = factor(pd_sub$day, levels = day_list_sub)

library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_sub))
names(day_color_plate_2) = day_list_sub

p = pd_sub %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.8) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    theme(legend.position="none") + 
    #scale_color_brewer(palette = "Set2") +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    ggsave(paste0("~/share/", name, ".celltype_day.png"), width = 6, height = 6, dpi = 300)


my_plot_cells(cds, genes = c("Otx2","Hoxb4","Hoxd4","Hoxc6"), how_many_rows = 1, cell_size = 0.5) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_7.png", dpi = 300, height = 2, width = 8)

my_plot_cells(cds, genes = c("Pax6","Reln","Nkx6-1","Slit1"), how_many_rows = 1, cell_size = 0.5) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_8.png", dpi = 300, height = 2, width = 8)




### making an easy plot to show the time series between different astrocytes

pd = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))

pd$day = factor(pd$day, levels = day_list)
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

p = x %>% 
    ggplot(aes(x=day, y=frac, fill = day)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(major_trajectory)) + 
    labs(x='',y='% of cells') +
    scale_fill_manual(values=day_color_plate_2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0("~/share/Neurons_density_2.pdf"), 4, 4)
p
dev.off()
