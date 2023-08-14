
######################################
### Section - 3, Kidney_mesenchyme ###
######################################

###################################################################################
### Checking two subpopulations which are spatially mapping to kidney from the 
### lateral plate & intermediate mesoderm
###################################################################################

### Supplementary Figure 13

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell.rds"))

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

celltype_include = c("Posterior intermediate mesoderm",
                     "Anterior intermediate mesoderm",
                     "Metanephric mesenchyme",
                     "Renal stromal cells",
                     "Renal pericytes and mesangial cells",
                     "Splanchnic mesoderm")

pd_sub = pd_all[pd_all$celltype_update %in% celltype_include | 
                    (pd_all$lateral_plate_mesoderm_sub_clustering %in% celltype_include & !is.na(pd_all$lateral_plate_mesoderm_sub_clustering)),]
pd_x_1 = pd_sub[is.na(pd_sub$lateral_plate_mesoderm_sub_clustering),]
pd_x_2 = pd_sub[!is.na(pd_sub$lateral_plate_mesoderm_sub_clustering),]
pd_x_2$celltype_update = as.vector(pd_x_2$lateral_plate_mesoderm_sub_clustering)
pd_sub = rbind(pd_x_1, pd_x_2)
rownames(pd_sub) = as.vector(pd_sub$cell_id)
### n = 206,908 cells

gene_count = doExtractData(pd_sub, mouse_gene_sub)

name = "Renal_big"
Matrix::writeMM(t(gene_count), paste0(work_path, name, ".gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, name, ".df_cell.csv"))
write.csv(mouse_gene_sub, paste0(work_path, name, ".df_gene.csv"))

### python Embedding_Renal.py

df = read.csv(paste0(work_path, name, ".obs.csv"), header=T, row.names=1, as.is=T)
day_list = names(day_color_plate)
df$day = factor(df$day, levels = day_list[day_list %in% df$day])

p = df %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, name, "_celltype.png"), width = 6, height = 6, dpi = 300)

pd = df
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
    ggsave(paste0(work_path, name, "_day.png"), width = 6, height = 6, dpi = 300)



#########################################################
### Only co-embedding renal pericytes and stromal cells


source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell.rds"))

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

celltype_include = c("Renal stromal cells",
                     "Renal pericytes and mesangial cells")

pd_sub = pd_all[pd_all$celltype_update %in% celltype_include | 
                    (pd_all$lateral_plate_mesoderm_sub_clustering %in% celltype_include & !is.na(pd_all$lateral_plate_mesoderm_sub_clustering)),]
pd_x_1 = pd_sub[is.na(pd_sub$lateral_plate_mesoderm_sub_clustering),]
pd_x_2 = pd_sub[!is.na(pd_sub$lateral_plate_mesoderm_sub_clustering),]
pd_x_2$celltype_update = as.vector(pd_x_2$lateral_plate_mesoderm_sub_clustering)
pd_sub = rbind(pd_x_1, pd_x_2)
rownames(pd_sub) = as.vector(pd_sub$cell_id)
### n = 39,468 cells

gene_count = doExtractData(pd_sub, mouse_gene_sub)

name = "Renal_pericytes_stromal"
Matrix::writeMM(t(gene_count), paste0(work_path, name, ".gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, name, ".df_cell.csv"))
write.csv(mouse_gene_sub, paste0(work_path, name, ".df_gene.csv"))

cds <- new_cell_data_set(gene_count,
                         cell_metadata = pd_sub,
                         gene_metadata = mouse_gene_sub)
saveRDS(cds, paste0(work_path, "Renal_pericytes_stromal_cds.rds"))

### python Embedding_Renal.py

df = read.csv(paste0(work_path, name, ".obs.csv"), header=T, row.names=1, as.is=T)
day_list = names(day_color_plate)
df$day = factor(df$day, levels = day_list[day_list %in% df$day])

p = df %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_manual(values=renal_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, name, "_celltype.png"), width = 6, height = 6, dpi = 300)

pd = df
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
    ggsave(paste0(work_path, name, "_day.png"), width = 6, height = 6, dpi = 300)



#################################################
### Plotting Foxd1 expression across timepoints

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

cds = readRDS(paste0(work_path, "Renal_pericytes_stromal_cds.rds"))
exp = exprs(cds)
exp_norm = t(t(exp)/as.vector(cds$Size_Factor))

df = data.frame(pData(cds))
df$Foxd1_exp = as.vector(exp["ENSMUSG00000078302",])
df$Foxd1_exp_norm = as.vector(exp_norm["ENSMUSG00000078302",])
day_list = names(day_color_plate)
df$day = factor(df$day, levels = rev(day_list[day_list %in% df$day]))

### mean expression of Foxd1
p1 = df %>% filter(lateral_plate_mesoderm_sub_clustering == "Renal pericytes and mesangial cells") %>% 
    group_by(day) %>% summarise(mean_exp = mean(Foxd1_exp_norm)) %>%
    ggplot(aes(x=day, y=mean_exp, fill = day)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=day_color_plate) +
    theme(legend.position="none") +
    labs(x="", y="", title="") +
    coord_flip()

p2 = df %>% filter(lateral_plate_mesoderm_sub_clustering == "Renal stromal cells") %>% 
    group_by(day) %>% summarise(mean_exp = mean(Foxd1_exp_norm)) %>% complete(day, fill = list(mean_exp = 0)) %>%
    ggplot(aes(x=day, y=mean_exp, fill = day)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=day_color_plate) +
    theme(legend.position="none") +
    labs(x="", y="", title="") +
    coord_flip()

