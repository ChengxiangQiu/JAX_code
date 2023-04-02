#3) UMAP of all inferred progenitors + first 500 cells (for dl1-5 only). 
#Color by dl assignment as well as by time in two different versions. Then color 
#by expression of each of the 18 TFs. I'd like to see how things cluster out and 
#also how these patterns look (esp. given how clustered the inferred progenitors are in Fig. 4j)

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

### progenitor cells only

dat = readRDS(paste0(work_path, "/Neurons/IN_backbone_prog.rds"))
dat_n = dat[["dat_n"]]
pd = dat[["pd_sub"]]
gene_count = dat[["gene_count"]]

day_value = as.vector(pd$somite_count)
day_value[pd$day == "E10.25"] = "35 somites"
day_value[pd$day == "E10.5"]  = "36 somites"
day_value[pd$day == "E10.75"] = "37 somites"
day_value[pd$day == "E11.0"]  = "38 somites"
day_value[pd$day == "E11.25"] = "39 somites"
day_value[pd$day == "E11.5"]  = "40 somites"
day_value[pd$day == "E11.75"] = "41 somites"
day_value[pd$day == "E12.0"]  = "42 somites"
day_value[pd$day == "E12.25"] = "43 somites"
day_value[pd$day == "E12.5"]  = "44 somites"
day_value[pd$day == "E12.75"] = "45 somites"
pd$somite_count = as.vector(day_value)

pd = pd[,c("cell_id","day","somite_count","celltype_update")]
pd$group = "prog"
pd$dI1_prog = if_else(pd$cell_id %in% as.vector(dat_n$cell_id[dat_n$celltype_derive == "Spinal dI1 interneurons"]), "Yes", "No")
pd$dI2_prog = if_else(pd$cell_id %in% as.vector(dat_n$cell_id[dat_n$celltype_derive == "Spinal dI2 interneurons"]), "Yes", "No")
pd$dI3_prog = if_else(pd$cell_id %in% as.vector(dat_n$cell_id[dat_n$celltype_derive == "Spinal dI3 interneurons"]), "Yes", "No")
pd$dI4_prog = if_else(pd$cell_id %in% as.vector(dat_n$cell_id[dat_n$celltype_derive == "Spinal dI4 interneurons"]), "Yes", "No")
pd$dI5_prog = if_else(pd$cell_id %in% as.vector(dat_n$cell_id[dat_n$celltype_derive == "Spinal dI5 interneurons"]), "Yes", "No")

Matrix::writeMM(t(gene_count), paste0(work_path, "/Neurons/IN_prog.gene_count.mtx"))
write.csv(pd, paste0(work_path, "/Neurons/IN_prog.df_cell.csv"))

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

Idents(obj) = as.vector(obj$dI1_prog)
res = FindMarkers(obj, ident.1 = "Yes", ident.2 = "No", only.pos = T)
res_1 = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>% filter(p_val_adj < 0.05)
res_1$group = "dI1"

Idents(obj) = as.vector(obj$dI2_prog)
res = FindMarkers(obj, ident.1 = "Yes", ident.2 = "No", only.pos = T)
res_2 = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>% filter(p_val_adj < 0.05)
res_2$group = "dI2"

Idents(obj) = as.vector(obj$dI3_prog)
res = FindMarkers(obj, ident.1 = "Yes", ident.2 = "No", only.pos = T)
res_3 = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>% filter(p_val_adj < 0.05)
res_3$group = "dI3"

Idents(obj) = as.vector(obj$dI4_prog)
res = FindMarkers(obj, ident.1 = "Yes", ident.2 = "No", only.pos = T)
res_4 = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>% filter(p_val_adj < 0.05)
res_4$group = "dI4"

Idents(obj) = as.vector(obj$dI5_prog)
res = FindMarkers(obj, ident.1 = "Yes", ident.2 = "No", only.pos = T)
res_5 = res %>% mutate(gene_ID = rownames(res)) %>% left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>% filter(p_val_adj < 0.05)
res_5$group = "dI5"

res = rbind(res_1, res_2, res_3, res_4, res_5)
write.csv(res, paste0(work_path, "/Neurons/IN_prog_DEGs.csv"))

res = read.csv("~/work/jax/rna_seq/Neurons/tmp/IN_prog/IN_prog_DEGs.csv", row.names=1, as.is=T)
df = res %>% mutate(log2_fc_pct = log2(pct.1/pct.2), log10_qval = -log10(p_val_adj))

p = df %>%
    ggplot() +
    geom_point(aes(x = log2_fc_pct, y = log10_qval), size=2.5) +
    geom_point(aes(x = log2_fc_pct, y = log10_qval, color = group), size=2.3) +
    theme_classic(base_size = 10) +
    #scale_color_manual(values=neuron_color_plate_sub) +
    scale_color_brewer(palette = "Set2") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    geom_text(data=subset(df, log2_fc_pct >= 4 | log10_qval >= 100),
              aes(log2_fc_pct,log10_qval,label=gene_short_name), hjust = 0, nudge_x = 0.1, color = "black", size = 2.5)

pdf("~/work/jax/rna_seq/Neurons/tmp/IN_prog/Vol_plot.pdf", 8, 5)
print(p)
dev.off()

df = df %>% group_by(group) %>% arrange(desc(log2_fc_pct), .by_group=TRUE)
write.csv(df[,c("gene_ID", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene_short_name", "group")], 
          "~/work/jax/rna_seq/Neurons/tmp/IN_prog/IN_prog_DEGs.output.csv")


### adding first 500 cells

obj = readRDS(paste0(work_path,"/Neurons/IN_early_cells.rds"))

gene_count_2 = GetAssayData(obj, slot = "counts")
pd_2 = data.frame(obj[[]])[,c("cell_id","day","somite_count","celltype_update")]
pd_2$group = "early"
pd_2$dI1_prog = NA
pd_2$dI2_prog = NA
pd_2$dI3_prog = NA
pd_2$dI4_prog = NA
pd_2$dI5_prog = NA

gene_count_merge = cbind(gene_count, gene_count_2)
pd_merge = rbind(pd, pd_2)

Matrix::writeMM(t(gene_count_merge), paste0(work_path, "/Neurons/IN_prog_early.gene_count.mtx"))
write.csv(pd_merge, paste0(work_path, "/Neurons/IN_prog_early.df_cell.csv"))

### checking the result ###

example_i = "IN_prog"

pd = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.csv"), row.names=1, as.is=T)

obj = CreateSeuratObject(gene_count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Neurons/", example_i, ".cds.rds"))

p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI1_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1.2) +
    geom_point(data = pd %>% filter(dI1_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#636EFA", size=1) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI1_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI2_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1.2) +
    geom_point(data = pd %>% filter(dI2_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#EF553B", size=1) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI2_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI3_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1.2) +
    geom_point(data = pd %>% filter(dI3_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#00CC96", size=1) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI3_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI4_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1.2) +
    geom_point(data = pd %>% filter(dI4_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#AB63FA", size=1) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI4_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI5_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1.2) +
    geom_point(data = pd %>% filter(dI5_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#FFA15A", size=1) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI5_prog.png"), width = 4, height = 4, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=1) +
    theme_void() +
    scale_color_manual(values=neuron_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_update.png"), width = 4, height = 4, dpi = 300)


day_value = as.vector(pd$somite_count)
day_value[pd$day == "E10.25"] = "E10.25"
day_value[pd$day == "E10.5"]  = "E10.5" 
day_value[pd$day == "E10.75"] = "E10.75"
day_value[pd$day == "E11.0"]  = "E11.0" 
day_value[pd$day == "E11.25"] = "E11.25"
day_value[pd$day == "E11.5"]  = "E11.5" 
day_value[pd$day == "E11.75"] = "E11.75"
day_value[pd$day == "E12.0"]  = "E12.0" 
day_value[pd$day == "E12.25"] = "E12.25"
day_value[pd$day == "E12.5"]  = "E12.5" 
day_value[pd$day == "E12.75"] = "E12.75"
pd$somite_count = as.vector(day_value)

day_list_sub = names(table(pd$somite_count))
library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_sub))
names(day_color_plate_2) = day_list_sub
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=1) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite_count.png"), width = 4, height = 4, dpi = 300)
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=1) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite_count.legend.png"), width = 6, height = 4, dpi = 300)



day_list_sub = names(table(pd$day))
day_list_sub = day_list[day_list %in% day_list_sub]
pd$day = factor(pd$day, levels = day_list_sub)
library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_sub))
names(day_color_plate_2) = day_list_sub
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=1) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".day.png"), width = 4, height = 4, dpi = 300)
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=1) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".day.legend.png"), width = 5, height = 4, dpi = 300)

example_i = "IN_prog"
cds = readRDS(paste0("~/work/jax/rna_seq/Neurons/", example_i, ".cds.rds"))
gene_list = c("Barhl1","Ebf3","E2f1","Pou3f2","Ebf2","Nr2f1","Neurog2","Nhlh1","Hoxc4","Hoxb6","Prdm13","Pax3","Npas3","Tfap2b","Casz1","Onecut2")
my_plot_cells(cds, genes = gene_list, how_many_rows = 3, cell_size = 0.6)  + 
    ggsave(paste0("./tmp/IN_prog/", example_i, ".gene_exp.png"), width = 12, height = 7, dpi = 300)

gene_list = c("Barhl1", "Neurog1", "Hoxa5", "Nphs1", "Grm4")
my_plot_cells(cds, genes = gene_list, how_many_rows = 1, cell_size = 0.8) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave(paste0("~/work/jax/rna_seq/Neurons/tmp/IN_prog/", example_i, ".gene_expression_9.png"), width = 10, height = 2, dpi = 300)






### After adding the first 500 cells

example_i = "IN_prog_early"

pd = read.csv(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.csv"), row.names=1, as.is=T)

obj = CreateSeuratObject(gene_count_merge, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = doObjectTransform(obj, transform_to = "monocle3")
cds = preprocess_cds(cds, use_genes = genes_include)
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/Neurons/", example_i, ".cds.rds"))

p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(group == "prog"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(group == "prog"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "blue", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".group_prog.png"), width = 4, height = 4, dpi = 300)

p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(group == "early"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(group == "early"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "red", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".group_early.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI1_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(dI1_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#636EFA", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI1_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI2_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(dI2_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#EF553B", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI2_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI3_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(dI3_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#00CC96", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI3_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI4_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(dI4_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#AB63FA", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI4_prog.png"), width = 4, height = 4, dpi = 300)


p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.6) +
    geom_point(data = pd %>% filter(dI5_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "black", size=1) +
    geom_point(data = pd %>% filter(dI5_prog == "Yes"), aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "#FFA15A", size=0.8) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".dI5_prog.png"), width = 4, height = 4, dpi = 300)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.8) +
    theme_void() +
    scale_color_manual(values=neuron_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".celltype_update.png"), width = 4, height = 4, dpi = 300)


day_value = as.vector(pd$somite_count)
day_value[pd$day == "E10.25"] = "E10.25"
day_value[pd$day == "E10.5"]  = "E10.5" 
day_value[pd$day == "E10.75"] = "E10.75"
day_value[pd$day == "E11.0"]  = "E11.0" 
day_value[pd$day == "E11.25"] = "E11.25"
day_value[pd$day == "E11.5"]  = "E11.5" 
day_value[pd$day == "E11.75"] = "E11.75"
day_value[pd$day == "E12.0"]  = "E12.0" 
day_value[pd$day == "E12.25"] = "E12.25"
day_value[pd$day == "E12.5"]  = "E12.5" 
day_value[pd$day == "E12.75"] = "E12.75"
pd$somite_count = as.vector(day_value)

day_list_sub = names(table(pd$somite_count))
library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_sub))
names(day_color_plate_2) = day_list_sub
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.8) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite_count.png"), width = 4, height = 4, dpi = 300)
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.8) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".somite_count.legend.png"), width = 6, height = 4, dpi = 300)



day_list_sub = names(table(pd$day))
day_list_sub = day_list[day_list %in% day_list_sub]
pd$day = factor(pd$day, levels = day_list_sub)
library(RColorBrewer)
day_color_plate_2 = rev(brewer.pal(11,"Spectral"))
day_color_plate_2 = colorRampPalette(day_color_plate_2)(length(day_list_sub))
names(day_color_plate_2) = day_list_sub
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.8) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".day.png"), width = 4, height = 4, dpi = 300)
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.8) +
    theme_void() +
    scale_color_manual(values=day_color_plate_2) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #theme(legend.position="none") + 
    ggsave(paste0("~/share/", example_i, ".day.legend.png"), width = 5, height = 4, dpi = 300)

example_i = "IN_prog_early"
cds = readRDS(paste0("~/work/jax/rna_seq/Neurons/", example_i, ".cds.rds"))
gene_list = c("Barhl1","Ebf3","E2f1","Pou3f2","Ebf2","Nr2f1","Neurog2","Nhlh1","Hoxc4","Hoxb6","Prdm13","Pax3","Npas3","Tfap2b","Casz1","Onecut2")
my_plot_cells(cds, genes = gene_list, how_many_rows = 3, cell_size = 0.6)  + 
    ggsave(paste0("./tmp/", example_i, ".gene_exp.png"), width = 12, height = 7, dpi = 300)

my_plot_cells(cds, genes = gene_list, how_many_rows = 4, cell_size = 0.5) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("gene_expression_9.png", dpi = 300, height = 2, width = 6)





#4) ONLY if it is easy, can you literally repeat the same analysis in dl1-5 
#but with all genes rather than only TFs? I am wondering if we can find other 
#genes that might be involved in early specification (we assume it's a TF code 
#but has anyone ever systematically checked?)

