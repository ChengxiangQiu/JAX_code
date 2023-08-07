

####################################
### Section - 7, Key TFs & genes ###
####################################

##################################################
### Anterior primitive streak -> Def. endoderm ###
##################################################

### Can we estimate pseudotime of cells corresponding to Def.endoderm developement 
### during transition and plot key TF/genes expression as a function of pseudotime?

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")

system_i = "Gastrulation"

edges = read.table(paste0(work_path, "edges.txt"), header=F, as.is=T, sep="\t")
names(edges) = c("system","x","y","x_name","y_name","edge_type")
edges = edges[edges$system == "Gastrulation",]

obj = readRDS(paste0(work_path, "obj_Early_PS.rds"))
pd = readRDS(paste0(work_path, system_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
pd_all = readRDS(paste0(work_path, "df_cell_graph.rds"))
pd$celltype_new = NULL
pd = pd %>% left_join(pd_all[,c("cell_id", "celltype_new")], by = "cell_id")

if ("meta_group" %in% names(pd)) {pd$meta_group = NULL}
pd_x = pd %>% left_join(nodes %>% filter(system == system_i) %>% select(celltype_new, meta_group))
pd$meta_group = as.vector(pd_x$meta_group)
gene_count = GetAssayData(obj, slot = "counts")

y = data.frame(i = 1:nrow(pd),
               j = 1:nrow(pd),
               meta_group = as.vector(pd$meta_group), stringsAsFactors = FALSE)

x = readRDS(paste0(work_path, system_i, ".MNN.rds"))
x_rev = x
x_rev$i = x$j
x_rev$j = x$i
x_rev = rbind(x, x_rev)
dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j")

cnt = 1
print(edges[cnt,])

xx = as.vector(edges$x)[cnt]
yy = as.vector(edges$y)[cnt]

dat_cnt = dat[dat$meta_group.x == xx & dat$meta_group.y == yy,]
group_1_MNN = as.vector(dat_cnt$i)
group_2_MNN = as.vector(dat_cnt$j)

coor = c(1:nrow(pd))

while(length(unique(group_1_MNN)) < 200){
    num_1 = length(unique(group_1_MNN))
    tmp = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
    group_1_MNN = c(group_1_MNN, tmp)
    num_2 = length(unique(group_1_MNN))
    if(num_1 == num_2) {break}
}

group_1_close = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
while(length(unique(group_1_close)) < 200){
    num_1 = length(unique(group_1_close))
    tmp = intersect(as.vector(dat$j[dat$i %in% group_1_close]), coor[pd$meta_group == xx])
    group_1_close = c(group_1_close, tmp)
    group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
    num_2 = length(unique(group_1_close))
    if(num_1 == num_2) {break}
}

while(length(unique(group_2_MNN)) < 200){
    num_1 = length(unique(group_2_MNN))
    tmp = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
    group_2_MNN = c(group_2_MNN, tmp)
    num_2 = length(unique(group_2_MNN))
    if(num_1 == num_2) {break}
}

group_2_close = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
while(length(unique(group_2_close)) < 200){
    num_1 = length(unique(group_2_close))
    tmp = intersect(as.vector(dat$j[dat$i %in% group_2_close]), coor[pd$meta_group == yy])
    group_2_close = c(group_2_close, tmp)
    group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
    num_2 = length(unique(group_2_close))
    if(num_1 == num_2) {break}
}

group = rep("other", nrow(pd))
group[coor %in% group_1_close] = "group_1"
group[coor %in% group_1_MNN] = "group_2"
group[coor %in% group_2_MNN] = "group_3"
group[coor %in% group_2_close] = "group_4"
pd$group = as.vector(group)

pd_sub = pd[pd$group != "other",]
group_table = table(pd_sub$group)

gene_count_sub = gene_count[,as.vector(pd_sub$cell_id)]
obj_sub = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)
obj_sub = FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj_sub)

cds = doObjectTransform(obj_sub, transform_to = "monocle")

pd_x = read.csv(paste0(work_path, "pijuan_obs.csv"), row.names=1, as.is=T)
pd_x = pd_x[colnames(cds),]

cds$batch_1 = as.vector(pd_x$batch)
cds$batch_2 = as.vector(pd_x$group)

cds = preprocess_cds(cds, use_genes = genes_include)
cds = align_cds(cds, alignment_group = "batch_1")
cds = reduce_dimension(cds)

saveRDS(cds, paste0(work_path, "cds_Def_endoderm.rds"))


#####################
### Making plots ####
#####################

plot_cells(cds, color_cells_by = "celltype_new", cell_size = 1)

cds = cluster_cells(cds)
cds = learn_graph(cds)
cds = order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1)

cds$pseudotime = cds@principal_graph_aux[["UMAP"]]$pseudotime

df = data.frame(pData(cds)) 
boxplot(df$pseudotime~factor(df$group))

df$UMAP_1 = reducedDims(cds)$UMAP[,1]
df$UMAP_2 = reducedDims(cds)$UMAP[,2]

### Fig. S16e

p = ggplot() +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=2, color = "black") +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2, color = celltype_new), size=1.8) +
    theme_void() +
    scale_color_manual(values=gastrulation_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Def_endoderm_UMAP_celltype.png"), width = 6, height = 4, dpi = 300)

p = ggplot() +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=2, color = "black") +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime), size=1.8) +
    theme_void() +
    scale_color_viridis(discrete=F) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Def_endoderm_UMAP_pseudotime.png"), width = 6, height = 4, dpi = 300)


###################################
### making gene expression plot ###
###################################

gene_count = exprs(cds)
gene_count = t(t(gene_count) / colSums(gene_count)) * 100000

target_genes = c("Sox17", "Elf3", "Sall4", "Hesx1", "Lin28a", "Ovol2",
                 "Cer1", "Slc25a4", "Cd24a",  "Slc2a3", "Lrpap1", "Krt18")

mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]

gene_count_x = gene_count[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

dat = data.frame(exp = c(t(as.matrix(gene_count_x))),
                 gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                 pseudotime = rep(as.vector(cds$pseudotime), nrow(gene_count_x)),
                 pseudotime_rank = rep(rank(as.vector(cds$pseudotime)), nrow(gene_count_x)), stringsAsFactors = F)
dat$gene = factor(dat$gene, levels = target_genes)

p = ggplot() +
    geom_smooth(data = dat, aes(pseudotime, exp, color = gene), method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Paired")


dat_add_offset = NULL
for(gene_i in target_genes){
    dat_sub = dat %>% filter(gene == gene_i)
    p_tmp = qplot(pseudotime, exp, data=dat_sub) + stat_smooth(method = loess, se = FALSE)
    dat_tmp = ggplot_build(p_tmp)$data[[2]][,c("x","y")]
    dat_tmp$y = dat_tmp$y - dat_tmp$y[dat_tmp$x == 0]
    dat_tmp$gene = gene_i
    
    dat_add_offset = rbind(dat_add_offset, dat_tmp)
}
dat_add_offset$gene = factor(dat_add_offset$gene, levels = target_genes)

p = ggplot() +
    geom_line(data = dat_add_offset, aes(x, y, color = gene), size = 1) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Paired")

### Fig. S16f

pdf(paste0(work_path, "Def_endoderm_gene_expression_pseudotime_add_offset.pdf"), 5, 3)
print(p)
dev.off()


