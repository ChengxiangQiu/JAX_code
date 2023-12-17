
#####################################
### Section - 6, Development tree ###
#####################################

##################################################################################################################
### Here, we processed data of cells from four previous studies, corresponding to pre- and gastrulation stages ###
##################################################################################################################

### E0 - E3.0; Xue et al.; GSE44183
### E3.5, E4.5, E5.5, E6.5; Mohammed et al.; GSE100597
### E5.25, E5.5, E6.25, E6.5; Cheng et al.; GSE109071
### E6.5 - E8.5; Pijuan-Sala et al.; E-MTAB-6967
### See more details in the Supplementary Table 17

##################################################################################
### First, combining cells from E3.5 to E8.5 which have been processed before. ###
##################################################################################

### In our previous project TOME (https://www.nature.com/articles/s41588-022-01018-x), we have reprocessed
### most of those published datasets, so we directly used those reprocessed profiles instead of performing
### extra analysis here.

### The data can be downloaded from https://tome.gs.washington.edu/
### Cells are organized into their initial timepoint, and then reprcessed & annotated.
### e.g., seurat_object_E6.5.rds

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

df_gene = read.csv(paste0(work_path, "df_gene_all.csv"), header=T, as.is=T, row.names=1)
obj_x = readRDS(paste0(work_path, "seurat_object_", "E6.75", ".rds"))
gene_overlap = as.vector(df_gene$gene_id[df_gene$gene_id %in% rownames(obj_x)])

time_point = c("E3.5", "E4.5", "E5.25", "E5.5", "E6.25", "E6.5", "E6.75", "E7.0", 
               "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5a")
obj = NULL

for(i in 1:length(time_point)){
    time_i = time_point[i]
    print(time_i)
    
    obj_i = readRDS(paste0(work_path, "seurat_object_", time_i, ".rds"))
    print(table(obj_i$group))
    
    if(i %in% c(1:5)){
        group_i = rep("Cheng", ncol(obj_i))
        group_i[obj_i$group == "351"] = "Mohammed"
        obj_i$group = as.vector(group_i)
    } else {
        group_i = rep("Pijuan_Sala", ncol(obj_i))
        group_i[obj_i$group == "351"] = "Mohammed"
        group_i[obj_i$group == "352"] = "Cheng"
        obj_i$group = as.vector(group_i)
    }
    print(table(obj_i$group))
    
    obj_i$cell_id = as.vector(colnames(obj_i))
    obj_i$sample = NULL
    obj_i$day = time_i
    
    obj_i = CreateSeuratObject(GetAssayData(obj_i, slot = "counts")[gene_overlap,], 
                               meta.data = data.frame(obj_i[[]]))
    
    if(is.null(obj)){
        obj = obj_i
    } else {
        obj = merge(obj, obj_i)
    }
}

saveRDS(obj, paste0(work_path, "obj_Early_PS.rds"))


##################################################################################################
### Second, integrating cells from P-S gastrulation dataset and E8-E8.5 subset of the new data ###
##################################################################################################

### Of note, due to the technology burden (10X vs. sci-), 
### we used seurat based integration method to correct batch effect
### For more details, please visit https://satijalab.org/seurat/articles/integration_rpca.html

library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

df_gene = read.csv(paste0(work_path, "df_gene_all.csv"), header=T, as.is=T, row.names=1)
obj_x = readRDS(paste0(work_path, "seurat_object_", "E6.75", ".rds"))
gene_overlap = as.vector(df_gene$gene_id[df_gene$gene_id %in% rownames(obj_x)])

obj_1 = readRDS(paste0(work_path, "obj_Early_PS.rds"))
obj_1 = subset(obj_1, subset = group == "Pijuan_Sala")

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)
x = as.vector(pd_all$day)
x[pd_all$day == "E8.0-E8.5"] = "E8.5"
pd_all$day = as.vector(x)

pd_all = pd_all[pd_all$day == "E8.5",]
pd_all$group = "Jax"
pd_all$cell_state = as.vector(pd_all$major_trajectory)
pd_all$cell_type = as.vector(pd_all$celltype_update)
pd_sub = pd_all[,c("day","group","cell_state","cell_type","cell_id")]

gene_count = doExtractData(pd_sub, mouse_gene[gene_overlap,])
obj = CreateSeuratObject(gene_count, pd)

obj_2 = CreateSeuratObject(gene_count, meta.data = pd_sub)

obj = merge(obj_1, obj_2)
print(table(obj$group))

obj.list <- SplitObject(obj, split.by = "group")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                  dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 2)

obj.integrated <- FindNeighbors(object = obj.integrated, dims = 1:30, reduction = "pca")
obj.integrated <- FindClusters(object = obj.integrated, resolution = 1)

saveRDS(obj.integrated, paste0(work_path, "obj_processed_PS_JaxE85.rds"))

pd = data.frame(obj.integrated[[]])
emb = Embeddings(obj.integrated, reduction = "umap")
pd$UMAP_1 = as.vector(emb[,1])
pd$UMAP_2 = as.vector(emb[,2])

saveRDS(pd, paste0(work_path, "pd_processed_PS_JaxE85.rds"))

emb = Embeddings(obj.integrated, reduction = "pca")
saveRDS(emb, paste0(work_path, "pd_processed_PS_JaxE85.PCs.rds"))

### Based on the co-embedding, we reannotated each individual cell clusters
### Updated cell type annotation is saved in "celltype_update" column in PS_JAXE8.5_integration.obs.rds

### Extended Data Fig. 11e

pd = readRDS(paste0(work_path, "PS_JAXE8.5_integration.obs.rds"))
ggplot() +
    geom_point(data = pd[sample(1:nrow(pd)),], aes(x = UMAP_1, y = UMAP_2), size=0.5, color = "black") +
    geom_point(data = pd[sample(1:nrow(pd)),], aes(x = UMAP_1, y = UMAP_2, color = celltype_update), size=0.35) +
    theme_void() +
    scale_color_manual(values=gastrulation_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "PS_JaxE85_celltype_update.png"), width = 8, height = 8, dpi = 300)

### Extended Data Fig. 11f

color_plate = c("Pijuan_Sala" = "#cc545e", "Jax" = "#9970c1")

for(i in names(color_plate)){
    try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "grey90") +
        geom_point(data = pd %>% filter(group == i), aes(x = UMAP_1, y = UMAP_2, color = group), size=0.5) +
        theme_void() +
        scale_color_manual(values=color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0(work_path, "pd_processed_PS_JaxE85.", i,".png"), width = 6, height = 6, dpi = 300))
}



###########################################################################################################################
### Third, making two heatmap, to present the consistent between the new annotation to previous P-S annotation or Jax's ###
###########################################################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "PS_JAXE8.5_integration.obs.rds"))

celltype_rank = c("Epiblast",
                  "Embryonic visceral endoderm",
                  "Extraembryonic visceral endoderm",
                  "Extraembryonic ectoderm",
                  "Parietal endoderm",
                  "Primitive streak",
                  "Nascent mesoderm",
                  "Extraembryonic mesoderm",
                  "Allantois",
                  "Amniotic mesoderm",
                  "Gut mesenchyme",
                  "Cardiopharyngeal mesoderm",
                  "Lateral plate mesoderm",
                  "Intermediate mesoderm",
                  "Paraxial mesoderm (Tbx6-)",
                  "Paraxial mesoderm (Tbx6+)",
                  "Cardiogenic mesoderm",
                  "CLE and NMPs",
                  "Definitive ectoderm",
                  "Floor plate",
                  "Spinal cord",
                  "Forebrain",
                  "Midbrain",
                  "Hindbrain",
                  "Neural crest",
                  "Surface ectoderm",
                  "Amniotic ectoderm",
                  "Primordial germ cells",
                  "Anterior primitive streak",
                  "Definitive endoderm",
                  "Gut",
                  "Notochord",
                  "Hematoendothelial progenitors",
                  "Endothelium",
                  "Blood progenitors",
                  "Primitive erythroid cells")

#########################
### comparing to P-S

dat = pd %>% filter(group == "Pijuan_Sala") %>% 
    group_by(celltype_update, pre_celltype) %>% tally()
exp = dat %>%
    dcast(pre_celltype ~ celltype_update)
rownames(exp) = as.vector(exp[,1])
exp = exp[,-1]
exp[is.na(exp)] = 0

dat$celltype_update = factor(dat$celltype_update, levels = celltype_rank)
dat_x = dat %>% group_by(pre_celltype) %>% slice_max(order_by = n, n = 1) %>% arrange(celltype_update)
exp = exp[as.vector(dat_x$pre_celltype), celltype_rank]

### Extended Data Fig. 11g

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "overlap_with_P-S.pdf"), 8, 8)
heatmap.2(as.matrix(exp), 
          col=Colors, 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(5,5))
dev.off()

############################
### comparing to JAX E8.5

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
pd_all = pd_all[pd_all$cell_id %in% pd$cell_id,]
pd_x = pd %>% filter(group != "Pijuan_Sala") %>% left_join(pd_all[,c("cell_id","celltype_new")])
pd_x$cell_type = as.vector(pd_x$celltype_new)

celltype_include = table(pd_x$cell_type)
celltype_include = c(names(celltype_include)[celltype_include >= 150], "Primordial germ cells")

dat = pd_x %>% filter(cell_type %in% celltype_include) %>% 
    group_by(celltype_update, cell_type) %>% tally() 

exp = dat %>%
    dcast(cell_type ~ celltype_update)
rownames(exp) = as.vector(exp[,1])
exp = exp[,-1]
exp[is.na(exp)] = 0
names_tmp = colnames(exp)
exp = cbind(exp, rep(0, nrow(exp)))
exp = cbind(exp, rep(0, nrow(exp)))
exp = cbind(exp, rep(0, nrow(exp)))
exp = cbind(exp, rep(0, nrow(exp)))
colnames(exp) = c(names_tmp, c("Embryonic visceral endoderm","Extraembryonic ectoderm","Definitive endoderm","Blood progenitors"))

dat$celltype_update = factor(dat$celltype_update, levels = celltype_rank)
dat_x = dat %>% group_by(cell_type) %>% slice_max(order_by = n, n = 1) %>% arrange(celltype_update)
exp = exp[as.vector(dat_x$cell_type), celltype_rank]

### Extended Data Fig. 11h

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "overlap_with_Jax.pdf"), 8, 8)
heatmap.2(as.matrix(exp), 
          col=Colors, 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(5,5))
dev.off()





####################################################################
### Making the graph from gastrulation stages based on MNN pairs ###
####################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "pd_processed_PS_JaxE85.rds"))
emb = readRDS(paste0(work_path, "pd_processed_PS_JaxE85.PCs.rds"))

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")
nodes = nodes[nodes$system == "Gastrulation",]

system_i = "Gastrulation"

pd_sub = pd[pd$group == "Pijuan_Sala",]
emb_sub = emb[pd$group == "Pijuan_Sala",]

pd_sub$celltype_new = as.vector(pd_sub$celltype_update)

pd_x = pd_sub %>% left_join(nodes %>% select(meta_group, celltype_new), by = "celltype_new")
pd_sub$meta_group = as.vector(pd_x$meta_group)

saveRDS(pd_sub[,c("cell_id","celltype_new","meta_group","day")], paste0(work_path, system_i, "_adata_scale.obs.rds"))

nodes_x = pd_sub %>% group_by(meta_group) %>% tally() %>% rename(celltype_num = n)
nodes_x = nodes %>% left_join(nodes_x)
nodes$celltype_num = as.vector(nodes_x$celltype_num)

kk = 10
k.param = kk + 1; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = emb_sub,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)

nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked[,-1]

x = data.frame(i = rep(1:nrow(nn_matrix), ncol(nn_matrix)),
               j = c(nn_matrix), stringsAsFactors = FALSE)

dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                           j = as.numeric(as.vector(x$j)),
                           x = 1)

dat_t = t(dat) + dat
x = data.frame(summary(dat_t))
x = x[x$x == 2 & x$i > x$j,]
x$x = NULL   ### x saves the MNN pairs

system_i = "Gastrulation"
saveRDS(x, paste0(work_path, system_i, ".MNN.rds"))

### the real observation

y = data.frame(i = 1:nrow(pd_sub),
               j = 1:nrow(pd_sub),
               meta_group = as.vector(pd_sub$meta_group), stringsAsFactors = FALSE)

dat = x %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j") %>%
    group_by(meta_group.x, meta_group.y) %>% tally() 

obs = dcast(dat, meta_group.x~meta_group.y)
rownames(obs) = as.vector(obs[,1])
obs = obs[,-1]
obs[is.na(obs)] = 0
obs = as.matrix(obs)

diag(obs) = 0

obs_x = obs + t(obs)
obs_y = as.vector(obs_x[upper.tri(obs_x)])

group = NULL
for(i in 2:nrow(obs_x)){
    print(i)
    for(j in 1:(i-1)){
        group = rbind(group, data.frame(x = colnames(obs_x)[j],
                                        y = rownames(obs_x)[i], stringsAsFactors = F))
    }
}

group$edge_num = obs_y

group = group %>% 
    left_join(nodes %>% rename(x = meta_group) %>% select(x, celltype_new, celltype_num), by = "x") %>% rename(x_name = celltype_new, x_size = celltype_num) %>%
    left_join(nodes %>% rename(y = meta_group) %>% select(y, celltype_new, celltype_num), by = "y") %>% rename(y_name = celltype_new, y_size = celltype_num)

group$min_size = if_else(group$x_size < group$y_size, group$x_size, group$y_size)
group$edge_num_norm = group$edge_num/log2(10*group$min_size)
group$min_size = NULL

system_i = "Gastrulation"
saveRDS(group, paste0(work_path, system_i, ".edges_new.rds"))

### output MNN pairs for review

edges = group

edges_2 = edges
edges_2$x = as.vector(edges$y); edges_2$y = as.vector(edges$x)
edges_2$x_size = as.vector(edges$y_size); edges_2$y_size = as.vector(edges$x_size)
edges_2$x_name = as.vector(edges$y_name); edges_2$y_name = as.vector(edges$x_name)

dat = rbind(edges, edges_2) %>% as.data.frame() %>%
    filter(edge_num != 0) %>% rename(MNN_pairs = edge_num, MNN_pairs_normalized = edge_num_norm) %>%
    group_by(x) %>% arrange(desc(MNN_pairs), .by_group = T) %>%
    as.data.frame()

write.table(dat, paste0(work_path, system_i, ".MNN_pairs.txt"), row.names=F, sep="\t", quote=F)
write.table(nodes[,c("celltype_new","celltype_num","meta_group","system")], paste0(work_path, system_i, ".nodes.txt"), row.names=F, sep="\t", quote=F)


