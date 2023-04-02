### Re Jay's comment - Are encoding TFs detectable in the inferred progenitors of the patterned neuroectoderm?

### For some selected spinal interneurons, can you plot the TF expressions in that 
### subset of cells as a function of time (as proxy for pseudotime). If they are turning 
### on late, that would be a clear story re: why it's not working that well

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

obj = readRDS(paste0(work_path, "/Neurons/Neurons_PCA_obj.rds"))
obj = CreateSeuratObject(GetAssayData(obj, slot = "counts"), meta = data.frame(obj[[]]))
Idents(obj) = as.vector(obj$celltype_sub_clustering)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

result = read.csv(paste0(work_path, "/Neurons/Neurons_interneurons_top_TFs.csv"))

result_sub = result %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 3)
result_sub = data.frame(result_sub, stringsAsFactors = F)

x = as.vector(obj$day)
x[obj$day == "E9"] = "E9.0"
x[obj$day == "E10"] = "E10.0"
x[obj$day == "E11"] = "E11.0"
x[obj$day == "E12"] = "E12.0"
x[obj$day == "E13"] = "E13.0"
x[obj$day == "E14"] = "E14.0"
x[obj$day == "E15"] = "E15.0"
x[obj$day == "E16"] = "E16.0"
x[obj$day == "E17"] = "E17.0"
x[obj$day == "E18"] = "E18.0"
obj$day = factor(as.vector(x), levels = day_list[day_list %in% x])

celltype_list = names(table(obj$celltype_sub_clustering))

for(celltype_i in celltype_list){
    print(celltype_i)
    
    result_i = result_sub %>% filter(cluster == celltype_i)
    
    obj_sub = subset(obj, subset = celltype_sub_clustering == celltype_i)
    
    name = gsub("Spinal ", "", celltype_i)
    name = gsub(" interneurons", "", name)
    
    for(i in c(1:3)){
        gene_ID = as.vector(result_i$gene_ID[i])
        gene_name = as.vector(result_i$gene_short_name[i])
        try(VlnPlot(object = obj_sub, features = gene_ID, group.by = 'day', pt.size = 0.1) +
            labs(title = paste0(name, " : ", gene_name)) +
            scale_fill_viridis(discrete=TRUE)  + NoLegend() +
            ggsave(paste0(work_path, "/Neurons/", name, "_", gene_name, ".2D_UMAP.day.png"), width = 6, height = 6, dpi = 300), silent = T)
    }
    
}


### Rather than having UMAP of the whole thing, can we just look at expression of the 
### TFs in the inferred backbone progenitors? For example is Isl1 expression highest 
### in dl3 inferred progreniors relative to other inferred progenitors? Hard to really 
### assess this in UMAP b/c this is such a small minority fo the cells. Violin or box plots better.


example_i = "Neurons_plus_big"

pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))


dat = readRDS(paste0(work_path, "/Neurons/", example_i, ".MNN.rds"))
pd_1 = dat[["pd_1"]]
pd_2 = dat[["pd_2"]]
nn_matrix = dat[["nn_matrix"]]
nn_matrix_2 = dat[["nn_matrix_2"]]

k.param = 10

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(pd_2$cell_id)[as.vector(nn_matrix[,i])])
}
datA = data.frame(A = rep(pd_1$cell_id, k.param),
                  B = c(resultA), 
                  value = 1, stringsAsFactors = F)

resultB = NULL
for(i in 1:k.param){
    print(i)
    resultB = cbind(resultB, as.vector(pd_1$cell_id)[as.vector(nn_matrix_2[,i])])
}
datB = data.frame(A = c(resultB),
                  B = rep(pd_2$cell_id, k.param), 
                  value = 2, stringsAsFactors = F)


dat = datB %>% left_join(datA, by = c("A" = "A", "B" = "B")) %>%
    filter(!is.na(value.y)) %>% select(A, B) %>%
    left_join(pd_1 %>% mutate(A = cell_id) %>% select(A, celltype_A = celltype), by = "A") %>%
    left_join(pd_2 %>% mutate(B = cell_id) %>% select(B, celltype_B = celltype), by = "B")



row_names_order = c("Spinal dI1 interneurons",
                    "Spinal V3 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons",
                    "Spinal V0 interneurons",
                    "Spinal V1 interneurons",
                    "Spinal V2a interneurons",
                    "Spinal V2b interneurons")

dat_m = dat %>% filter(celltype_A %in% row_names_order) %>% group_by(B, celltype_A) %>% tally() %>%
    group_by(B) %>% slice_max(order_by = n, n = 1, with_ties = F) %>%
    rename(cell_id = B, celltype_derive = celltype_A, freq = n)

pd_sub = pd %>% filter(cell_id %in% as.vector(dat_m$cell_id))


### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = names(table(as.vector(pd_sub$embryo_id)))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE])
}

gene_count = gene_count[,rownames(pd_sub)]

pd_x = pd_sub %>% left_join(dat_m, by = "cell_id")

pd_sub$celltype_derive = as.vector(pd_x$celltype_derive)

obj = CreateSeuratObject(gene_count, meta = pd_sub)
Idents(obj) = as.vector(obj$celltype_derive)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

result = read.csv(paste0(work_path, "/Neurons/Neurons_interneurons_top_TFs.csv"))
result_sub = result %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 3)
result_sub = data.frame(result_sub, stringsAsFactors = F)

for(i in 1:nrow(result_sub)){
    print(paste0(i,"/",nrow(result_sub)))
    
    celltype_i = as.vector(result_sub$cluster[i])
    gene_ID = as.vector(result_sub$gene_ID[i])
    gene_name = as.vector(result_sub$gene_short_name[i])

    name = gsub("Spinal ", "", celltype_i)
    name = gsub(" interneurons", "", name)
    
    try(VlnPlot(object = obj, features = gene_ID, group.by = "celltype_derive", pt.size = 1) +
            labs(title = paste0(name, " : ", gene_name)) +
            scale_fill_viridis(discrete=TRUE)  + NoLegend() +
            ggsave(paste0("~/share/", name, "_", gene_name, ".backbone.png"), width = 6, height = 6, dpi = 300), silent = T)


}


### we could take one or several interneurons, plotting all their key TFs, not just 
### the top three, over time. Ideally, we will see some TFs are early, some are late


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

obj = readRDS(paste0(work_path, "/Neurons/Neurons_PCA_obj.rds"))
obj = CreateSeuratObject(GetAssayData(obj, slot = "counts"), meta = data.frame(obj[[]]))
Idents(obj) = as.vector(obj$celltype_sub_clustering)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

x = as.vector(obj$day)
x[obj$day == "E9"] = "E9.0"
x[obj$day == "E10"] = "E10.0"
x[obj$day == "E11"] = "E11.0"
x[obj$day == "E12"] = "E12.0"
x[obj$day == "E13"] = "E13.0"
x[obj$day == "E14"] = "E14.0"
x[obj$day == "E15"] = "E15.0"
x[obj$day == "E16"] = "E16.0"
x[obj$day == "E17"] = "E17.0"
x[obj$day == "E18"] = "E18.0"
obj$day = factor(as.vector(x), levels = day_list[day_list %in% x])

data = GetAssayData(obj, slot = "data")

result = read.csv(paste0(work_path, "/Neurons/Neurons_interneurons_top_TFs.csv"))
celltype_list = names(table(result$cluster))

celltype_list = paste0("Spinal dI", c(1:5), " interneurons")

early_list = NULL

for(celltype_i in celltype_list){
    
    result_sub = result %>% filter(cluster == celltype_i)
    
    day_list_sub = names(table(as.vector(obj$day[obj$celltype_sub_clustering == celltype_i])))
    day_list_sub = as.vector(day_list[day_list %in% day_list_sub])
    
    exp = NULL
    for(day_i in day_list_sub){
        print(day_i)
        exp = cbind(exp, Matrix::rowMeans(data[rownames(data) %in% as.vector(result_sub$gene_ID),obj$celltype_sub_clustering == celltype_i &
                                                   obj$day == day_i, drop=FALSE]))
    }
    
    exp = exp[as.vector(result_sub$gene_ID),]
    rownames(exp) = as.vector(result_sub$gene_short_name)
    colnames(exp) = day_list_sub
    
    exp_row_max = sort(apply(exp, 1, which.max))
    exp = exp[names(exp_row_max),]
    
    early_list = rbind(early_list, data.frame(celltype = celltype_i,
                                              gene = names(exp_row_max)[exp_row_max == 1]))
    
    name = gsub("Spinal ", "", celltype_i)
    name = gsub(" interneurons", "", name)
    
    library("gplots")
    library(RColorBrewer)
    Colors=rev(brewer.pal(11,"Spectral"))
    Colors=colorRampPalette(Colors)(120)
    pdf(paste0("~/share/", name, "_heatmap.pdf"), 8, 8)
    heatmap.2(as.matrix(exp), 
              col=Colors, 
              scale="none", 
              Rowv = F, 
              Colv = F, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 0.5, 
              cexCol = 1,
              margins = c(5,5))
    dev.off()
    
    
}


### checking Slc32a1, Gad1 (GABA neurons) and Slc17a6 (Gluta neurons) expression across interneurons

try(VlnPlot(object = obj, features = "ENSMUSG00000070880", group.by = 'celltype_sub_clustering', pt.size = 0) +
        scale_color_manual(values=neuron_color_plate) + NoLegend() +
        ggsave("~/share/Gad1.png", width = 6, height = 6, dpi = 300), silent = T)

try(VlnPlot(object = obj, features = "ENSMUSG00000037771", group.by = 'celltype_sub_clustering', pt.size = 0) +
        scale_color_manual(values=neuron_color_plate) + NoLegend() +
        ggsave("~/share/Slc32a1.png", width = 6, height = 6, dpi = 300), silent = T)

try(VlnPlot(object = obj, features = "ENSMUSG00000030500", group.by = 'celltype_sub_clustering', pt.size = 0) +
        scale_color_manual(values=neuron_color_plate) + NoLegend() +
        ggsave("~/share/Slc17a6.png", width = 6, height = 6, dpi = 300), silent = T)


### I think we need to look at one more heatmap, which is of the expression fo 
### the earliest TFs (e.g. top rows of each heatmap, perhaps focusing on dl1-dl5) 
### and switch over to their inferred progenitors and ask whether they are high there? 
### question being is whether we can claim that the inferred progenitors 'choose' that 
### route b/c the TF combination is already on.


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

early_list = readRDS(paste0(work_path, "/Neurons/Neurons_interneurons_early_TFs.rds"))
early_list_gene_unique = table(as.vector(early_list$gene))
early_list = early_list[early_list$gene %in% names(early_list_gene_unique)[early_list_gene_unique == 1],]

example_i = "Neurons_plus_big"

pd = readRDS(paste0(work_path, "/Neurons/", example_i, "_adata_scale.obs.rds"))

dat = readRDS(paste0(work_path, "/Neurons/", example_i, ".MNN.rds"))
pd_1 = dat[["pd_1"]]
pd_2 = dat[["pd_2"]]
nn_matrix = dat[["nn_matrix"]]
nn_matrix_2 = dat[["nn_matrix_2"]]

k.param = 10

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(pd_2$cell_id)[as.vector(nn_matrix[,i])])
}
datA = data.frame(A = rep(pd_1$cell_id, k.param),
                  B = c(resultA), 
                  value = 1, stringsAsFactors = F)

resultB = NULL
for(i in 1:k.param){
    print(i)
    resultB = cbind(resultB, as.vector(pd_1$cell_id)[as.vector(nn_matrix_2[,i])])
}
datB = data.frame(A = c(resultB),
                  B = rep(pd_2$cell_id, k.param), 
                  value = 2, stringsAsFactors = F)


dat = datB %>% left_join(datA, by = c("A" = "A", "B" = "B")) %>%
    filter(!is.na(value.y)) %>% select(A, B) %>%
    left_join(pd_1 %>% mutate(A = cell_id) %>% select(A, celltype_A = celltype), by = "A") %>%
    left_join(pd_2 %>% mutate(B = cell_id) %>% select(B, celltype_B = celltype), by = "B")

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

### a progenitor cell from backbone is mutually exclusive to a given neuron type
dat_m = dat %>% filter(celltype_A %in% row_names_order) %>% group_by(B, celltype_A) %>% tally() %>%
    group_by(B) %>% slice_max(order_by = n, n = 1, with_ties = F) %>%
    rename(cell_id = B, celltype_derive = celltype_A, freq = n)

### a progenitor cell from backbone is not required to be exclusive to a given neuron type
dat_n = dat %>% filter(celltype_A %in% row_names_order) %>%
    rename(cell_id = B, celltype_derive = celltype_A)

pd_sub = pd %>% filter(cell_id %in% as.vector(dat_m$cell_id))

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

embryo_list = names(table(as.vector(pd_sub$embryo_id)))
gene_count = NULL
for(i in embryo_list){
    print(paste0(i, "/", length(embryo_list)))
    gene_count_tmp = readRDS(paste0(work_path, "/embryo/", i, "_gene_count.rds"))
    gene_count = cbind(gene_count, 
                       gene_count_tmp[rownames(gene_count_tmp) %in% rownames(mouse_gene_sub), 
                                      colnames(gene_count_tmp) %in% pd_sub$cell_id, drop=FALSE])
}

gene_count = gene_count[,rownames(pd_sub)]

saveRDS(list(gene_count = gene_count, pd_sub = pd_sub,
             dat_m = dat_m, dat_n = dat_n), paste0(work_path, "/Neurons/IN_backbone_prog.rds"))

### Let's think about the approach which not require progenitor to be exclusive - dat_n

obj = CreateSeuratObject(gene_count, meta = pd_sub)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
data = GetAssayData(obj, slot = "data")

gene_include = as.vector(mouse_gene$gene_ID[mouse_gene$gene_short_name %in% as.vector(early_list$gene)])
exp = NULL
for(celltype_i in row_names_order){
    cell_include = as.vector(dat_n %>% filter(celltype_derive == celltype_i) %>% pull(cell_id))
    exp = cbind(exp, Matrix::rowMeans(data[rownames(data) %in% gene_include,
                                           colnames(data) %in% cell_include, 
                                           drop=F]))
}

colnames(exp) = row_names_order
mouse_gene_sub = mouse_gene[rownames(exp),]
rownames(exp) = as.vector(mouse_gene_sub$gene_short_name)

gene_order = unique(as.vector(early_list$gene))
exp = exp[gene_order,]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0("~/share/", "progenitor_heatmap.pdf"), 8, 8)
heatmap.2(as.matrix(exp), 
          col=Colors, 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

obj$group = "IN_prog"
obj_processed = doClusterSeurat(obj)

dat_x = data.frame(obj_processed[[]]) %>% left_join(dat_m[,c("cell_id","celltype_derive")], by = "cell_id")
obj_processed$celltype_derive = as.vector(dat_x$celltype_derive)

pdf("~/share/IN_prog.pdf")
DimPlot(obj_processed, label=T)+NoLegend()
dev.off()

pdf("~/share/IN_prog_derive.pdf")
DimPlot(obj_processed, group.by = "celltype_derive", label=T)+NoLegend()
dev.off()

pdf("~/share/IN_prog_celltype.pdf")
DimPlot(obj_processed, group.by = "celltype_update", label=T)+NoLegend()
dev.off()


