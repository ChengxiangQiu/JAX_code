

### Re: #1, what if you instead scaled each TF across all its values for all timepoints 
### of all 11 spinal interneurons (even ones for which it is not nominated)? That way, 
### we would be better able to distinguish between TFs that start "on" and maybe dip 
### a little vs. those that start on and really shut off (while still accounting for
### the fact that some level of gene-specific normalization is necessary) 

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
celltype_list_all = names(table(result$cluster))
celltype_list = paste0("Spinal dI", c(1:5), " interneurons")
gene_list = as.vector(unique(result$gene_ID))

dat = data.frame(obj[[]]) %>% group_by(celltype_sub_clustering, day) %>% tally()
dat$col_id = paste0("col_", c(1:nrow(dat)))

exp = NULL
for(i in 1:nrow(dat)){
    print(paste0(i,"/",nrow(dat)))
    exp = cbind(exp, Matrix::rowMeans(data[rownames(data) %in% gene_list,
                                           obj$celltype_sub_clustering == as.vector(dat$celltype_sub_clustering)[i] &
                                           obj$day == as.vector(dat$day)[i], drop=FALSE]))
}

mouse_gene_sub = mouse_gene[rownames(exp),]
rownames(exp) = as.vector(mouse_gene_sub$gene_short_name)
colnames(exp) = as.vector(dat$col_id)

exp_scale = t(scale(t(exp)))

for(celltype_i in celltype_list_all){
    
    result_sub = result %>% filter(cluster == celltype_i)
    dat_sub = dat %>% filter(celltype_sub_clustering == celltype_i)
    
    exp_sub = exp_scale[as.vector(result_sub$gene_short_name), as.vector(dat_sub$col_id)]
    colnames(exp_sub) = as.vector(dat_sub$day)
    
    exp_row_max = sort(apply(exp_sub, 1, which.max))
    exp_sub = exp_sub[names(exp_row_max),]
    
    name = gsub("Spinal ", "", celltype_i)
    name = gsub(" interneurons", "", name)
    
    library("gplots")
    library(RColorBrewer)
    Colors=rev(brewer.pal(11,"Spectral"))
    Colors=colorRampPalette(Colors)(120)
    pdf(paste0("~/share/", name, "_heatmap.pdf"), 8, 8)
    heatmap.2(as.matrix(exp_sub), 
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



### my own idea
### We select their progenitor cells, to check if they are really separated by reclustering, and then identify the top highly expressed TFs for each dI1-5
### For each dI1-5, identify the TF which expresses at the early stages
### Checking the two list, if any overlap


#####################
### first, focusing the prog. of individual dI1-5 neurons
#####################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

dat = readRDS(paste0(work_path, "/Neurons/IN_backbone_prog.rds"))

dat_n = dat[["dat_n"]]

obj = CreateSeuratObject(dat[["gene_count"]], meta.data = dat[["pd_sub"]])
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

mouse_TF = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/TF_list/AnimalTFDB/Mus_musculus_TF.txt", sep="\t", header=T, as.is=T)
mouse_TF_overlap = intersect(as.vector(mouse_TF$Ensembl), rownames(obj))

res_all = NULL
for(celltype_i in row_names_order){
    print(celltype_i)
    prog_i = dat_n %>% filter(celltype_derive == celltype_i) %>% pull(cell_id)
    obj_i = obj
    Idents(obj_i) = if_else(obj_i$cell_id %in% prog_i, "yes", "no")
    print(table(Idents(obj_i)))
    res_i = FindMarkers(obj_i, ident.1 = "yes", ident.2 = "no", features = mouse_TF_overlap, only.pos = T)
    res_i = res_i %>% mutate(gene_ID = rownames(res_i), cluster = celltype_i) %>%
        left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>%
        filter(p_val_adj < 0.05)
    res_all = rbind(res_all, res_i)
}

write.csv(res_all, paste0(work_path, "/Neurons/Neurons_interneurons_prog_TFs.csv"))


######################
### second, focusing on dI1-5 neurons, which TF expresses at the earliest stages (top 500 cells)
######################

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

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

dat = readRDS(paste0(work_path, "/Neurons/", "Neurons_plus_big", ".MNN.rds"))
pd_1 = dat[["pd_1"]]
pd_sub = pd_1[pd_1$celltype_update %in% row_names_order,]

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

obj_early = CreateSeuratObject(gene_count, meta.data = pd_sub)
obj_early$group = "early"
saveRDS(obj_early, paste0(work_path,"/Neurons/IN_early_cells.rds"))

obj$celltype_update = as.vector(obj$celltype_sub_clustering)
obj$keep = obj$celltype_update %in% row_names_order & 
           (!obj$cell_id %in% pd_sub$cell_id)
obj_sub = subset(obj, subset = keep)
obj_sub$group = "late"

obj_merge = merge(obj_early, obj_sub)
obj_merge = NormalizeData(obj_merge, normalization.method = "LogNormalize", scale.factor = 10000)

mouse_TF = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/TF_list/AnimalTFDB/Mus_musculus_TF.txt", sep="\t", header=T, as.is=T)
mouse_TF_overlap = intersect(as.vector(mouse_TF$Ensembl), rownames(obj_merge))

res_all = NULL
for(celltype_i in row_names_order){
    print(celltype_i)
    obj_i = subset(obj_merge, subset = celltype_update == celltype_i)
    Idents(obj_i) = as.vector(obj_i$group)
    res_i = FindMarkers(obj_i, ident.1 = "early", ident.2 = "late", features = mouse_TF_overlap, only.pos = T)
    res_i = res_i %>% mutate(gene_ID = rownames(res_i), cluster = celltype_i) %>%
        left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>%
        filter(p_val_adj < 0.05)
    res_all = rbind(res_all, res_i)
}

write.csv(res_all, paste0(work_path, "/Neurons/Neurons_interneurons_early_TFs.csv"))


###################
### Third, checking if any overlaps
###################

res_top = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_top_TFs.csv"), header=T, as.is=T)
res_early = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_early_TFs.csv"), header=T, as.is=T)
res_prog = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_prog_TFs.csv"), header=T, as.is=T)

res_overlap = NULL
for(celltype_i in row_names_order){
    x1 = res_top %>% filter(cluster == celltype_i) %>% pull(gene_short_name)
    x2 = res_early %>% filter(cluster == celltype_i) %>% pull(gene_short_name)
    x3 = res_prog %>% filter(cluster == celltype_i) %>% pull(gene_short_name)
    
    x = intersect(x1, intersect(x2, x3))
    res_overlap = rbind(res_overlap, data.frame(cluster = celltype_i,
                                                gene_short_name = x, stringsAsFactors = F))
}

res_sub = res_top[res_top$cluster %in% row_names_order,]

write.csv(res_sub, paste0(work_path,"/Neurons/Neurons_interneurons_top_TFs_sub.csv"))

write.csv(res_overlap, paste0(work_path,"/Neurons/Neurons_interneurons_overlap_TFs.csv"))

##############
### Trying Jay's suggestion 

#I am not sure if I think we should go with violin plots. We may not even need to 
#do anything other than provide a list. But can we try these things quickly?
    
#1) Make a heatmap like Fig. 4f, but with these 18 TFs and only using the first 
#500 cells from each of dl1-dl5 for each column

#2) Same thing but use the inferred progenitors for dl1-dl5 for each column



work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

res_overlap = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_overlap_TFs.csv"))
res_x = res_overlap %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_short_name")
res_overlap$gene_ID = as.vector(res_x$gene_ID)

### first, making heatmap with prog cells

dat = readRDS(paste0(work_path, "/Neurons/IN_backbone_prog.rds"))

dat_n = dat[["dat_n"]]

obj = CreateSeuratObject(dat[["gene_count"]], meta.data = dat[["pd_sub"]])
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
data = GetAssayData(obj, slot = "data")

exp = NULL
for(celltype_i in row_names_order){
    cell_include = as.vector(dat_n %>% filter(celltype_derive == celltype_i) %>% pull(cell_id))
    exp = cbind(exp, Matrix::rowMeans(data[rownames(data) %in% as.vector(res_overlap$gene_ID),
                                           colnames(data) %in% cell_include, 
                                           drop=F]))
}

colnames(exp) = row_names_order
exp = exp[as.vector(res_overlap$gene_ID),]
rownames(exp) = as.vector(res_overlap$gene_short_name)

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0("~/share/", "progenitor_heatmap.pdf"), 8, 8)
heatmap.2(as.matrix(t(exp)), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

### second, plot the gene exp across timepoints


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

res_overlap = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_overlap_TFs.csv"))
res_x = res_overlap %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_short_name")
res_overlap$gene_ID = as.vector(res_x$gene_ID)

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

obj_early = readRDS(paste0(work_path,"/Neurons/IN_early_cells.rds"))
obj_early = NormalizeData(obj_early, normalization.method = "LogNormalize", scale.factor = 10000)

obj$celltype_update = as.vector(obj$celltype_sub_clustering)
obj$keep = obj$celltype_update %in% row_names_order & 
    (!obj$cell_id %in% obj_early$cell_id)
obj_sub = subset(obj, subset = keep)
obj_sub$group = "late"
obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)

dat_early = GetAssayData(obj_early, slot = "data")[as.vector(res_overlap$gene_ID),]
dat_late = GetAssayData(obj_sub, slot = "data")[as.vector(res_overlap$gene_ID),]
pd_sub = data.frame(obj_sub[[]])

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
pd_x = pd_sub %>% left_join(pd_all[,c("cell_id","somite_count")], by = "cell_id")
pd_sub$somite_count = as.vector(pd_x$somite_count)

day_value = as.vector(pd_sub$somite_count)
day_value[pd_sub$day == "E10.25"] = "35 somites"
day_value[pd_sub$day == "E10.5"]  = "36 somites"
day_value[pd_sub$day == "E10.75"] = "37 somites"
day_value[pd_sub$day == "E11.0"]  = "38 somites"
day_value[pd_sub$day == "E11.25"] = "39 somites"
day_value[pd_sub$day == "E11.5"]  = "40 somites"
day_value[pd_sub$day == "E11.75"] = "41 somites"
day_value[pd_sub$day == "E12.0"]  = "42 somites"
day_value[pd_sub$day == "E12.25"] = "43 somites"
day_value[pd_sub$day == "E12.5"]  = "44 somites"
day_value[pd_sub$day == "E12.75"] = "45 somites"
pd_sub$somite_count = as.vector(day_value)

day_x = data.frame(somite_count = paste0(c(27:45), " somites"),
                   day_value = c(1:length(c(27:45))))
day_y = pd_sub %>% left_join(day_x, by = "somite_count")
pd_sub$day_value = as.vector(day_y$day_value)

exp = NULL
for(i in 1:nrow(res_overlap)){
    print(i)
    celltype_i = as.vector(res_overlap$cluster)[i]
    gene_i = as.vector(res_overlap$gene_ID)[i]
    
    exp_i = NULL
    exp_i = c(exp_i, mean(as.vector(dat_early[gene_i, obj_early$celltype_update == celltype_i])))
    
    pd_i = pd_sub[pd_sub$celltype_update == celltype_i,]
    print(table(pd_i$day_value))
    pd_i = pd_i[order(pd_i$day_value),]
    bin_size = floor(nrow(pd_i)/20)+1
    for(j in 1:20){
        
        bin_coor = ((j-1)*bin_size + 1):(j*bin_size)
        
        if(j == 20){
            bin_coor = ((j-1)*bin_size + 1):nrow(pd_i)
        }
        
        pd_j = pd_i[bin_coor,]
        exp_i = c(exp_i, mean(as.vector(dat_late[gene_i, colnames(dat_late) %in% as.vector(pd_j$cell_id)])))
    }
    exp = cbind(exp, exp_i)
}

rownames(exp) = c("Early 500", paste0("bin_", 1:20))
colnames(exp) = as.vector(res_overlap$gene_short_name)

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0("~/share/", "timepoint_heatmap.pdf"), 8, 8)
heatmap.2(as.matrix(exp), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 1,
          margins = c(5,5))
dev.off()






### For genes


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

res_overlap = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_prog_genes.csv"))
res_x = res_overlap %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_short_name")
res_overlap$gene_ID = as.vector(res_x$gene_ID)

### first, making heatmap with prog cells

dat = readRDS(paste0(work_path, "/Neurons/IN_backbone_prog.rds"))

dat_n = dat[["dat_n"]]

obj = CreateSeuratObject(dat[["gene_count"]], meta.data = dat[["pd_sub"]])
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
data = GetAssayData(obj, slot = "data")

exp = NULL
for(celltype_i in row_names_order){
    cell_include = as.vector(dat_n %>% filter(celltype_derive == celltype_i) %>% pull(cell_id))
    exp = cbind(exp, Matrix::rowMeans(data[rownames(data) %in% as.vector(res_overlap$gene_ID),
                                           colnames(data) %in% cell_include, 
                                           drop=F]))
}

colnames(exp) = row_names_order
exp = exp[as.vector(res_overlap$gene_ID),]
rownames(exp) = as.vector(res_overlap$gene_short_name)

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0("~/share/", "progenitor_heatmap_2.pdf"), 8, 8)
heatmap.2(as.matrix(t(exp)), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

### second, plot the gene exp across timepoints


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")
save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_data"

row_names_order = c("Spinal dI1 interneurons",
                    "Spinal dI2 interneurons",
                    "Spinal dI3 interneurons",
                    "Spinal dI4 interneurons",
                    "Spinal dI5 interneurons")

res_overlap = read.csv(paste0(work_path,"/Neurons/Neurons_interneurons_prog_genes.csv"))
res_x = res_overlap %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_short_name")
res_overlap$gene_ID = as.vector(res_x$gene_ID)

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

obj$celltype_update = as.vector(obj$celltype_sub_clustering)
obj$keep = obj$celltype_update %in% row_names_order
obj_sub = subset(obj, subset = keep)
obj_sub$group = "late"
obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)

dat_late = GetAssayData(obj_sub, slot = "data")[as.vector(res_overlap$gene_ID),]
pd_sub = data.frame(obj_sub[[]])

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
pd_x = pd_sub %>% left_join(pd_all[,c("cell_id","somite_count")], by = "cell_id")
pd_sub$somite_count = as.vector(pd_x$somite_count)

day_value = as.vector(pd_sub$somite_count)
day_value[pd_sub$day == "E10.25"] = "35 somites"
day_value[pd_sub$day == "E10.5"]  = "36 somites"
day_value[pd_sub$day == "E10.75"] = "37 somites"
day_value[pd_sub$day == "E11.0"]  = "38 somites"
day_value[pd_sub$day == "E11.25"] = "39 somites"
day_value[pd_sub$day == "E11.5"]  = "40 somites"
day_value[pd_sub$day == "E11.75"] = "41 somites"
day_value[pd_sub$day == "E12.0"]  = "42 somites"
day_value[pd_sub$day == "E12.25"] = "43 somites"
day_value[pd_sub$day == "E12.5"]  = "44 somites"
day_value[pd_sub$day == "E12.75"] = "45 somites"
pd_sub$somite_count = as.vector(day_value)

day_x = data.frame(somite_count = paste0(c(27:45), " somites"),
                   day_value = c(1:length(c(27:45))))
day_y = pd_sub %>% left_join(day_x, by = "somite_count")
pd_sub$day_value = as.vector(day_y$day_value)

bin_num = 10
exp = NULL
for(i in 1:nrow(res_overlap)){
    print(i)
    celltype_i = as.vector(res_overlap$cluster)[i]
    gene_i = as.vector(res_overlap$gene_ID)[i]
    
    exp_i = NULL

    pd_i = pd_sub[pd_sub$celltype_update == celltype_i,]
    print(table(pd_i$day_value))
    pd_i = pd_i[order(pd_i$day_value),]
    bin_size = floor(nrow(pd_i)/bin_num)+1
    for(j in 1:bin_num){
        
        bin_coor = ((j-1)*bin_size + 1):(j*bin_size)
        
        if(j == bin_num){
            bin_coor = ((j-1)*bin_size + 1):nrow(pd_i)
        }
        
        pd_j = pd_i[bin_coor,]
        exp_i = c(exp_i, mean(as.vector(dat_late[gene_i, colnames(dat_late) %in% as.vector(pd_j$cell_id)])))
    }
    exp = cbind(exp, exp_i)
}

rownames(exp) = paste0("bin_", 1:bin_num)
colnames(exp) = as.vector(res_overlap$gene_short_name)

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0("~/share/", "timepoint_heatmap_2.pdf"), 8, 8)
heatmap.2(as.matrix(exp), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

