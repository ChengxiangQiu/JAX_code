
#################################
### Section - 8, Birth series ###
#################################

########################################################################################
### Identifying DEGs between E18.75 vs. P0 for the top 20 cell types shown in Fig.7b ###
########################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

top_20_celltypes = c("Hepatocytes",
                     "Alveolar Type 1 cells",
                     "Adipocyte cells (Cyp2e1+)",
                     "Glomerular endothelial cells",
                     "Mast cells",
                     "Liver sinusoidal endothelial cells",
                     "Dorsal root ganglion neurons",
                     "Proximal tubule cells",
                     "Endothelium",
                     "Definitive early erythroblasts (CD36-)",
                     "Brown adipocyte cells",
                     "Airway club cells",
                     "Alveolar Type 2 cells",
                     "Kupffer cells",
                     "Myofibroblasts",
                     "Border-associated macrophages",
                     "Midgut/Hindgut epithelial cells",
                     "Ependymal cells",
                     "Lymphatic vessel endothelial cells",
                     "Adipocyte progenitor cells")

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y")),]

res_all = NULL
for(kk in 1:20){
    
    celltype_i = top_20_celltypes[kk]
    print(celltype_i)
    
    pd = pd_all[pd_all$celltype_update == celltype_i & pd_all$day %in% c("E18.75","P0"),]
    min_num = min(table(pd$day))
    pd_sub = pd %>% group_by(day) %>% sample_n(min_num)
    pd = pd[pd$cell_id %in% as.vector(pd_sub$cell_id),]
    print(min_num)
    
    gene_count = doExtractData(pd, mouse_gene_sub)
    rownames(pd) = as.vector(pd$cell_id)
    
    obj = CreateSeuratObject(gene_count, meta.data = pd)
    obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
    genes_include = VariableFeatures(obj)
    Idents(obj) = as.vector(obj$day)
    
    res = FindMarkers(obj, ident.1 = "P0", ident.2 = "E18.75", features = genes_include)
    res = res %>% mutate(gene_ID = rownames(res)) %>% 
        left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>% as.data.frame() %>% filter(p_val_adj < 0.05)
    
    if(nrow(res) > 0){
        res$comparing = if_else(res$avg_logFC > 0, "up_in_P0", "down_in_P0")
        res$celltype_update = celltype_i
        res$p_val = NULL
        res_all = rbind(res_all, res)
    }
    
}

### Supplementary Table 20
write.csv(res_all, paste0(work_path, "top20_DEGs.csv"))


######################################################################################################################
### For the birth-series dataset, running regression to find genes which are correlated with C-section time series ###
######################################################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

trajectory_list = c("Hepatocytes", "Lung_and_airway", "Adipocytes")

for(kk in 1:3){
    trajectory_i = trajectory_list[kk]
    print(trajectory_i)
    
    pd = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
    mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y")),]
    gene_count = doExtractData(pd, mouse_gene_sub)
    
    cds = new_cell_data_set(gene_count,
                            cell_metadata = pd,
                            gene_metadata = mouse_gene_sub)
    
    if(trajectory_i == "Hepatocytes"){
        celltype_list = c("Hepatocytes")
    }
    
    if(trajectory_i == "Lung_and_airway"){
        celltype_list = c("Airway club cells", "Alveolar Type 1 cells", "Alveolar Type 2 cells")
    }
    
    if(trajectory_i == "Adipocytes"){
        celltype_list = c("Adipocyte cells (Cyp2e1+)", "Brown adipocyte cells")
    }
    
    for(celltype_i in celltype_list){
        
        print(celltype_i)
        cds_sub = cds[,cds$anno_subclustering == celltype_i]
        print(dim(cds_sub))
        
        obj = doObjectTransform(cds_sub, transform_to = "seurat")
        obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
        gene_include = VariableFeatures(obj)
        
        cds_sub = preprocess_cds(cds_sub, use_genes = gene_include)
        
        if(ncol(cds_sub) > 25000){
            x = data.frame(pData(cds_sub)) %>% group_by(day) %>% sample_n(5000)
            cds_sub = cds_sub[,pData(cds_sub)$cell_id %in% as.vector(x$cell_id)]
        }
        
        cds_test = cds_sub[fData(cds_sub)$gene_ID %in% gene_include,]
        gene_fits <- fit_models(cds_sub, model_formula_str = "~day")
        fit_coefs <- coefficient_table(gene_fits)
        emb_time_terms <- fit_coefs %>% filter(term == "day")
        res = emb_time_terms %>% filter (q_value < 0.05) %>%
            select(gene_short_name, term, q_value, estimate) %>% arrange(desc(estimate))
        
        celltype_name = doSimpleName(celltype_i)
        
        write.csv(res, paste0(work_path, trajectory_i, "_", celltype_name, "_monocle3.csv"))
        
    }
}


### searching if any overlap between the regression result and the DEGS bewteen P0 vs. E18.75

res_orig = read.csv(paste0(work_path, "top20_DEGs.csv"), row.names=1, as.is=T)

for(trajectory_i in trajectory_list){
    
    if(trajectory_i == "Hepatocytes"){
        celltype_list = c("Hepatocytes")
    }
    
    if(trajectory_i == "Lung_and_airway"){
        celltype_list = c("Airway club cells","Alveolar Type 1 cells","Alveolar Type 2 cells")
    }
    
    if(trajectory_i == "Adipocytes"){
        celltype_list = c("Adipocyte cells (Cyp2e1+)","Brown adipocyte cells")
    }
    
    for(celltype_i in celltype_list){
        res_sub = res_orig %>% filter(celltype_update == celltype_i)
        
        celltype_name = doSimpleName(celltype_i)
        
        res = read.csv(paste0(work_path, trajectory_i, "_", celltype_name, "_monocle3.csv"), row.names=1, as.is=T)
        
        res_up = res %>% filter(estimate > 0) %>% mutate(replicated_in_original = if_else(gene_short_name %in% as.vector(res_sub$gene_short_name[res_sub$comparing == "up_in_P0"]), "yes", "no"))
        res_down = res %>% filter(estimate < 0) %>% mutate(replicated_in_original = if_else(gene_short_name %in% as.vector(res_sub$gene_short_name[res_sub$comparing == "down_in_P0"]), "yes", "no"))
        
        res = rbind(res_up, res_down)
        write.csv(res, paste0(work_path, trajectory_i, "_", celltype_name, "_monocle3.csv"))
    }
}


###################################################################################
### Plotting gene expression for selected genes in those three major cell types ###
###################################################################################

### Fig. 6f

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

### plotting genes for the birth-series dataset

target_genes = c("Setdb2", "Pck1", "G6pc", "Lepr", "Ppargc1a", "Got1", "Vwf")

trajectory_i = "Hepatocytes"

pd = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
pd_sub = pd
mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]
gene_count_x = doExtractData(pd_sub, mouse_gene_sub)

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(t(as.matrix(gene_count_x))),
                gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                timepoint = rep(as.vector(pd_sub$day), nrow(gene_count_x)), stringsAsFactors = F)

### Fig. 6f (bottom left panel)
df %>%
    ggplot(aes(timepoint, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 


target_genes = c("Nr4a3","Irf4","Irs2","Ppargc1a","Ucp1","Acer2")

trajectory_i = "Adipocytes"

pd = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
pd_sub = pd[pd$anno_subclustering == "Brown adipocyte cells",]
mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]
gene_count_x = doExtractData(pd_sub, mouse_gene_sub)

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(t(as.matrix(gene_count_x))),
                gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                timepoint = rep(as.vector(pd_sub$day), nrow(gene_count_x)), stringsAsFactors = F)

### Fig. 6f (bottom middle panel)
df %>%
    ggplot(aes(timepoint, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 



target_genes = c("Pmvk","Acer2","Cap2","Tox")

trajectory_i = "Lung_and_airway"

pd = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
pd_sub = pd[pd$anno_subclustering == "Alveolar Type 1 cells",]
mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]
gene_count_x = doExtractData(pd_sub, mouse_gene_sub)

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(t(as.matrix(gene_count_x))),
                gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                timepoint = rep(as.vector(pd_sub$day), nrow(gene_count_x)), stringsAsFactors = F)

### Fig. 6f (bottom right panel)
df %>%
    ggplot(aes(timepoint, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 


### plotting genes in the original dataset, between E18.75 vs. P0
### Fig. 6f top three panels

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
pd_sub = pd_all[pd_all$celltype_update %in% c("Adipocyte cells (Cyp2e1+)", "Airway club cells", "Alveolar Type 1 cells",
                                              "Alveolar Type 2 cells", "Brown adipocyte cells", "Hepatocytes") &
                    pd_all$day %in% c("E18.75","P0"),]
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y")),]
gene_count = doExtractData(pd_sub, mouse_gene_sub)
gene_count = t(t(gene_count) / colSums(gene_count)) * 100000

### Repeat the below script by changing celltype_i to 1) Hepatocytes; 2) Brown adipocyte cells; 3) Alveolar Type 1 cells

celltype_i = "Hepatocytes"
target_genes = c("Setdb2", "Pck1", "G6pc", "Lepr", "Ppargc1a", "Got1", "Vwf")

celltype_i = "Brown adipocyte cells"
target_genes = c("Nr4a3","Irf4","Irs2","Ppargc1a","Ucp1","Acer2")

celltype_i = "Alveolar Type 1 cells"
target_genes = c("Pmvk","Acer2","Cap2","Tox")

mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]
gene_count_x = gene_count[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

gene_count_early = gene_count_x[,colnames(gene_count_x) %in% as.vector(pd_sub$cell_id)[pd_sub$day == "E18.75" & pd_sub$celltype_update == celltype_i]]
gene_count_late = gene_count_x[,colnames(gene_count_x) %in% as.vector(pd_sub$cell_id)[pd_sub$day == "P0" & pd_sub$celltype_update == celltype_i]]

### BARPLOT with mean value
df = rbind(data.frame(exp = Matrix::rowMeans(gene_count_early),
                      gene = rownames(gene_count_early), 
                      group = "E18.75", stringsAsFactors = F),
           data.frame(exp = Matrix::rowMeans(gene_count_late),
                      gene = rownames(gene_count_late), 
                      group = "P0", stringsAsFactors = F))

day_group_color_plate = c("E18.75" = "#7aa457",
                          "P0" = "#cb6a49")

df %>%
    ggplot(aes(x = gene, y = exp, fill = group)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    labs(x="", y="exp") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values=day_group_color_plate) +
    theme(axis.text.x = element_text(color="black", angle = 60), axis.text.y = element_text(color="black"))


