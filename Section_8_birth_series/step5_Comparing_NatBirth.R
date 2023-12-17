
#################################
### Section - 8, Birth series ###
#################################

##########################################################################################################
### 2D UMAP of subclustering results for adipocytes and lung & airway (adding Natural birthed samples) ###
##########################################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

trajectory_list = c("Hepatocytes", "Adipocytes", "Lung_and_airway")

### Extended Data Fig. 12f

for(trajectory_i in trajectory_list){
    
    print(trajectory_i)
    
    df = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections_NatBirth.obs.rds"))

    rep_id = as.vector(df$day)
    rep_id[df$embryo_id == "embryo_76"] = "NatBirth_rep1"
    rep_id[df$embryo_id == "embryo_77"] = "NatBirth_rep2"
    rep_id[df$embryo_id == "embryo_78"] = "NatBirth_rep3"
    df$rep_id = as.vector(rep_id)
    
    for(i in names(table(df$rep_id))){
        
        try(ggplot() +
            geom_point(data = df, aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=1.5, alpha = 0.3) +
            geom_point(data = df %>% filter(rep_id == i),
                       aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=1.5) +
            theme_void() +
            scale_color_manual(values=birth_color_plate) +
            theme(legend.position="none") + 
            ggsave(paste0(work_path, "birth_", trajectory_i, "_NatBirth_", i, ".png"), width = 4, height = 4, dpi = 300))
    }
    
}


###############################################################################
### Identifying the neighbors for each NatBirth samples in the co-embedding ###
###############################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

trajectory_list = c("Hepatocytes", "Adipocytes", "Lung_and_airway")

### Extended Data Fig. 12e

for(trajectory_i in trajectory_list){
    
    print(trajectory_i)
    
    df = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections_NatBirth.obs.rds"))
    
    rep_id = as.vector(df$day)
    rep_id[df$embryo_id == "embryo_76"] = "NatBirth_rep1"
    rep_id[df$embryo_id == "embryo_77"] = "NatBirth_rep2"
    rep_id[df$embryo_id == "embryo_78"] = "NatBirth_rep3"
    df$rep_id = as.vector(rep_id)
    
    emb = read.csv(paste0(work_path, "adata_", trajectory_i, "_NatBirth.PCs.csv"), header=F, as.is=T)
    colnames(emb) = paste0("PC_", 1:30)
    rownames(emb) = rownames(df) = as.vector(df$cell_id)
    emb = as.matrix(emb)
    
    result = list()
    for(kk in 1:3){
        rep_i = paste0("NatBirth_rep", kk); print(rep_i)
        
        df_1 = df[df$rep_id == rep_i,]
        emb_1 = emb[df$rep_id == rep_i,]
        
        df_2 = df[df$rep_id != rep_i,]
        df_2_x = df_2 %>% group_by(rep_id) %>% sample_n(5000)
        df_2 = df_2[df_2$cell_id %in% df_2_x$cell_id,]
        emb_2 = emb[as.vector(df_2$cell_id),]
        
        k.param = 10; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
        nn.ranked = Seurat:::NNHelper(
            data = emb_2,
            query = emb_1,
            k = k.param,
            method = nn.method,
            searchtype = "standard",
            eps = nn.eps,
            metric = annoy.metric)
        nn.ranked = Indices(object = nn.ranked)
        nn_matrix = nn.ranked
        
        resultA = NULL
        for(i in 1:k.param){
            print(i)
            resultA = cbind(resultA, as.vector(df_2$rep_id)[as.vector(nn_matrix[,i])])
        }
        
        result[[rep_i]] = table(c(resultA))
    }
    
    dat_1 = data.frame(target_id = names(result[[1]]), num = as.vector(result[[1]]))
    dat_1 = rbind(dat_1, data.frame(target_id = "NatBirth_rep1", num = 0))
    dat_1$rep_id = "NatBirth_rep1"
    
    dat_2 = data.frame(target_id = names(result[[2]]), num = as.vector(result[[2]]))
    dat_2 = rbind(dat_2, data.frame(target_id = "NatBirth_rep2", num = 0))
    dat_2$rep_id = "NatBirth_rep2"
    
    dat_3 = data.frame(target_id = names(result[[3]]), num = as.vector(result[[3]]))
    dat_3 = rbind(dat_3, data.frame(target_id = "NatBirth_rep3", num = 0))
    dat_3$rep_id = "NatBirth_rep3"
    
    dat = rbind(dat_1, dat_2, dat_3)
    
    dat$target_id = factor(dat$target_id, levels = c(paste0("NatBirth_rep",1:3), paste0("Csection_",c(0,20,40,60,80),"m")))
    
    p = dat %>%
        ggplot(aes(x=target_id, y = num, fill = target_id))+
        geom_bar(stat="identity") + facet_grid(rep_id ~ .) +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        labs(x="",y="# of kNNs") +
        scale_fill_manual(values=birth_color_plate) +
        theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color="black"))
    
}


##################################################################
###  DEGs between NatBirth and C-section samples (20m and 40m) ###
##################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "pd_birth.rds"))
pd$anno = as.vector(pd$major_trajectory)

pd_sub = pd %>% filter(day %in% c("NatBirth","Csection_20m","Csection_40m"))
pd_sub$embryo_group = paste0(pd_sub$anno, "_", pd_sub$embryo_group)
x_table = table(pd_sub$embryo_group)
pd_sub_1 = pd_sub[pd_sub$embryo_group %in% names(x_table)[x_table > 10000],]
pd_sub_2 = pd_sub[pd_sub$embryo_group %in% names(x_table)[x_table <= 10000],]
pd_sub_1_x = pd_sub_1 %>% group_by(embryo_group) %>% sample_n(10000)
pd_sub_1 = pd_sub_1[pd_sub_1$cell_id %in% pd_sub_1_x$cell_id,]
df = rbind(pd_sub_1, pd_sub_2)

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y")),]
gene_count = doExtractData(df, mouse_gene_sub)
obj = CreateSeuratObject(gene_count, meta.data = df)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

anno_list = names(table(pd$anno))
res_all = NULL

for(i in 1:length(anno_list)){
    anno_i = anno_list[i]; print(anno_i)
    obj_sub = subset(obj, subset = anno == anno_i)
    Idents(obj_sub) = as.vector(obj_sub$day)
    obj_sub = FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = 5000)
    genes_include = VariableFeatures(obj_sub)
    
    res = FindMarkers(obj_sub, ident.1 = "NatBirth", features = genes_include)
    res = res %>% mutate(gene_ID = rownames(res)) %>% 
        left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>% as.data.frame() %>% filter(p_val_adj < 0.05)
    res$high_in_which = if_else(res$avg_logFC > 0, "Up_in_NatBirth", "Down_in_NatBirth")
    res$p_val = NULL
    res$major_cell_cluster = anno_i
    
    res_all = rbind(res_all, res)
}

write.csv(res_all, paste0(work_path, "adata_major_cell_cluster_NatBirth_DEGs.csv"))






