#############################################################
### extracting data used to perform spatial analysis using Tangram with stereo-seq data 
#############################################################


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mosta/sc_data"
source("/net/gs/vol1/home/cxqiu/work/scripts/tome/utils.R")

pd_all = readRDS(paste0(work_path, "/mtx/adata_scale.obs.rds"))
pd = readRDS(paste0(work_path, "/mtx/example/LPM_adata_scale.obs.rds"))

rownames(pd_all) = as.vector(pd_all$cell_id)
rownames(pd) = as.vector(pd$cell_id)

pd_all = pd_all[rownames(pd),]
pd$embryo_id = as.vector(pd_all$embryo_id)

day_include = list()
day_include[["E95"]] = c("E9.25","E9.5","E9.75")
day_include[["E105"]] = c("E10.25","E10.5","E10.75")
day_include[["E115"]] = c("E11.25","E11.5","E11.75")
day_include[["E125"]] = c("E12.25","E12.5","E12.75")
day_include[["E135"]] = c("E13.25","E13.5","E13.75")
day_include[["E145"]] = c("E14.25","E14.333","E14.75")
day_include[["E155"]] = c("E15.25","E15.5","E15.75")
day_include[["E165"]] = c("E16.25","E16.5","E16.75")

day_list = names(day_include)

for(i in 1:length(day_list)){
    day_i = day_list[i]; print(day_i)
    pd_sub = pd[pd$day %in% day_include[[day_i]],]
    
    ### downsampling each cell type to 5000 cells 
    ### & removing celltypes with less than 200 cells
    x_table = table(pd_sub$celltype_sub_clustering)
    x_1 = names(x_table)[x_table >= 200 & x_table <= 5000]
    x_2 = names(x_table)[x_table > 5000]
    
    pd_sub_1 = pd_sub[pd_sub$celltype_sub_clustering %in% x_1,]
    pd_sub_2 = pd_sub[pd_sub$celltype_sub_clustering %in% x_2,] %>% 
        group_by(celltype_sub_clustering) %>%
        sample_n(5000)
    pd_sub_2 = pd_sub[pd_sub$cell_id %in% pd_sub_2$cell_id,]
    pd_sub = rbind(pd_sub_1, pd_sub_2)
    
    embryo_list = names(table(pd_sub$embryo_id))
    print(embryo_list)
    gene_count = NULL
    for(j in embryo_list){
        count_j = readRDS(paste0(work_path, "/embryo/", j, "_gene_count.rds"))
        gene_count = cbind(gene_count, count_j[,colnames(count_j) %in% rownames(pd_sub)])
    }
    pd_sub = pd_sub[colnames(gene_count),]
    pd_sub = pd_sub[,c("day", "batch", "embryo_id", "celltype_sub_clustering")]
    colnames(pd_sub) = c("day", "sequencing_batch", "embryo_id", "celltype")
    
    gene_keep = mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA') & 
        mouse_gene$chr %in% paste0("chr", c(1:19, "M"))
    gene_count = gene_count[rownames(gene_count) %in% rownames(mouse_gene)[gene_keep],]
    
    df_gene = mouse_gene[rownames(gene_count),]
    df_gene$mean_exp = Matrix::rowMeans(gene_count)
    df_gene_x = df_gene %>% group_by(gene_short_name) %>% slice_max(order_by = mean_exp, n = 1, with_ties = F)
    gene_count = gene_count[rownames(gene_count) %in% df_gene_x$gene_ID,]
    df_gene = mouse_gene[rownames(gene_count),]
    row.names(gene_count) = rownames(df_gene) = as.vector(df_gene$gene_short_name)
    
    Matrix::writeMM(t(gene_count), paste0(save_path, "/", day_i, ".gene_count.mtx"))
    write.csv(df_gene, paste0(save_path, "/", day_i, ".df_gene.csv"))
    write.csv(pd_sub, paste0(save_path, "/", day_i, ".df_cell.csv"))
}





