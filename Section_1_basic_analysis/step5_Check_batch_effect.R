
###################################
### Section - 1, Basic analysis ###
###################################

################################################################
### checking potential batch effects (Supplementary Fig. 3b) ###
################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")

work_path = "./"

pd = readRDS(paste0(work_path, "df_cell.rds"))
### n = 11,441,407 cells

x = as.vector(pd$day)
x[pd$day == "E8.0-E8.5"] = "E8.5"
pd$day = as.vector(x)

pd_sub_1 = pd[pd$sequencing_batch == "run_22" &
                  pd$group == "E14",]
pd_sub_2 = pd[pd$sequencing_batch == "run_22" &
                  pd$group == "E17",]
pd_sub_2 = pd_sub_2[sample(1:nrow(pd_sub_2), 500000),]
pd_sub_3 = pd[pd$sequencing_batch %in% c("run_18", "run_13") &
                  pd$group == "E14",]
pd_sub_3 = pd_sub_3[sample(1:nrow(pd_sub_3), 500000),]
pd_sub = rbind(pd_sub_2, pd_sub_3)

neighbors <- get.knnx(pd_sub[,c("UMAP_1","UMAP_2","UMAP_3")], pd_sub_1[,c("UMAP_1","UMAP_2","UMAP_3")], k = 10)$nn.index

res_1 = neighbors > 500000
res_2 = apply(res_1, 1, sum)/10
res_3 = 1 - res_2

df = rbind(data.frame(pct = res_2, resource = "same_time_window"),
           data.frame(pct = res_3, resource = "same_batch"))

p = ggplot(df, aes(x=pct, fill=resource)) +
    geom_histogram(position="dodge")+
    theme(legend.position="top") +
    scale_fill_brewer(palette="Set2") +
    theme_classic()

pd_sub_1 = pd[pd$sequencing_batch == "run_19" &
                  pd$group == "E13",]
pd_sub_2 = pd[pd$sequencing_batch == "run_19" &
                  pd$group == "E10",]
pd_sub_2 = pd_sub_2[sample(1:nrow(pd_sub_2), 500000),]
pd_sub_3 = pd[pd$sequencing_batch %in% c("run_13") &
                  pd$group == "E13",]
pd_sub_3 = pd_sub_3[sample(1:nrow(pd_sub_3), 500000),]
pd_sub = rbind(pd_sub_2, pd_sub_3)

neighbors <- get.knnx(pd_sub[,c("UMAP_1","UMAP_2","UMAP_3")], pd_sub_1[,c("UMAP_1","UMAP_2","UMAP_3")], k = 10)$nn.index

res_1 = neighbors > 500000
res_2 = apply(res_1, 1, sum)/10
res_3 = 1 - res_2

df = rbind(data.frame(pct = res_2, resource = "same_time_window"),
           data.frame(pct = res_3, resource = "same_batch"))

p = ggplot(df, aes(x=pct, fill=resource)) +
    geom_histogram(position="dodge")+
    theme(legend.position="top") +
    scale_fill_brewer(palette="Set2") +
    theme_classic()


##################################################
### Integrating cells from adjacent timepoints ###
##################################################

### Here, to check if any potential batch effects, we integrate cells from adjacent 
### timepoints but profiled by different sci-RNA-seq3 experiments.

### Supplementary Figure 4a-b

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "/df_cell.rds"))

mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

day_list = c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
             "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
             "E12.25", "E12.5", "E12.75", "E13.0", "E13.25", "E13.5", "E13.75", 
             "E14.0", "E14.25", "E14.333", "E14.75", "E15.0", "E15.25", "E15.5", 
             "E15.75", "E16.0", "E16.25", "E16.5", "E16.75", "E17.0", "E17.25", 
             "E17.5", "E17.75", "E18.0", "E18.25", "E18.5", "E18.75", "P0")

batch_list = c()

for(kk in 1:(length(day_list)-1)){
    day_1 = day_list[kk]
    day_2 = day_list[kk+1]
    print(paste0(day_1, ":", day_2))
    pd_sub = pd_all[pd_all$day %in% c(day_1, day_2),]
    if(length(table(pd_sub$sequencing_batch)) != 1){
        
        ### downsampling if necessary to save computational time
        pd_sub = pd_sub %>% group_by(sequencing_batch) %>% slice_sample(n = 100000) %>% as.data.frame()
        rownames(pd_sub) = as.vector(pd_sub$cell_id)
        
        gene_count = doExtractData(pd_sub, mouse_gene_sub)
        
        name = paste0(day_1, "_", day_2)
        Matrix::writeMM(t(gene_count), paste0(work_path, name, ".gene_count.mtx"))
        write.csv(pd_sub, paste0(work_path, name, ".df_cell.csv"))
        write.csv(mouse_gene_sub, paste0(work_path, name, ".df_gene.csv"))
        
        batch_list = c(batch_list, name)
    }
}

### run_23_A and run_23_B

pd_sub = pd_all[pd_all$sequencing_batch %in% c("run_23", "run_27"),]
pd_sub = pd_sub %>% group_by(sequencing_batch) %>% slice_sample(n = 100000) %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

gene_count = doExtractData(pd_sub, mouse_gene_sub)

name = "run_23_A_B"
Matrix::writeMM(t(gene_count), paste0(work_path, name, ".gene_count.mtx"))
write.csv(pd_sub, paste0(work_path, name, ".df_cell.csv"))
write.csv(mouse_gene_sub, paste0(work_path, name, ".df_gene.csv"))

batch_list = c(batch_list, name)

write.table(batch_list, paste0(work_path, "batch_list.txt"), row.names=F, col.names=F, quote=F, sep="\t")

### Performing co-embedding using Scanpy
### python Integrating_adjacent_timepoints.py

for(batch_id in batch_list){
    print(batch_id)
    pd = read.csv(paste0(work_path, batch_id, "_adata_scale.obs.csv"))
    
    tmp = as.vector(pd$sequencing_batch)
    tmp[pd$sequencing_batch == "run_23"] = "run_23_A"
    tmp[pd$sequencing_batch == "run_27"] = "run_23_B"
    pd$sequencing_batch = as.vector(tmp)
    
    run_list = names(table(pd$sequencing_batch))
    
    run_color = c("#cb5362","#6ea84e","#9d6cc1","#bb873c", "#3fadaf")
    
    for(i in 1:length(run_list)){
        run_i = run_list[i]; print(run_i)
        try(ggplot() +
                geom_point(data = pd[pd$sequencing_batch != run_i,], aes(x = UMAP_2d_1, y = UMAP_2d_2), color = "grey80", size=0.1) +
                geom_point(data = pd[pd$sequencing_batch == run_i,], aes(x = UMAP_2d_1, y = UMAP_2d_2), color = run_color[i], size=0.1) +
                theme_void() +
                theme(legend.position="none") + 
                ggsave(paste0(work_path, batch_id, "_", run_i, ".png"), width = 5, height = 5, dpi = 300))
    }
}



#################################################################################
### Ambient noise (e.g. as might be due to transcript leakage) was assessed by
### examining hemoglobin and collagen transcripts.
#################################################################################

### Supplementary Figure 5

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "/df_cell.rds"))

select_genes = c("Hbb-y", "Hba-x", "Col1a1", "Col2a1", "Hbb-bt", "Hbb-bs")
mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% select_genes,]

pd_sub = pd_all
gene_count = doExtractData(pd_sub, mouse_gene_sub)
x = mouse_gene[rownames(gene_count),]
rownames(gene_count) = as.vector(x$gene_short_name)
saveRDS(gene_count, paste0(work_path, "Npm1_signature.rds"))

df = pd_all
gene_count = readRDS(paste0(work_path, "Npm1_signature.rds"))

for(i in 1:length(select_genes)){
    gene_i = select_genes[i]; print(gene_i)
    
    if(gene_i %in% c("Hbb-y", "Hba-x")){
        df$group = if_else(df$major_trajectory %in% c("Primitive_erythroid"),
                           "Primitive erythroid", "Other celltypes")
        df$group = factor(df$group, levels = c("Primitive erythroid", "Other celltypes"))
    } else if (gene_i %in% c("Hbb-bt", "Hbb-bs")) {
        df$group = if_else(df$major_trajectory %in% c("Definitive_erythroid"),
                           "Definitive erythroid", "Other celltypes")
        df$group = factor(df$group, levels = c("Definitive erythroid", "Other celltypes"))
    } else if (gene_i == "Col1a1") {
        df$group = if_else(df$celltype_update %in% c("Pre-osteoblasts (Sp7+)"),
                           "Pre-osteoblasts", "Other celltypes")
        df$group = factor(df$group, levels = c("Pre-osteoblasts", "Other celltypes"))
    } else {
        df$group = if_else(df$celltype_update %in% c("Early chondrocytes"),
                           "Early chondrocytes", "Other celltypes")
        df$group = factor(df$group, levels = c("Early chondrocytes", "Other celltypes"))
    }
    
    df$exp = exp = as.vector(gene_count[gene_i,])
    print(df %>% group_by(group) %>% summarise(mean_exp = mean(exp)))
    df$exp[exp > 10] = 10
    df_x = df %>% group_by(group, exp) %>% tally() %>%
        left_join(df %>% group_by(group) %>% tally() %>% rename(total_n = n), by = "group") %>%
        mutate(pct = 100*n/total_n)
    
    p = ggplot(df_x, aes(exp, pct, fill = group)) +
        geom_bar(stat="identity", position=position_dodge()) +
        labs(x = paste0("# of reads from ", gene_i), y = "% of cells") +
        theme_classic() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        theme(text=element_text(size=18,  family="Helvetica")) +
        scale_fill_brewer(palette="Set1") +
        scale_x_continuous(breaks=seq(0, 10, 1))
    
    pdf(paste0(work_path, gene_i, "_pct.pdf"), 10, 6)
    print(p)
    dev.off()
    
}









