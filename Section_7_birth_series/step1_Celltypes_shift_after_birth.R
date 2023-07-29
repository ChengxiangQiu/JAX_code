
#################################
### Section - 7, Birth series ###
#################################

################################################################################################
### Re-embedded 2D UMAP of cells from three major cell clusters before E18.75, E18.75, or P0 ###
################################################################################################

### Fig. 7a

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

day_group_color_plate = c("Early" = "#a46cb7",
                          "E18.75" = "#7aa457",
                          "P0" = "#cb6a49",
                          "Other" = "grey90")

for(example_i in c("Hepatocytes", "Adipocytes", "Lung_and_airway")){
    
    example_i = "Renal"; print(example_i)
    
    pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))
    
    day_group = rep("Early", nrow(pd))
    day_group[pd$day == "E18.75"] = "E18.75"
    day_group[pd$day == "P0"] = "P0"
    pd$day_group = as.vector(day_group)
    
    pd$tmp = if_else(pd$day_group == "Early", "Early", "Other")
    p = pd %>%
        ggplot() +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
        geom_point(data = subset(pd, tmp == 'Early'),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
        theme_void() +
        scale_color_manual(values=day_group_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0(work_path, example_i, "_1.png"), width = 6, height = 6, dpi = 300)
    
    pd$tmp = if_else(pd$day_group == "E18.75", "E18.75", "Other")
    p = pd %>%
        ggplot() +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
        geom_point(data = subset(pd, tmp == 'E18.75'),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
        theme_void() +
        scale_color_manual(values=day_group_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0(work_path, example_i, "_2.png"), width = 6, height = 6, dpi = 300)
    
    pd$tmp = if_else(pd$day_group == "P0", "P0", "Other")
    p = pd %>%
        ggplot() +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
        geom_point(data = subset(pd, tmp == 'P0'),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=0.5) +
        theme_void() +
        scale_color_manual(values=day_group_color_plate) +
        theme(legend.position="none") + 
        ggsave(paste0(work_path, example_i, "_3.png"), width = 6, height = 6, dpi = 300)
    
}


#################################################################################################################
### To systematically identify which cell types exhibit abrupt transcriptional changes before vs. after birth ###
#################################################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
x = as.vector(pd_all$day)
x[pd_all$day == "E8.0-E8.5"] = "E8.5"
pd_all$day = as.vector(x)

pd_sub = pd_all[pd_all$day %in% c("E16.0", "E16.25", "E16.5", "E16.75", "E17.0", "E17.25", 
                                  "E17.5", "E17.75", "E18.0", "E18.25", "E18.5", "E18.75", "P0"),]
dat = pd_sub %>% group_by(celltype_update, day) %>% tally() %>% filter(n >= 200)
celltype_1 = dat %>% filter(day == "P0") %>% pull(celltype_update)
celltype_2 = dat %>% filter(day != "P0") %>% group_by(celltype_update) %>% tally() %>% filter(n >= 5) %>% pull(celltype_update)
celltype_x = intersect(celltype_1, celltype_2)

celltype_include = data.frame(celltype_update = as.vector(celltype_x),
                              celltype_name = doSimpleName(as.vector(celltype_x)), 
                              stringsAsFactors = F)

write.table(celltype_include, paste0(work_path, "celltype_include.txt"), row.names=F, col.names=F, sep="\t", quote=F)

### Running embedding on individual cell types
### python Embedding_individual_celltype.py

for(kk in 1:nrow(celltype_include)){
    
    celltype_update_i = celltype_include$celltype_update[kk]
    celltype_name_i = celltype_include$celltype_name[kk]
    print(celltype_update_i)
    
    pd = read.csv(paste0(work_path, celltype_name_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
    
    emb = read.csv(paste0(work_path, celltype_name_i, "_adata_scale.PCs.csv"), header=F, as.is=T)
    colnames(emb) = paste0("PC_", 1:30)
    rownames(emb) = rownames(pd) = as.vector(pd$cell_id)
    emb = as.matrix(emb)
    
    x_table = pd %>% group_by(day) %>% tally() %>% filter(n >= 200)
    x_table_median = median(x_table$n)
    pd_1 = pd %>% filter(day %in% as.vector(x_table$day[x_table$n <= x_table_median])) %>% as.data.frame()
    pd_2 = pd %>% filter(day %in% as.vector(x_table$day[x_table$n > x_table_median])) %>% group_by(day) %>% sample_n(x_table_median) %>% as.data.frame()
    pd_sub = rbind(pd_1, pd_2)
    rownames(pd_sub) = as.vector(pd_sub$cell_id)
    emb_sub = emb[as.vector(pd_sub$cell_id),]
    
    k.param = floor(log2(x_table_median)) + 1 + 1; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
    nn.ranked = Seurat:::NNHelper(
        data = emb_sub,
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
        resultA = cbind(resultA, as.vector(pd_sub$day)[as.vector(nn_matrix[,i])])
    }
    
    resultB = NULL
    for(i in 2:k.param){
        print(i)
        resultB = cbind(resultB, as.vector(resultA[,i] != resultA[,1]))
    }
    
    res = data.frame(day = resultA[,1],
                     pct = apply(resultB, 1, sum)/ncol(resultB))
    
    print(res %>% group_by(day) %>% summarise(mean_pct = mean(pct)) %>% as.data.frame())
    
    saveRDS(res, paste0(work_path, celltype_name_i, "_res.rds"))
    
}


###################################
### Summarizing the kNN results ###
###################################

df = NULL
for(kk in 1:nrow(celltype_include)){
    celltype_update_i = celltype_include$celltype_update[kk]
    celltype_name_i = celltype_include$celltype_name[kk]
    print(celltype_update_i)
    
    res_i = readRDS(paste0(work_path, celltype_name_i, "_res.rds"))
    res_i = res_i %>% group_by(day) %>% summarise(mean_pct = mean(pct)) %>% mutate(celltype_update = celltype_update_i) %>% as.data.frame()
    df = rbind(df, res_i)
}

df_order = df %>% filter(day == "P0") %>% arrange(mean_pct)
df$celltype_update = factor(df$celltype_update, levels = rev(as.vector(df_order$celltype_update)))
df$mean_pct = 100*df$mean_pct

day_list = c("E16.0", "E16.25", "E16.5", "E16.75", "E17.0", "E17.25", 
             "E17.5", "E17.75", "E18.0", "E18.25", "E18.5", "E18.75", "P0")
day_color = c("#5dae46",
              "#855ecd",
              "#b2b044",
              "#c94ca4",
              "#55a574",
              "#d53f63",
              "#4cbad2",
              "#c95534",
              "#617fc6",
              "#d89248",
              "#bd80c4",
              "#7f702f",
              "#bc6476")
names(day_color) = day_list

### Fig. 7b
### size 10 X 6

p = ggplot() +
    geom_point(data = df %>% filter(day == "P0"), aes(x = mean_pct, y = celltype_update), color = "black", size = 3) +
    geom_point(data = df, aes(x = mean_pct, y = celltype_update, color = day), size = 2) +
    scale_color_manual(values=day_color) +
    labs(x = "Mean % of the nearest neighboring cells from different timepoints", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))







