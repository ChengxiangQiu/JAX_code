
###############################
##### Update celltype names ###
###############################

source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"

file_list = paste0("E", c(95,105,115,125,135,145,155,165))

celltype_list = NULL
for(i in file_list){
    pd_i = read.csv(paste0(work_path, "/mosta/sc_data/", i, ".df_cell.csv"), header=T, row.names=1, as.is=T)
    celltype_list = c(celltype_list, as.vector(pd_i$celltype))
}
celltype_list = as.vector(names(table(celltype_list)))

pd = readRDS(paste0(work_path, "/LPM/LPM_adata_scale.obs.rds"))
celltype_list_update = as.vector(names(table(pd$celltype_sub_clustering)))

celltype_list = c("Airway smooth muscle cells" = "Airway smooth muscle cells",
                  "Allantois" = "Allantois",
                  "Amniotic mesoderm" = "Amniotic mesoderm",
                  "Cardiopharyngeal mesoderm (Tbx1+)" = "Cardiopharyngeal mesoderm (Tbx1+)",
                  "Extraembryonic mesoderm" = "Extraembryonic mesoderm",
                  "Foregut mesenchyme" = "Foregut mesenchyme",
                  "Gastrointestinal smooth muscle cells" = "Gastrointestinal smooth muscle cells",
                  "Gonad progenitor cells" = "Gonad progenitor cells",
                  "Gut mesenchyme" = "Gut mesenchyme",
                  "Hepatic mesenchyme" = "Hepatic mesenchyme",
                  "Intermediate mesoderm" = "Renal capsule",
                  "LPM-derived cells (Cfh+)" = "Renal stromal cells",
                  "Lung mesenchyme" = "Lung mesenchyme",
                  "Mesothelial cells" = "Mesothelial cells",
                  "Migrated vascular smooth muscle cells" = "Meninges",
                  "Proepicardium (Tbx18+)" = "Proepicardium (Tbx18+)",
                  "Somatic mesoderm" = "Somatic mesoderm",
                  "Splanchnic mesoderm" = "Splanchnic mesoderm",
                  "Vascular smooth muscle cells" = "Vascular smooth muscle cells",
                  "Vascular smooth muscle cells (Pparg+)" = "Vascular smooth muscle cells (Pparg+)")

celltype_list = data.frame(old_name = names(celltype_list),
                           new_name = as.vector(celltype_list), stringsAsFactors = F)

write_name = gsub("[(|)|+]", "", as.vector(celltype_list$new_name))
write_name = gsub(" ", "_", write_name)

celltype_list$write_name = write_name

col_name = gsub("[+]", "", as.vector(celltype_list$old_name))
col_name = gsub("[(|-]", ".", col_name)
col_name = gsub("[)]", "..", col_name)
col_name = gsub(" ", ".", col_name)

celltype_list$col_name = as.vector(col_name)

saveRDS(celltype_list, paste0(work_path, "/mosta/celltype_update.rds"))



############################################
### plot every cell type X every section ###
############################################

source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
save_path = "/net/shendure/vol10/www/content/members/cxqiu/public/nobackup/jax/spatial_mapping"

file_list = gsub(".MOSTA.h5ad", "", as.vector(read.table(paste0(work_path, "/mosta/file_list.txt"))$V1))

file_list_x = c("E9.5_E1S1", "E9.5_E2S1", "E9.5_E2S2", "E9.5_E2S3", "E9.5_E2S4", "E10.5_E2S1")

celltype_list = readRDS(paste0(work_path, "/mosta/celltype_update.rds"))

### Of note, after performing Tangram, a matrix with mapping prob between pairwise cell and voxel is returned
### for individual section. In this matrix, the sum of mapping prob across voxels for individual cell is 1.
### Then, we simply aggregated mapping prob for cells within each cell type (i.e. subpopulation of LPM), which will
### used as mapping prob (across voxels) for each given cell type (XXX.result.csv, row is voxel, col is cell type)

### In the below code, for individual section X individual cell type, we smooth the mapping prob for each voxel
### by average its k-nearest neighboring voxels (k = log2(total voxels)), and then rescale it to 0->1.

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

cnt = 0
for(i in 1:length(file_list)){
    cnt = cnt + 1
    file_i = file_list[i]; print(paste0(i, "/", file_i))
    
    if(file_i %in% file_list_x){
        t1 = 1.8; t2 = 1.6
    } else {
        t1 = 0.8; t2 = 0.6
    }
    
    pd = read.csv(paste0(work_path, "/mosta/result/", file_i, ".result.coor.csv"), header=F, as.is=T)
    colnames(pd) = c("coor_1", "coor_2")
    pd$coor_2 = 0 - pd$coor_2
    
    anno = read.csv(paste0(work_path, "/mosta/result/", file_i, ".result.csv"), header=T, as.is=T)
    anno_list = colnames(anno)
    anno_list[1] = gsub("X..", "", anno_list[1])
    
    kk = round(log(nrow(pd)))
    k.param = kk + 1; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
    nn.ranked = Seurat:::NNHelper(
        data = pd,
        k = k.param,
        method = nn.method,
        searchtype = "standard",
        eps = nn.eps,
        metric = annoy.metric)
    
    nn.ranked = Indices(object = nn.ranked)
    nn_matrix = nn.ranked
    
    for(j in 1:length(anno_list)){
        anno_j = anno_list[j]
        name = celltype_list$write_name[celltype_list$col_name == anno_j]
        print(paste0(cnt, "/", file_i, ":", name))
        
        df = cbind(pd, as.vector(anno[,j]))
        names(df) = c("coor_1", "coor_2", "mapping_prob")
        
        resA = NULL
        for(i in 1:ncol(nn_matrix)){
            resA = cbind(resA, df$mapping_prob[as.vector(nn_matrix[,i])])
        }
        
        df$mapping_prob_smooth = apply(resA, 1, mean)
        
        x_max = max(df$mapping_prob_smooth); x_min = min(df$mapping_prob_smooth)
        df$mapping_prob_smooth_norm = (df$mapping_prob_smooth - x_min)/(x_max - x_min)
        
        try(print(df %>%
            ggplot() +
            geom_point(aes(x = coor_1, y = coor_2, color = mapping_prob_smooth_norm), size=t2) +
            sc +
            theme_void() +
            theme(legend.position="none") +
                ggsave(paste0(paste0(save_path, "/", file_i, "/", file_i, "_", name, ".png")),
                       dpi = 300,
                       height  = 8, 
                       width = 6)), silent = TRUE)
    }

}
    
    
    





