
######################################
### Section - 3, Kidney_mesenchyme ###
######################################

#######################################################
### Spatial mapping using Tangram + Mosta (Fig. 3e) ###
#######################################################

### The Mosta (spatial transcriptome) data was downloaded from https://db.cngb.org/stomics/mosta/download/
### Stereo-seq data from each mouse embryo in MOSAT database, including bin 50 and segmented single-cell.

###############################################
### Plotting the backbone for each sections ###
###############################################

### This step is not necessary unless you want to re-plot the spatial annotations
### Make sure you have extracted the spatial coordinates using the first part of Spatial_mapping.py

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

file_list = gsub(".MOSTA.h5ad", "", as.vector(read.table(paste0(work_path, "Mosta_file_list.txt"))$V1))

dat = NULL

for(i in 1:length(file_list)){
    file_i = file_list[i]
    x = read.csv(paste0(work_path, "annotation/", file_i, ".spatial_color.csv"), header=F)
    y = read.csv(paste0(work_path, "annotation/", file_i, ".spatial_color_id.csv"), header=F)
    dat = rbind(dat, data.frame(color = as.vector(x$V1),
                                color_id = as.vector(y$V1)))
}
dat = unique(dat)
dat$color = as.vector(dat$color)
dat$color_id = as.vector(dat$color_id)
dat$color = substr(dat$color,1,nchar(dat$color)-2)

mosta_color_plate = as.vector(dat$color)
names(mosta_color_plate) = as.vector(dat$color_id)

df = data.frame(anno_id = names(mosta_color_plate),
                x = sample(1:4, length(mosta_color_plate), replace = T),
                y = sample(1:4, length(mosta_color_plate), replace = T))

p = df %>% 
    ggplot(aes(x, y, color = anno_id)) + geom_point(size = 5) + 
    theme_classic(base_size = 10) +
    scale_color_manual(values=mosta_color_plate)

mosta_color_plate = c(mosta_color_plate, "other" = "grey90")

file_list_x = c("E9.5_E1S1", "E9.5_E2S1", "E9.5_E2S2", "E9.5_E2S3", "E9.5_E2S4", "E10.5_E2S1")

for(i in 1:length(file_list)){
    
    if(file_i %in% file_list_x){
        t1 = 1.8; t2 = 1.6
    } else {
        t1 = 0.8; t2 = 0.6
    }
    
    file_i = file_list[i]; print(paste0(i, "/", file_i))
    pd = read.csv(paste0(work_path, "annotation/", file_i, ".spatial_coor.csv"), header=F, as.is=T)
    colnames(pd) = c("coor_1", "coor_2")
    pd$annotation = as.vector(read.csv(paste0(work_path, "annotation/", file_i, ".spatial_anno.csv"), header=F, as.is=T)$V1)
    
    pd$coor_2 = 0 - pd$coor_2
    
    mosta_color_plate_sub = mosta_color_plate[names(mosta_color_plate) %in% pd$annotation]
    p = pd %>%
        ggplot() +
        geom_point(aes(x = coor_1, y = coor_2), size=t1) +
        geom_point(aes(x = coor_1, y = coor_2, color = annotation), size=t2) +
        scale_color_manual(values=mosta_color_plate_sub) +
        theme_void()
    
    pdf(paste0(work_path, "annotation/", file_i, ".spatial_images.pdf"))
    print(p)
    dev.off()
    
    p = pd %>%
        ggplot() +
        geom_point(aes(x = coor_1, y = coor_2), size=t1) +
        geom_point(aes(x = coor_1, y = coor_2, color = annotation), size=t2) +
        scale_color_manual(values=mosta_color_plate) +
        theme_void() + 
        theme(legend.position="none")
    
    try(p + ggsave(paste0(work_path, "annotation/", file_i, ".spatial_images.png"),
                   dpi = 300,
                   height  = 8, 
                   width = 6), silent = TRUE)
    
    anno_list = names(table(pd$annotation))
    for(ii in anno_list){
        pd$tmp = if_else(pd$annotation == ii, ii, "other")
        output_name = gsub(" ", "_", ii)
        
        p = pd %>%
            ggplot() +
            geom_point(aes(x = coor_1, y = coor_2), size=t1) +
            geom_point(aes(x = coor_1, y = coor_2, color = tmp), size=t2) +
            scale_color_manual(values=mosta_color_plate) +
            theme_void()
        
        try(p + theme(legend.position="none") +
                ggsave(paste0(work_path, "annotation/", file_i, "_", output_name, ".spatial_images.png"),
                       dpi = 300,
                       height  = 8, 
                       width = 6), silent = TRUE)
    }
    
}


############################################################################################
### Preparing the subset of sc-RNA-seq data from LPM to perform spatial mapping analysis ###
############################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

dir.create(paste0(work_path, "sc_data"))

pd = readRDS(paste0(work_path, "LPM_adata_scale.obs.rds"))

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
    x_table = table(pd_sub$lateral_plate_mesoderm_sub_clustering)
    x_1 = names(x_table)[x_table >= 200 & x_table <= 5000]
    x_2 = names(x_table)[x_table > 5000]
    
    pd_sub_1 = pd_sub[pd_sub$lateral_plate_mesoderm_sub_clustering %in% x_1,]
    pd_sub_2 = pd_sub[pd_sub$lateral_plate_mesoderm_sub_clustering %in% x_2,] %>% 
        group_by(lateral_plate_mesoderm_sub_clustering) %>%
        sample_n(5000)
    pd_sub_2 = pd_sub[pd_sub$cell_id %in% pd_sub_2$cell_id,]
    pd_sub = rbind(pd_sub_1, pd_sub_2)
    pd_sub = pd_sub[,c("cell_id", "day", "lateral_plate_mesoderm_sub_clustering")]
    rownames(pd_sub) = as.vector(pd_sub$cell_id)
    
    fd_sub = mouse_gene[mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA') & 
        mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
    
    gene_count = doExtractData(pd_sub, fd_sub)

    df_gene = mouse_gene[rownames(gene_count),]
    df_gene$mean_exp = Matrix::rowMeans(gene_count)
    df_gene_x = df_gene %>% group_by(gene_short_name) %>% slice_max(order_by = mean_exp, n = 1, with_ties = F)
    gene_count = gene_count[rownames(gene_count) %in% df_gene_x$gene_ID,]
    df_gene = mouse_gene[rownames(gene_count),]
    row.names(gene_count) = rownames(df_gene) = as.vector(df_gene$gene_short_name)
    
    Matrix::writeMM(t(gene_count), paste0(work_path, 'sc_data/', day_i, ".gene_count.mtx"))
    write.csv(df_gene, paste0(work_path, 'sc_data/', day_i, ".df_gene.csv"))
    write.csv(pd_sub, paste0(work_path, 'sc_data/', day_i, ".df_cell.csv"))
}



############################################
### Plotting the spatial mapping results ###
############################################

### After running Part-2 and Part-3 of Spatial_mapping.py

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

file_list = paste0("E", c(95,105,115,125,135,145,155,165))

celltype_list = NULL
for(i in file_list){
    pd_i = read.csv(paste0(work_path, "sc_data/", i, ".df_cell.csv"), header=T, row.names=1, as.is=T)
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





### Of note, after performing Tangram, a matrix with mapping prob between pairwise cell and voxel is returned
### for individual section. In this matrix, the sum of mapping prob across voxels for individual cell is 1.
### Then, we simply aggregated mapping prob for cells within each cell type (i.e. subpopulation of LPM), which will
### used as mapping prob (across voxels) for each given cell type (XXX.result.csv, row is voxel, col is cell type)

### In the below script, for individual section X individual cell type, we smooth the mapping prob for each voxel
### by average its k-nearest neighboring voxels (k = log2(total voxels)), and then rescale it to 0->1.

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"


pd = readRDS(paste0(work_path, "LPM_adata_scale.obs.rds"))
celltype_list = data.frame(old_name = names(table(pd$lateral_plate_mesoderm_sub_clustering)),
                           new_name = names(table(pd$lateral_plate_mesoderm_sub_clustering)), stringsAsFactors = F)

write_name = gsub("[(|)|+]", "", as.vector(celltype_list$new_name))
write_name = gsub(" ", "_", write_name)

celltype_list$write_name = write_name

col_name = gsub("[+]", "", as.vector(celltype_list$old_name))
col_name = gsub("[(|-]", ".", col_name)
col_name = gsub("[)]", "..", col_name)
col_name = gsub(" ", ".", col_name)

celltype_list$col_name = as.vector(col_name)

file_list = gsub(".MOSTA.h5ad", "", as.vector(read.table(paste0(work_path, "MOSTA_file_list.txt"))$V1))
file_list_x = c("E9.5_E1S1", "E9.5_E2S1", "E9.5_E2S2", "E9.5_E2S3", "E9.5_E2S4", "E10.5_E2S1")

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
    
    pd = read.csv(paste0(work_path, "result/", file_i, ".result.coor.csv"), header=F, as.is=T)
    colnames(pd) = c("coor_1", "coor_2")
    pd$coor_2 = 0 - pd$coor_2
    
    anno = read.csv(paste0(work_path, "result/", file_i, ".result.csv"), header=T, as.is=T)
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
                      ggsave(paste0(paste0(work_path, "result/", file_i, "/", file_i, "_", name, ".png")),
                             dpi = 300,
                             height  = 8, 
                             width = 6)), silent = TRUE)
    }
    
}











