
source("~/work/scripts/tome/utils.R")
setwd("/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq/mosta")

file_list = gsub(".MOSTA.h5ad", "", as.vector(read.table("file_list.txt")$V1))

# dat = NULL
# 
# for(i in 1:length(file_list)){
#     file_i = file_list[i]
#     x = read.csv(paste0("./annotation/", file_i, ".spatial_color.csv"), header=F)
#     y = read.csv(paste0("./annotation/", file_i, ".spatial_color_id.csv"), header=F)
#     dat = rbind(dat, data.frame(color = as.vector(x$V1),
#                                 color_id = as.vector(y$V1)))
# }
# dat = unique(dat)
# dat$color = as.vector(dat$color)
# dat$color_id = as.vector(dat$color_id)
# dat$color = substr(dat$color,1,nchar(dat$color)-2)

# mosta_color_plate = as.vector(dat$color)
# names(mosta_color_plate) = as.vector(dat$color_id)

# saveRDS(mosta_color_plate, "mosta_color_plate.rds")

mosta_color_plate = readRDS("mosta_color_plate.rds")

df = data.frame(anno_id = names(mosta_color_plate),
                x = sample(1:4, length(mosta_color_plate), replace = T),
                y = sample(1:4, length(mosta_color_plate), replace = T))

p = df %>% 
    ggplot(aes(x, y, color = anno_id)) + geom_point(size = 5) + 
    theme_classic(base_size = 10) +
    scale_color_manual(values=mosta_color_plate)

pdf(paste0("mosta_color_plate.pdf"))
print(p)
dev.off()

mosta_color_plate = c(mosta_color_plate, "other" = "grey90")

file_list_x = c("E9.5_E1S1", "E9.5_E2S1", "E9.5_E2S2", "E9.5_E2S3", "E9.5_E2S4", "E10.5_E2S1")

for(i in 1:length(file_list)){
    
    if(file_i %in% file_list_x){
        t1 = 1.8; t2 = 1.6
    } else {
        t1 = 0.8; t2 = 0.6
    }
    
    file_i = file_list[i]; print(paste0(i, "/", file_i))
    pd = read.csv(paste0("./annotation/", file_i, ".spatial_coor.csv"), header=F, as.is=T)
    colnames(pd) = c("coor_1", "coor_2")
    pd$annotation = as.vector(read.csv(paste0("./annotation/", file_i, ".spatial_anno.csv"), header=F, as.is=T)$V1)
    
    pd$coor_2 = 0 - pd$coor_2
    
    mosta_color_plate_sub = mosta_color_plate[names(mosta_color_plate) %in% pd$annotation]
    p = pd %>%
        ggplot() +
        geom_point(aes(x = coor_1, y = coor_2), size=t1) +
        geom_point(aes(x = coor_1, y = coor_2, color = annotation), size=t2) +
        scale_color_manual(values=mosta_color_plate_sub) +
        theme_void()
    
    pdf(paste0("./annotation/", file_i, ".spatial_images.pdf"))
    print(p)
    dev.off()
    
    p = pd %>%
        ggplot() +
        geom_point(aes(x = coor_1, y = coor_2), size=t1) +
        geom_point(aes(x = coor_1, y = coor_2, color = annotation), size=t2) +
        scale_color_manual(values=mosta_color_plate) +
        theme_void() + 
        theme(legend.position="none")
  
    try(p +
            ggsave(paste0("./annotation/", file_i, ".spatial_images.png"),
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
                ggsave(paste0("./annotation/", file_i, "_", output_name, ".spatial_images.png"),
                       dpi = 300,
                       height  = 8, 
                       width = 6), silent = TRUE)
    }

}


