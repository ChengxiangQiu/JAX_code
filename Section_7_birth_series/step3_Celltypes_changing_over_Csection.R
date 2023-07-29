
#################################
### Section - 7, Birth series ###
#################################

###########################################################
### which cell type is changing across C-section series ###
###########################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "pd_birth.rds"))
pd$anno = as.vector(pd$major_trajectory)
emb = readRDS(paste0(work_path, "Birth_series.PCs.rds"))
emb = as.matrix(emb)

index = pd$day %in% paste0("Csection_", c(0,20,40,60,80), "m")
pd_sub = pd[index,]
emb_sub = emb[index,]

k.param = 11; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = emb_sub,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked

saveRDS(nn_matrix, paste0(work_path, "nn_matrix.rds"))

day_value = gsub("Csection_", "", as.vector(pd_sub$day))
day_value = gsub("m", "", day_value) 
day_value = as.numeric(day_value)

nn_res = NULL
for(i in 1:ncol(nn_matrix)){
    nn_res = cbind(nn_res, day_value[as.vector(nn_matrix[,i])])
}
dat = data.frame(org_day = nn_res[,1], 
                 nn_day = apply(nn_res[,2:ncol(nn_res)], 1, mean),
                 anno = as.vector(pd_sub$anno))

celltype_list = names(table(dat$anno))
res = NULL
for(celltype_i in celltype_list){
    res = rbind(res,
                data.frame(anno = celltype_i,
                           cor = cor.test(dat$org_day[dat$anno == celltype_i], dat$nn_day[dat$anno == celltype_i])$estimate, stringsAsFactors = FALSE))
}
res = res[order(res$cor),]
res$anno = factor(res$anno, levels = as.vector(res$anno))

p = res %>%
    ggplot(aes(anno, cor, fill = anno)) + 
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_manual(values = major_trajectory_color_plate) +
    labs(x = "", y = "") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
