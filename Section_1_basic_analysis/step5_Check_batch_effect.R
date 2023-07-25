
###################################
### Section - 1, Basic analysis ###
###################################

################################################################
### checking potential batch effects (Supplementary Fig. 3b) ###
################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")

pd = readRDS("df_cell.rds")
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