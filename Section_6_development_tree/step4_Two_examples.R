
#####################################
### Section - 6, Development tree ###
#####################################

#########################################################################################################
### We present two examples in Fig. 5 a-e for illustration of basis for edge inference heuristic

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell_graph.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")

####################################
### example_1: suppressor_cells ####
####################################

celltype_include = c("Hematopoietic stem cells (Cd34+)",
                     "Hematopoietic stem cells (Mpo+)",
                     "Monocytic myeloid-derived suppressor cells",
                     "PMN myeloid-derived suppressor cells")
example_i = "suppressor_cells"

pd_sub = pd_all[pd_all$celltype_new %in% celltype_include,]
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count = doExtractData(pd_sub, mouse_gene_sub)

writeMM(t(gene_count), paste0(work_path, example_i, ".gene_count.mtx"))
write.csv(pd, paste0(work_path, example_i, ".df_cell.csv"))

### python Two_examples.py
pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)

nn_matrix = read.csv(paste0(work_path, example_i, "_adata_scale.kNN_15.csv"), as.is=T, header=F)
nn_matrix = as.matrix(nn_matrix)
nn_matrix = nn_matrix + 1 ### python and R using different start index

x = data.frame(i = rep(1:nrow(nn_matrix), ncol(nn_matrix)),
               j = c(nn_matrix), stringsAsFactors = FALSE)

dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                           j = as.numeric(as.vector(x$j)),
                           x = 1)

dat_t = t(dat) + dat
x = data.frame(summary(dat_t))
x = x[x$x == 2 & x$i > x$j,]
x$x = NULL   ### x saves the MNN pairs

y = data.frame(i = 1:nrow(pd_x),
               j = 1:nrow(pd_x),
               meta_group = as.vector(pd_x$meta_group), stringsAsFactors = FALSE)

dat = x %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j") %>%
    group_by(meta_group.x, meta_group.y) %>% tally() 

library(reshape2)
obs = dcast(dat, meta_group.x~meta_group.y)
rownames(obs) = as.vector(obs[,1])
obs = obs[,-1]
obs[is.na(obs)] = 0
obs = as.matrix(obs)

diag(obs) = 0

obs_x = obs + t(obs)
obs_y = as.vector(obs_x[upper.tri(obs_x)])

group = NULL
for(i in 2:nrow(obs_x)){
    print(i)
    for(j in 1:(i-1)){
        group = rbind(group, data.frame(x = colnames(obs_x)[j],
                                        y = rownames(obs_x)[i], stringsAsFactors = F))
    }
}

group$edge_num = obs_y
nodes = pd_x %>% group_by(meta_group, celltype_new) %>% tally() %>% rename(celltype_num = n)
group = group %>% 
    left_join(nodes %>% rename(x = meta_group) %>% select(x, celltype_new, celltype_num), by = "x") %>% rename(x_name = celltype_new, x_size = celltype_num) %>%
    left_join(nodes %>% rename(y = meta_group) %>% select(y, celltype_new, celltype_num), by = "y") %>% rename(y_name = celltype_new, y_size = celltype_num)
group$min_size = if_else(group$x_size < group$y_size, group$x_size, group$y_size)
group$edge_num_norm = group$edge_num/log2(15*group$min_size)
group$min_size = NULL

x_rev = x
x_rev$i = x$j
x_rev$j = x$i
x_rev = rbind(x, x_rev)
dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j") %>%
    filter(meta_group.x != meta_group.y)
dat_1 = dat[dat$meta_group.x == "B_M11" & dat$meta_group.y == "B_M12",]
dat_2 = dat[dat$meta_group.x == "B_M12" & dat$meta_group.y == "B_M19",]
dat_3 = dat[dat$meta_group.x == "B_M19" & dat$meta_group.y == "B_M23",]
dat = rbind(dat_1, dat_2, dat_3)

group = rep("other", nrow(pd_x))
group[c(1:nrow(pd_x)) %in% as.vector(dat_1$i)] = "group_1"
group[c(1:nrow(pd_x)) %in% as.vector(dat_1$j)] = "group_2"
group[c(1:nrow(pd_x)) %in% as.vector(dat_2$i)] = "group_3"
group[c(1:nrow(pd_x)) %in% as.vector(dat_2$j)] = "group_4"
group[c(1:nrow(pd_x)) %in% as.vector(dat_3$i)] = "group_5"
group[c(1:nrow(pd_x)) %in% as.vector(dat_3$j)] = "group_6"
pd_x$group = as.vector(group)

color_plate = c("#b69340", "#8d70c9", "#67a64e", "#c8588c", "#49adad", "#cb5b42")
names(color_plate) = paste0("group_", 1:6)

p = ggplot() +
    geom_point(data = pd_x, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6, color = "grey80") +
    geom_point(data = pd_x[pd_x$group != "other",], aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(data = pd_x[pd_x$group != "other",], aes(x = UMAP_2d_1, y = UMAP_2d_2, color = group), size=0.8) +
    theme_void() +
    scale_color_manual(values=color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".MNN.png"), width = 6, height = 6, dpi = 300)

pd_all_day = pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n)

df = pd_x %>% filter(group != "other") %>%
    group_by(group, day) %>% tally() %>% left_join(pd_all_day) %>% mutate(frac = 100*n/total_n)
df$day = factor(df$day, levels = day_list[day_list %in% df$day])

p = df %>% 
    ggplot(aes(x=day, y=frac, fill = group)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(group)) + 
    labs(x='',y='% of cells') +
    scale_fill_manual(values=color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0(work_path, example_i, ".MNN.pdf"), 6, 9)
p
dev.off()




#######################
### example_2: lung ###
#######################

celltype_include = c("Gut",
                     "Lung progenitor cells",
                     "Alveolar Type 2 cells")
example_i = "lung"

pd_sub = pd_all[pd_all$celltype_new %in% celltype_include,]
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count = doExtractData(pd_sub, mouse_gene_sub)

writeMM(t(gene_count), paste0(work_path, example_i, ".gene_count.mtx"))
write.csv(pd, paste0(work_path, example_i, ".df_cell.csv"))

### python Two_examples.py
pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)

nn_matrix = read.csv(paste0(work_path, example_i, "_adata_scale.kNN_15.csv"), as.is=T, header=F)
nn_matrix = as.matrix(nn_matrix)
nn_matrix = nn_matrix + 1 ### python and R using different start index

x = data.frame(i = rep(1:nrow(nn_matrix), ncol(nn_matrix)),
               j = c(nn_matrix), stringsAsFactors = FALSE)

dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                           j = as.numeric(as.vector(x$j)),
                           x = 1)

dat_t = t(dat) + dat
x = data.frame(summary(dat_t))
x = x[x$x == 2 & x$i > x$j,]
x$x = NULL   ### x saves the MNN pairs

y = data.frame(i = 1:nrow(pd_x),
               j = 1:nrow(pd_x),
               meta_group = as.vector(pd_x$meta_group), stringsAsFactors = FALSE)

dat = x %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j") %>%
    group_by(meta_group.x, meta_group.y) %>% tally() 

library(reshape2)
obs = dcast(dat, meta_group.x~meta_group.y)
rownames(obs) = as.vector(obs[,1])
obs = obs[,-1]
obs[is.na(obs)] = 0
obs = as.matrix(obs)

diag(obs) = 0

obs_x = obs + t(obs)
obs_y = as.vector(obs_x[upper.tri(obs_x)])

group = NULL
for(i in 2:nrow(obs_x)){
    print(i)
    for(j in 1:(i-1)){
        group = rbind(group, data.frame(x = colnames(obs_x)[j],
                                        y = rownames(obs_x)[i], stringsAsFactors = F))
    }
}

group$edge_num = obs_y
nodes = pd_x %>% group_by(meta_group, celltype_new) %>% tally() %>% rename(celltype_num = n)
group = group %>% 
    left_join(nodes %>% rename(x = meta_group) %>% select(x, celltype_new, celltype_num), by = "x") %>% rename(x_name = celltype_new, x_size = celltype_num) %>%
    left_join(nodes %>% rename(y = meta_group) %>% select(y, celltype_new, celltype_num), by = "y") %>% rename(y_name = celltype_new, y_size = celltype_num)
group$min_size = if_else(group$x_size < group$y_size, group$x_size, group$y_size)
group$edge_num_norm = group$edge_num/log2(15*group$min_size)
group$min_size = NULL


x_rev = x
x_rev$i = x$j
x_rev$j = x$i
x_rev = rbind(x, x_rev)
dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j") %>%
    filter(meta_group.x != meta_group.y)
dat_1 = dat[dat$meta_group.x == "G_M8" & dat$meta_group.y == "G_M13",]
dat_2 = dat[dat$meta_group.x == "G_M13" & dat$meta_group.y == "G_M4",]
dat = rbind(dat_1, dat_2)


group = rep("other", nrow(pd_x))
group[c(1:nrow(pd_x)) %in% as.vector(dat_1$i)] = "group_1"
group[c(1:nrow(pd_x)) %in% as.vector(dat_1$j)] = "group_2"
group[c(1:nrow(pd_x)) %in% as.vector(dat_2$i)] = "group_3"
group[c(1:nrow(pd_x)) %in% as.vector(dat_2$j)] = "group_4"
pd_x$group = as.vector(group)

color_plate = c("#b69340", "#8d70c9", "#67a64e", "#c8588c")
names(color_plate) = paste0("group_", 1:4)

p = ggplot() +
    geom_point(data = pd_x, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6, color = "grey80") +
    geom_point(data = pd_x[pd_x$group != "other",], aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
    geom_point(data = pd_x[pd_x$group != "other",], aes(x = UMAP_2d_1, y = UMAP_2d_2, color = group), size=0.8) +
    theme_void() +
    scale_color_manual(values=color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".MNN.png"), width = 6, height = 6, dpi = 300)

pd_all_day = pd_all %>% group_by(day) %>% tally() %>% rename(total_n = n)

df = pd_x %>% filter(group != "other") %>%
    group_by(group, day) %>% tally() %>% left_join(pd_all_day) %>% mutate(frac = 100*n/total_n)
df$day = factor(df$day, levels = day_list[day_list %in% df$day])

p = df %>% 
    ggplot(aes(x=day, y=frac, fill = group)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(group)) + 
    labs(x='',y='% of cells') +
    scale_fill_manual(values=color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))

pdf(paste0(work_path, example_i, ".MNN.pdf"), 6, 7)
p
dev.off()



