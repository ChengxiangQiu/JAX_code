
#####################################
### Section - 6, Development tree ###
#####################################

#########################################################################################################
### The MNN approach used for graph construction is robust to subsampling and choice of the k parameter.

################################################################################
### First, we examined whether the MNNs that we identified between different cell types 
### were enriched for cells from the same embryo.

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
pd_all_graph = readRDS(work_path, "df_cell_graph.rds")
pd_embryo = pd_all %>% group_by(embryo_id, day) %>% tally()

early_day = c("E8.5","E8.75","E9.0","E9.25","E9.5","E9.75","E10.0")
select_day = c("E13.0", "E13.25", "E13.5", "E13.75")

embryo_early = as.vector(pd_embryo$embryo_id[pd_embryo$day %in% early_day])
embryo_select = as.vector(pd_embryo$embryo_id[pd_embryo$day %in% select_day])

system_list = c("Endothelium",                    
                "Epithelial_cells",               
                "Eye",                            
                "Gut",                            
                "Notochord",                      
                "PNS_glia",                       
                "PNS_neurons",                    
                "Renal",                          
                "Lateral_plate_mesoderm",         
                "Blood",     
                "Neuroectoderm",
                "Brain_spinal_cord",              
                "Mesoderm")      

res_1 = res_2 = res_3 = NULL

for(kk in 1:length(system_list)){
    system_i = system_list[kk]
    print(system_i)
    
    pd = read.csv(paste0(work_path, system_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
    x = readRDS(paste0(work_path, system_i, ".MNN.rds"))
    
    pd$meta_group = NULL
    pd = pd %>% left_join(pd_all_graph[,c("cell_id","meta_group")])
    
    y = data.frame(meta_group_i = as.vector(pd$meta_group)[as.vector(x$i)],
                   meta_group_j = as.vector(pd$meta_group)[as.vector(x$j)],
                   embryo_id_i = as.vector(pd$embryo_id)[as.vector(x$i)],
                   embryo_id_j = as.vector(pd$embryo_id)[as.vector(x$j)], stringsAsFactors = F)
    
    y_sub = y[y$meta_group_i != y$meta_group_j,]
    res_1 = rbind(res_1, data.frame(Diff = sum(y_sub$embryo_id_i != y_sub$embryo_id_j),
                                    Same = sum(y_sub$embryo_id_i == y_sub$embryo_id_j))) 
    
    y_sub_2 = y_sub[(y_sub$embryo_id_i %in% embryo_early) |
                        (y_sub$embryo_id_j %in% embryo_early),]
    res_2 = rbind(res_2, data.frame(Diff = sum(y_sub_2$embryo_id_i != y_sub_2$embryo_id_j),
                                    Same = sum(y_sub_2$embryo_id_i == y_sub_2$embryo_id_j))) 
    
    y_sub_3 = y_sub[(y_sub$embryo_id_i %in% embryo_select) |
                        (y_sub$embryo_id_j %in% embryo_select),]
    res_3 = rbind(res_3, data.frame(Diff = sum(y_sub_3$embryo_id_i != y_sub_3$embryo_id_j),
                                    Same = sum(y_sub_3$embryo_id_i == y_sub_3$embryo_id_j))) 
    
}

saveRDS(list(res_1, res_2, res_3), paste0(work_path, "MNN_same_diff_embryo.rds"))

res = data.frame(res_1, stringsAsFactors = F)
colnames(res) = c("Diff", "Same")
res$system = system_list
res = res[(res$Diff + res$Same) != 0,]
res$Diff_pct = res$Diff/(res$Diff + res$Same)
res$Same_pct = res$Same/(res$Diff + res$Same)
print(sum(res$Same)/sum(res$Same + res$Diff))

df = data.frame(system = c(res$system, res$system),
                pct = c(res$Diff_pct, res$Same_pct),
                group = rep(c("Diff_embryo", "Same_embryo"), each = nrow(res)))
df$system = factor(df$system, levels = system_list)

p = ggplot(df, aes(x = system, y = pct*100, fill = group)) +
    geom_bar(stat="identity", width = 0.7) +
    labs(x = "system", y = "percent", fill = "group") +
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black")) 

res = data.frame(res_2, stringsAsFactors = F)
colnames(res) = c("Diff", "Same")
res$system = system_list
res = res[(res$Diff + res$Same) != 0,]
res$Diff_pct = res$Diff/(res$Diff + res$Same)
res$Same_pct = res$Same/(res$Diff + res$Same)
print(sum(res$Same)/sum(res$Same + res$Diff))

df = data.frame(system = c(res$system, res$system),
                pct = c(res$Diff_pct, res$Same_pct),
                group = rep(c("Diff_embryo", "Same_embryo"), each = nrow(res)))
df$system = factor(df$system, levels = system_list)

p = ggplot(df, aes(x = system, y = pct*100, fill = group)) +
    geom_bar(stat="identity", width = 0.7) +
    labs(x = "system", y = "percent", fill = "group") +
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black")) 

res = data.frame(res_3, stringsAsFactors = F)
colnames(res) = c("Diff", "Same")
res$system = system_list
res = res[(res$Diff + res$Same) != 0,]
res$Diff_pct = res$Diff/(res$Diff + res$Same)
res$Same_pct = res$Same/(res$Diff + res$Same)
print(sum(res$Same)/sum(res$Same + res$Diff))

df = data.frame(system = c(res$system, res$system),
                pct = c(res$Diff_pct, res$Same_pct),
                group = rep(c("Diff_embryo", "Same_embryo"), each = nrow(res)))
df$system = factor(df$system, levels = system_list)

p = ggplot(df, aes(x = system, y = pct*100, fill = group)) +
    geom_bar(stat="identity", width = 0.7) +
    labs(x = "system", y = "percent", fill = "group") +
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black")) 


################################################################################
### Second, to assess the robustness of MNNs to cell sampling, we randomly subsampled 80% 
### of cells from each developmental system during organogenesis & fetal development

### python Graph_robust.py

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(work_path, "df_cell_graph.rds")

system_list = c("Endothelium",                    
                "Epithelial_cells",               
                "Eye",                            
                "Gut",                            
                "Notochord",                      
                "PNS_glia",                       
                "PNS_neurons",                    
                "Renal",                          
                "Lateral_plate_mesoderm",         
                "Blood",     
                "Neuroectoderm",
                "Brain_spinal_cord",              
                "Mesoderm")    

for(kk in 1:length(system_list)){
    
    system_i = system_list[kk]
    print(system_i)
    
    pd = read.csv(paste0(work_path, system_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
    rownames(pd) = as.vector(pd$cell_id)
    pd$celltype_new = pd$system = pd$meta_group = NULL
    pd = pd %>% left_join(pd_all[,c("cell_id","celltype_new","system","meta_group")], by = "cell_id") %>% as.data.frame()
    
    if(system_i == "Neuroectoderm"){
        pd$system = "Neuroectoderm"
    }
    
    result = list()
    
    for(cnt in c(1:100)){
        
        print(cnt)
        
        idx = read.csv(paste0(work_path, system_i, "_idx_", cnt, ".csv"), header = F)
        knn = read.csv(paste0(work_path, system_i, "_knn_", cnt, ".csv"), header = F)
        
        pd_x = pd[as.vector(idx$V1) + 1,]
        knn = as.matrix(knn)
        knn = knn + 1
        
        x = data.frame(i = rep(1:nrow(knn), ncol(knn)),
                       j = c(knn), stringsAsFactors = FALSE)
        
        if(max(x$j) != nrow(knn)){
            x = rbind(x, data.frame(i = nrow(knn),
                                    j = nrow(knn), stringsAsFactors = FALSE))
        }
        
        dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                                   j = as.numeric(as.vector(x$j)),
                                   x = 1)
        
        nodes = pd_x %>% group_by(meta_group) %>% tally() %>% rename(celltype_num = n)
        
        dat_t = t(dat) + dat
        x = data.frame(summary(dat_t))
        x = x[x$x == 2 & x$i > x$j,]
        x$x = NULL   
        
        y = data.frame(i = 1:nrow(pd_x),
                       j = 1:nrow(pd_x),
                       meta_group = as.vector(pd_x$meta_group), stringsAsFactors = FALSE)
        
        dat = x %>% left_join(y %>% select(i, meta_group), by = "i") %>%
            left_join(y %>% select(j, meta_group), by = "j") %>%
            group_by(meta_group.x, meta_group.y) %>% tally() 
        
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
            for(j in 1:(i-1)){
                group = rbind(group, data.frame(system = system_i,
                                                x = colnames(obs_x)[j],
                                                y = rownames(obs_x)[i], stringsAsFactors = F))
            }
        }
        
        group$edge_num = obs_y
        
        group = group %>% 
            left_join(nodes %>% rename(x = meta_group) %>% select(x, celltype_num), by = "x") %>% rename(x_size = celltype_num) %>%
            left_join(nodes %>% rename(y = meta_group) %>% select(y, celltype_num), by = "y") %>% rename(y_size = celltype_num)
        
        group$min_size = if_else(group$x_size < group$y_size, group$x_size, group$y_size)
        group$edge_num_norm = group$edge_num/log2(15*group$min_size)
        group$min_size = NULL
        
        result[[cnt]] = group
    }
    
    saveRDS(result, paste0(work_path, system_i, "_edges.rds"))
    
}


### Making plot 

df = NULL

for(kk in 1:length(system_list)){
    system_i = system_list[kk]
    print(system_i)
    
    ### *.edges_new.rds is calculated by step2_Late_stage_graph.R, which is based on the full dataset
    obs = readRDS(paste0(work_path, system_i, ".edges_new.rds"))
    perm = readRDS(paste0(work_path, system_i, "_edges.rds"))
    
    res = NULL
    for(cnt in 1:100){
        tmp = obs %>% rename(edge_num_norm_obs = edge_num_norm) %>%
            left_join(perm[[cnt]][,c("x","y","edge_num_norm")], by = c("x","y"))
        res = c(res, cor.test(tmp$edge_num_norm_obs, tmp$edge_num_norm, method = "spearman")$estimate)
    }
    
    df = rbind(df, data.frame(corr = as.vector(res),
                              system = system_i, stringsAsFactors = F))
}

df$system = factor(df$system, levels = system_list)
p = ggplot(df,aes(system, corr, fill = system)) +
    geom_boxplot() +
    labs(x="System", y="Spearman correlation of normalized # of MNNs between cell types in 80% subsampling and full size", title="") +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=subsystem_color_plate) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 90), axis.text.y = element_text(color="black")) 



##################################################################################
### Third, to determine the effect of k parameter choice on the MNNs identified 
### between cell types, we examined different k values (k = 5, 10, 20, 30, 40, 50)
##################################################################################

### python Graph_robust.py

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(work_path, "df_cell_graph.rds")

system_list = c("Endothelium",                    
                "Epithelial_cells",               
                "Eye",                            
                "Gut",                            
                "Notochord",                      
                "PNS_glia",                       
                "PNS_neurons",                    
                "Renal",                          
                "Lateral_plate_mesoderm",         
                "Blood",     
                "Neuroectoderm",
                "Brain_spinal_cord",              
                "Mesoderm")     

for(kk in 1:length(system_list)){
    
    system_i = system_list[kk]
    print(system_i)
    
    pd = readRDS(paste0(work_path, "/graph/system/", system_i, "_adata_scale.obs.rds"))
    rownames(pd) = as.vector(pd$cell_id)
    pd$celltype_new = pd$system = pd$meta_group = NULL
    pd = pd %>% left_join(pd_all[,c("cell_id","celltype_new","system","meta_group")], by = "cell_id") %>% as.data.frame()
    
    if(system_i == "Neuroectoderm"){
        pd$system = "Neuroectoderm"
    }
    
    result = list()
    
    for(cnt in c(5,10,20,30,40,50)){
        
        print(cnt)
        
        knn = read.csv(paste0(work_path, system_i, "_knn_", cnt, ".csv"), header = F)
        
        pd_x = pd
        knn = as.matrix(knn)
        knn = knn + 1
        
        x = data.frame(i = rep(1:nrow(knn), ncol(knn)),
                       j = c(knn), stringsAsFactors = FALSE)
        
        if(max(x$j) != nrow(knn)){
            x = rbind(x, data.frame(i = nrow(knn),
                                    j = nrow(knn), stringsAsFactors = FALSE))
        }
        
        dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                                   j = as.numeric(as.vector(x$j)),
                                   x = 1)
        
        nodes = pd_x %>% group_by(meta_group) %>% tally() %>% rename(celltype_num = n)
        
        dat_t = t(dat) + dat
        x = data.frame(summary(dat_t))
        x = x[x$x == 2 & x$i > x$j,]
        x$x = NULL   
        
        y = data.frame(i = 1:nrow(pd_x),
                       j = 1:nrow(pd_x),
                       meta_group = as.vector(pd_x$meta_group), stringsAsFactors = FALSE)
        
        dat = x %>% left_join(y %>% select(i, meta_group), by = "i") %>%
            left_join(y %>% select(j, meta_group), by = "j") %>%
            group_by(meta_group.x, meta_group.y) %>% tally() 
        
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
            for(j in 1:(i-1)){
                group = rbind(group, data.frame(system = system_i,
                                                x = colnames(obs_x)[j],
                                                y = rownames(obs_x)[i], stringsAsFactors = F))
            }
        }
        
        group$edge_num = obs_y
        
        group = group %>% 
            left_join(nodes %>% rename(x = meta_group) %>% select(x, celltype_num), by = "x") %>% rename(x_size = celltype_num) %>%
            left_join(nodes %>% rename(y = meta_group) %>% select(y, celltype_num), by = "y") %>% rename(y_size = celltype_num)
        
        group$min_size = if_else(group$x_size < group$y_size, group$x_size, group$y_size)
        group$edge_num_norm = group$edge_num/log2(15*group$min_size)
        group$min_size = NULL
        
        result[[cnt]] = group
    }
    
    saveRDS(result, paste0(work_path, system_i, "_edges.rds"))
    
}

### Making plot

df = NULL

for(kk in 1:length(system_list)){
    system_i = system_list[kk]
    print(system_i)
    
    ### *.edges_new.rds is calculated by step2_Late_stage_graph.R, which is based on k = 15
    obs = readRDS(paste0(work_path, system_i, ".edges_new.rds"))
    perm = readRDS(paste0(work_path, system_i, "_edges.rds"))
    
    res = NULL
    for(cnt in c(5,10,20,30,40,50)){
        tmp = obs %>% rename(edge_num_norm_obs = edge_num_norm) %>%
            left_join(perm[[cnt]][,c("x","y","edge_num_norm")], by = c("x","y"))
        res = c(res, cor.test(tmp$edge_num_norm_obs, tmp$edge_num_norm, method = "spearman")$estimate)
    }
    
    df = rbind(df, data.frame(corr = as.vector(res),
                              system = system_i, 
                              k = c(5,10,20,30,40,50), stringsAsFactors = F))
}

df$k = factor(df$k)
p = df %>% 
    ggplot(aes(k, corr, color = k)) + geom_point(size = 3.5) + 
    facet_wrap(~system, nrow = 2) +
    labs(x="k in kNN", y="Spearman correlation of normalized MNNs between cell types in variant k and original k", title="") +
    theme_classic(base_size = 10) +
    scale_color_viridis(discrete=TRUE) +
    #  theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 


