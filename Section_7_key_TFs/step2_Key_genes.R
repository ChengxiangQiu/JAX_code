

####################################
### Section - 7, Key TFs & genes ###
####################################

###########################################
### systemactially nominating key genes ###
###########################################


### systemactially nominating key Genes

# 	Focusing on A, Gene is highly expressed in MNN cells relative to the left cells.
# 	Focusing on A and B, Gene is highly expressed in MNN cells of B than MNN cells of A.
# 	Focusing on B, Gene is highly expressed in MNN cells relative to the left cells.


####################
### Gastrulation ###
####################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")

system_i = "Gastrulation"

edges = read.table(paste0(work_path, "edges.txt"), header=F, as.is=T, sep="\t")
names(edges) = c("system","x","y","x_name","y_name","edge_type")
edges = edges[edges$system == "Gastrulation",]

obj = readRDS(paste0(work_path, "obj_Early_PS.rds"))
pd = readRDS(paste0(work_path, system_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
pd_all = readRDS(paste0(work_path, "df_cell_graph.rds"))
pd$celltype_new = NULL
pd = pd %>% left_join(pd_all[,c("cell_id", "celltype_new")], by = "cell_id")

if ("meta_group" %in% names(pd)) {pd$meta_group = NULL}
pd_x = pd %>% left_join(nodes %>% filter(system == system_i) %>% select(celltype_new, meta_group))
pd$meta_group = as.vector(pd_x$meta_group)
gene_count = GetAssayData(obj, slot = "counts")

y = data.frame(i = 1:nrow(pd),
               j = 1:nrow(pd),
               meta_group = as.vector(pd$meta_group), stringsAsFactors = FALSE)

x = readRDS(paste0(work_path, system_i, ".MNN.rds"))
x_rev = x
x_rev$i = x$j
x_rev$j = x$i
x_rev = rbind(x, x_rev)
dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j")

res_all = NULL

for(cnt in 1:nrow(edges)){
    
    print(cnt)
    
    xx = as.vector(edges$x)[cnt]
    yy = as.vector(edges$y)[cnt]
    
    dat_cnt = dat[dat$meta_group.x == xx & dat$meta_group.y == yy,]
    group_1_MNN = as.vector(dat_cnt$i)
    group_2_MNN = as.vector(dat_cnt$j)
    
    coor = c(1:nrow(pd))
    
    while(length(unique(group_1_MNN)) < 200){
        num_1 = length(unique(group_1_MNN))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
        group_1_MNN = c(group_1_MNN, tmp)
        num_2 = length(unique(group_1_MNN))
        if(num_1 == num_2) {break}
    }
    
    group_1_close = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
    group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
    while(length(unique(group_1_close)) < 200){
        num_1 = length(unique(group_1_close))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_1_close]), coor[pd$meta_group == xx])
        group_1_close = c(group_1_close, tmp)
        group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
        num_2 = length(unique(group_1_close))
        if(num_1 == num_2) {break}
    }
    
    while(length(unique(group_2_MNN)) < 200){
        num_1 = length(unique(group_2_MNN))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
        group_2_MNN = c(group_2_MNN, tmp)
        num_2 = length(unique(group_2_MNN))
        if(num_1 == num_2) {break}
    }
    
    group_2_close = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
    group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
    while(length(unique(group_2_close)) < 200){
        num_1 = length(unique(group_2_close))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_2_close]), coor[pd$meta_group == yy])
        group_2_close = c(group_2_close, tmp)
        group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
        num_2 = length(unique(group_2_close))
        if(num_1 == num_2) {break}
    }
    
    
    group = rep("other", nrow(pd))
    group[coor %in% group_1_close] = "group_1"
    group[coor %in% group_1_MNN] = "group_2"
    group[coor %in% group_2_MNN] = "group_3"
    group[coor %in% group_2_close] = "group_4"
    pd$group = as.vector(group)
    
    pd_sub = pd[pd$group != "other",]
    group_table = table(pd_sub$group)
    
    saveRDS(group_table, paste0(work_path, system_i, "_", xx, "_", yy, ".group_table.rds"))
    
    if(sum(group_table >= 20) >= 2){
        
        group_oversize = names(group_table)[group_table > 2000]
        if(length(group_oversize) != 0){
            pd_sub_1 = pd_sub[pd_sub$group %in% group_oversize,]
            pd_sub_1_sub = pd_sub_1 %>% group_by(group) %>% sample_n(2000)
            pd_sub_1 = pd_sub_1[pd_sub_1$cell_id %in% pd_sub_1_sub$cell_id,]
            pd_sub_2 = pd_sub[!pd_sub$group %in% group_oversize,]
            pd_sub = rbind(pd_sub_1, pd_sub_2)
        }
        gene_count_sub = gene_count[,as.vector(pd_sub$cell_id)]
        
        obj_sub = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
        obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)
        Idents(obj_sub) = as.vector(obj_sub$group)
        
        
        ### key Gene
        gene_use = rownames(obj_sub)
        res_1 = NULL
        res_2 = NULL
        res_3 = NULL
        
        if(sum(pd_sub$group == "group_1") >= 20 & sum(pd_sub$group == "group_2") >= 20){
            res_1 = FindMarkers(obj_sub, ident.1 = "group_1", ident.2 = "group_2", logfc.threshold = 0, min.pct = 0, features = gene_use)
            if(nrow(res_1) != 0){
                res_1 = res_1 %>% mutate(gene_ID = rownames(res_1)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                    filter(p_val_adj < 0.05)
                if(nrow(res_1) != 0){
                    res_1$comparing = "group_1 vs. group_2"
                    res_1$high_in_which = if_else(res_1$avg_logFC > 0, "down in early transition", "up in early transition")
                    res_1$node_A = xx
                    res_1$node_B = yy
                    res_1$celltype_A = edges$x_name[cnt]
                    res_1$celltype_B = edges$y_name[cnt]
                    res_1$edge_type = edges$edge_type[cnt]
                    res_1 = res_1[order(res_1$high_in_which),]
                } else {res_1 = NULL}
            }
        }
        
        if(sum(pd_sub$group == "group_2") >= 20 & sum(pd_sub$group == "group_3") >= 20){
            res_2 = FindMarkers(obj_sub, ident.1 = "group_2", ident.2 = "group_3", logfc.threshold = 0, min.pct = 0, features = gene_use)
            if(nrow(res_2) != 0){
                res_2 = res_2 %>% mutate(gene_ID = rownames(res_2)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                    filter(p_val_adj < 0.05)
                if(nrow(res_2) != 0){
                    res_2$comparing = "group_2 vs. group_3"
                    res_2$high_in_which = if_else(res_2$avg_logFC > 0, "down in transition", "up in transition")
                    res_2$node_A = xx
                    res_2$node_B = yy
                    res_2$celltype_A = edges$x_name[cnt]
                    res_2$celltype_B = edges$y_name[cnt]
                    res_2$edge_type = edges$edge_type[cnt]
                    res_2 = res_2[order(res_2$high_in_which),]
                } else {res_2 = NULL}
            }
        }
        
        if(sum(pd_sub$group == "group_3") >= 20 & sum(pd_sub$group == "group_4") >= 20){
            res_3 = FindMarkers(obj_sub, ident.1 = "group_3", ident.2 = "group_4", logfc.threshold = 0, min.pct = 0, features = gene_use)
            if(nrow(res_3) != 0){
                res_3 = res_3 %>% mutate(gene_ID = rownames(res_3)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                    filter(p_val_adj < 0.05)
                if(nrow(res_3) != 0){
                    res_3$comparing = "group_3 vs. group_4"
                    res_3$high_in_which = if_else(res_3$avg_logFC > 0, "down in late transition", "up in late transition")
                    res_3$node_A = xx
                    res_3$node_B = yy
                    res_3$celltype_A = edges$x_name[cnt]
                    res_3$celltype_B = edges$y_name[cnt]
                    res_3$edge_type = edges$edge_type[cnt]
                    res_3 = res_3[order(res_3$high_in_which),]
                } else {res_3 = NULL}
            }
        }
        
        if(!is.null(res_1) | !is.null(res_2) | !is.null(res_3)){
            saveRDS(rbind(res_1, res_2, res_3), paste0(work_path, system_i, "_", xx, "_", yy, ".rds"))
            res_all = rbind(res_all, rbind(res_1, res_2, res_3))
        }
        
    }
    
}

saveRDS(res_all, paste0(work_path, system_i, ".keyGene.rds"))
write.csv(res_all, paste0(work_path, system_i, ".keyGene.csv"))




#########################################
### Organogenesis & fetal development ###
#########################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

system_list = c("Endothelium",             ### En
                "Epithelial_cells",        ### Ep
                "Eye",                     ### Ey
                "Gut",                     ### G
                "Notochord",               ### No
                "PNS_glia",                ### PG
                "PNS_neurons",             ### PN
                "Renal",                   ### R
                "Lateral_plate_mesoderm",  ### L
                "Blood",                   ### B
                "Brain_spinal_cord",       ### BS
                "Mesoderm")                ### M

for(kk in 1:length(system_list)){
    
    system_i = system_list[kk]
    print(system_i)
    
    nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")
    
    neuroectoderm_nodes = c("BS_M20", "BS_M40", "BS_M22", "BS_M54", "BS_M14", "BS_M25", "BS_M27", "BS_M15", "BS_M34", "BS_M23", "BS_M35", "BS_M2", "BS_M1")
    
    edges = read.table(paste0(work_path, "edges.txt"), header=F, as.is=T, sep="\t")
    names(edges) = c("system","x","y","x_name","y_name","edge_type")
    edges = edges[edges$system == system_i,]
    
    pd = readRDS(paste0(work_path, system_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
    pd$celltype_new = pd$celltype_update = NULL
    
    pd_all = readRDS(paste0(work_path, "df_cell_graph.rds"))
    pd = pd %>% left_join(pd_all[,c("cell_id", "celltype_new")], by = "cell_id")
    
    pd_all = readRDS(paste0(work_path, "df_cell.rds"))
    pd = pd %>% left_join(pd_all[,c("cell_id", "celltype_update")], by = "cell_id")
    
    if ("meta_group" %in% names(pd)) {pd$meta_group = NULL}
    pd_x = pd %>% left_join(nodes %>% filter(system == system_i) %>% select(celltype_new, meta_group))
    pd$meta_group = as.vector(pd_x$meta_group)
    
    y = data.frame(i = 1:nrow(pd),
                   j = 1:nrow(pd),
                   meta_group = as.vector(pd$meta_group), stringsAsFactors = FALSE)
    
    x = readRDS(paste0(work_path, system_i, ".MNN.rds"))
    x_rev = x
    x_rev$i = x$j
    x_rev$j = x$i
    x_rev = rbind(x, x_rev)
    dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
        left_join(y %>% select(j, meta_group), by = "j")
    
    res_all = NULL
    
    for(cnt in 1:nrow(edges)){
        
        print(cnt)
        
        xx = as.vector(edges$x)[cnt]
        yy = as.vector(edges$y)[cnt]
        
        if(xx %in% neuroectoderm_nodes & yy %in% neuroectoderm_nodes){
            next
        }
        
        dat_cnt = dat[dat$meta_group.x == xx & dat$meta_group.y == yy,]
        group_1_MNN = as.vector(dat_cnt$i)
        group_2_MNN = as.vector(dat_cnt$j)
        
        coor = c(1:nrow(pd))
        
        while(length(unique(group_1_MNN)) < 200){
            num_1 = length(unique(group_1_MNN))
            tmp = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
            group_1_MNN = c(group_1_MNN, tmp)
            num_2 = length(unique(group_1_MNN))
            if(num_1 == num_2) {break}
        }
        
        group_1_close = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
        group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
        while(length(unique(group_1_close)) < 200){
            num_1 = length(unique(group_1_close))
            tmp = intersect(as.vector(dat$j[dat$i %in% group_1_close]), coor[pd$meta_group == xx])
            group_1_close = c(group_1_close, tmp)
            group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
            num_2 = length(unique(group_1_close))
            if(num_1 == num_2) {break}
        }
        
        while(length(unique(group_2_MNN)) < 200){
            num_1 = length(unique(group_2_MNN))
            tmp = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
            group_2_MNN = c(group_2_MNN, tmp)
            num_2 = length(unique(group_2_MNN))
            if(num_1 == num_2) {break}
        }
        
        group_2_close = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
        group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
        while(length(unique(group_2_close)) < 200){
            num_1 = length(unique(group_2_close))
            tmp = intersect(as.vector(dat$j[dat$i %in% group_2_close]), coor[pd$meta_group == yy])
            group_2_close = c(group_2_close, tmp)
            group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
            num_2 = length(unique(group_2_close))
            if(num_1 == num_2) {break}
        }
        
        
        group = rep("other", nrow(pd))
        group[coor %in% group_1_close] = "group_1"
        group[coor %in% group_1_MNN] = "group_2"
        group[coor %in% group_2_MNN] = "group_3"
        group[coor %in% group_2_close] = "group_4"
        pd$group = as.vector(group)
        
        pd_sub = pd[pd$group != "other",]
        group_table = table(pd_sub$group)
        
        saveRDS(group_table, paste0(work_path, system_i, "_", xx, "_", yy, ".group_table.rds"))
        
        if(sum(group_table >= 20) >= 2){
            
            group_oversize = names(group_table)[group_table > 2000]
            if(length(group_oversize) != 0){
                pd_sub_1 = pd_sub[pd_sub$group %in% group_oversize,]
                pd_sub_1_sub = pd_sub_1 %>% group_by(group) %>% sample_n(2000)
                pd_sub_1 = pd_sub_1[pd_sub_1$cell_id %in% pd_sub_1_sub$cell_id,]
                pd_sub_2 = pd_sub[!pd_sub$group %in% group_oversize,]
                pd_sub = rbind(pd_sub_1, pd_sub_2)
            }
            
            mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M","X","Y")),]
            gene_count_sub = doExtractData(pd_sub, mouse_gene_sub)
            
            obj_sub = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
            obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)
            Idents(obj_sub) = as.vector(obj_sub$group)
            
            ### key Gene
            gene_use = rownames(obj_sub)
            res_1 = NULL
            res_2 = NULL
            res_3 = NULL
            
            if(sum(pd_sub$group == "group_1") >= 20 & sum(pd_sub$group == "group_2") >= 20){
                res_1 = FindMarkers(obj_sub, ident.1 = "group_1", ident.2 = "group_2", logfc.threshold = 0, min.pct = 0, features = gene_use)
                if(nrow(res_1) != 0){
                    res_1 = res_1 %>% mutate(gene_ID = rownames(res_1)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                        filter(p_val_adj < 0.05)
                    if(nrow(res_1) != 0){
                        res_1$comparing = "group_1 vs. group_2"
                        res_1$high_in_which = if_else(res_1$avg_logFC > 0, "down in early transition", "up in early transition")
                        res_1$node_A = xx
                        res_1$node_B = yy
                        res_1$celltype_A = edges$x_name[cnt]
                        res_1$celltype_B = edges$y_name[cnt]
                        res_1$edge_type = edges$edge_type[cnt]
                        res_1 = res_1[order(res_1$high_in_which),]
                    } else {res_1 = NULL}
                }
            }
            
            if(sum(pd_sub$group == "group_2") >= 20 & sum(pd_sub$group == "group_3") >= 20){
                res_2 = FindMarkers(obj_sub, ident.1 = "group_2", ident.2 = "group_3", logfc.threshold = 0, min.pct = 0, features = gene_use)
                if(nrow(res_2) != 0){
                    res_2 = res_2 %>% mutate(gene_ID = rownames(res_2)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                        filter(p_val_adj < 0.05)
                    if(nrow(res_2) != 0){
                        res_2$comparing = "group_2 vs. group_3"
                        res_2$high_in_which = if_else(res_2$avg_logFC > 0, "down in transition", "up in transition")
                        res_2$node_A = xx
                        res_2$node_B = yy
                        res_2$celltype_A = edges$x_name[cnt]
                        res_2$celltype_B = edges$y_name[cnt]
                        res_2$edge_type = edges$edge_type[cnt]
                        res_2 = res_2[order(res_2$high_in_which),]
                    } else {res_2 = NULL}
                }
            }
            
            if(sum(pd_sub$group == "group_3") >= 20 & sum(pd_sub$group == "group_4") >= 20){
                res_3 = FindMarkers(obj_sub, ident.1 = "group_3", ident.2 = "group_4", logfc.threshold = 0, min.pct = 0, features = gene_use)
                if(nrow(res_3) != 0){
                    res_3 = res_3 %>% mutate(gene_ID = rownames(res_3)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                        filter(p_val_adj < 0.05)
                    if(nrow(res_3) != 0){
                        res_3$comparing = "group_3 vs. group_4"
                        res_3$high_in_which = if_else(res_3$avg_logFC > 0, "down in late transition", "up in late transition")
                        res_3$node_A = xx
                        res_3$node_B = yy
                        res_3$celltype_A = edges$x_name[cnt]
                        res_3$celltype_B = edges$y_name[cnt]
                        res_3$edge_type = edges$edge_type[cnt]
                        res_3 = res_3[order(res_3$high_in_which),]
                    } else {res_3 = NULL}
                }
            }
            
            if(!is.null(res_1) | !is.null(res_2) | !is.null(res_3)){
                saveRDS(rbind(res_1, res_2, res_3), paste0(work_path, system_i, "_", xx, "_", yy, ".rds"))
                res_all = rbind(res_all, rbind(res_1, res_2, res_3))
            }
            
        }
        
    }
    
    saveRDS(res_all, paste0(work_path, system_i, ".keyGene.rds"))
    write.csv(res_all, paste0(work_path, system_i, ".keyGene.csv"))
    
    
}






#####################
### Neuroectoderm ###
#####################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

system_i = "Brain_spinal_cord"
print(system_i)

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")

neuroectoderm_nodes = c("BS_M20", "BS_M40", "BS_M22", "BS_M54", "BS_M14", "BS_M25", "BS_M27", "BS_M15", "BS_M34", "BS_M23", "BS_M35", "BS_M2", "BS_M1")

edges = read.table(paste0(work_path, "edges.txt"), header=F, as.is=T, sep="\t")
names(edges) = c("system","x","y","x_name","y_name","edge_type")
edges = edges[edges$system == system_i,]

pd = readRDS(paste0(work_path, "Neuroectoderm", "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
pd$celltype_new = pd$celltype_update = NULL

pd_all = readRDS(paste0(work_path, "df_cell_graph.rds"))
pd = pd %>% left_join(pd_all[,c("cell_id", "celltype_new")], by = "cell_id")

pd_all = readRDS(paste0(work_path, "df_cell.rds"))
pd = pd %>% left_join(pd_all[,c("cell_id", "celltype_update")], by = "cell_id")

if ("meta_group" %in% names(pd)) {pd$meta_group = NULL}
pd_x = pd %>% left_join(nodes %>% filter(system == system_i) %>% select(celltype_new, meta_group))
pd$meta_group = as.vector(pd_x$meta_group)

y = data.frame(i = 1:nrow(pd),
               j = 1:nrow(pd),
               meta_group = as.vector(pd$meta_group), stringsAsFactors = FALSE)

x = readRDS(paste0(work_path, "Neuroectoderm", ".MNN.rds"))
x_rev = x
x_rev$i = x$j
x_rev$j = x$i
x_rev = rbind(x, x_rev)
dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
    left_join(y %>% select(j, meta_group), by = "j")

res_all = NULL

for(cnt in 1:nrow(edges)){
    
    print(cnt)
    
    xx = as.vector(edges$x)[cnt]
    yy = as.vector(edges$y)[cnt]
    
    if(!xx %in% neuroectoderm_nodes | !yy %in% neuroectoderm_nodes){
        next
    }
    
    dat_cnt = dat[dat$meta_group.x == xx & dat$meta_group.y == yy,]
    group_1_MNN = as.vector(dat_cnt$i)
    group_2_MNN = as.vector(dat_cnt$j)
    
    coor = c(1:nrow(pd))
    
    while(length(unique(group_1_MNN)) < 200){
        num_1 = length(unique(group_1_MNN))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
        group_1_MNN = c(group_1_MNN, tmp)
        num_2 = length(unique(group_1_MNN))
        if(num_1 == num_2) {break}
    }
    
    group_1_close = intersect(as.vector(dat$j[dat$i %in% group_1_MNN]), coor[pd$meta_group == xx])
    group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
    while(length(unique(group_1_close)) < 200){
        num_1 = length(unique(group_1_close))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_1_close]), coor[pd$meta_group == xx])
        group_1_close = c(group_1_close, tmp)
        group_1_close = group_1_close[!group_1_close %in% group_1_MNN]
        num_2 = length(unique(group_1_close))
        if(num_1 == num_2) {break}
    }
    
    while(length(unique(group_2_MNN)) < 200){
        num_1 = length(unique(group_2_MNN))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
        group_2_MNN = c(group_2_MNN, tmp)
        num_2 = length(unique(group_2_MNN))
        if(num_1 == num_2) {break}
    }
    
    group_2_close = intersect(as.vector(dat$j[dat$i %in% group_2_MNN]), coor[pd$meta_group == yy])
    group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
    while(length(unique(group_2_close)) < 200){
        num_1 = length(unique(group_2_close))
        tmp = intersect(as.vector(dat$j[dat$i %in% group_2_close]), coor[pd$meta_group == yy])
        group_2_close = c(group_2_close, tmp)
        group_2_close = group_2_close[!group_2_close %in% group_2_MNN]
        num_2 = length(unique(group_2_close))
        if(num_1 == num_2) {break}
    }
    
    
    group = rep("other", nrow(pd))
    group[coor %in% group_1_close] = "group_1"
    group[coor %in% group_1_MNN] = "group_2"
    group[coor %in% group_2_MNN] = "group_3"
    group[coor %in% group_2_close] = "group_4"
    pd$group = as.vector(group)
    
    pd_sub = pd[pd$group != "other",]
    group_table = table(pd_sub$group)
    
    saveRDS(group_table, paste0(work_path, system_i, "_", xx, "_", yy, ".group_table.rds"))
    
    if(sum(group_table >= 20) >= 2){
        
        group_oversize = names(group_table)[group_table > 2000]
        if(length(group_oversize) != 0){
            pd_sub_1 = pd_sub[pd_sub$group %in% group_oversize,]
            pd_sub_1_sub = pd_sub_1 %>% group_by(group) %>% sample_n(2000)
            pd_sub_1 = pd_sub_1[pd_sub_1$cell_id %in% pd_sub_1_sub$cell_id,]
            pd_sub_2 = pd_sub[!pd_sub$group %in% group_oversize,]
            pd_sub = rbind(pd_sub_1, pd_sub_2)
        }
        
        mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M","X","Y")),]
        gene_count_sub = doExtractData(pd_sub, mouse_gene_sub)
        
        obj_sub = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
        obj_sub = NormalizeData(obj_sub, normalization.method = "LogNormalize", scale.factor = 10000)
        Idents(obj_sub) = as.vector(obj_sub$group)
        
        
        ### key Gene
        gene_use = rownames(obj_sub)
        res_1 = NULL
        res_2 = NULL
        res_3 = NULL
        
        if(sum(pd_sub$group == "group_1") >= 20 & sum(pd_sub$group == "group_2") >= 20){
            res_1 = FindMarkers(obj_sub, ident.1 = "group_1", ident.2 = "group_2", logfc.threshold = 0, min.pct = 0, features = gene_use)
            if(nrow(res_1) != 0){
                res_1 = res_1 %>% mutate(gene_ID = rownames(res_1)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                    filter(p_val_adj < 0.05)
                if(nrow(res_1) != 0){
                    res_1$comparing = "group_1 vs. group_2"
                    res_1$high_in_which = if_else(res_1$avg_logFC > 0, "down in early transition", "up in early transition")
                    res_1$node_A = xx
                    res_1$node_B = yy
                    res_1$celltype_A = edges$x_name[cnt]
                    res_1$celltype_B = edges$y_name[cnt]
                    res_1$edge_type = edges$edge_type[cnt]
                    res_1 = res_1[order(res_1$high_in_which),]
                } else {res_1 = NULL}
            }
        }
        
        if(sum(pd_sub$group == "group_2") >= 20 & sum(pd_sub$group == "group_3") >= 20){
            res_2 = FindMarkers(obj_sub, ident.1 = "group_2", ident.2 = "group_3", logfc.threshold = 0, min.pct = 0, features = gene_use)
            if(nrow(res_2) != 0){
                res_2 = res_2 %>% mutate(gene_ID = rownames(res_2)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                    filter(p_val_adj < 0.05)
                if(nrow(res_2) != 0){
                    res_2$comparing = "group_2 vs. group_3"
                    res_2$high_in_which = if_else(res_2$avg_logFC > 0, "down in transition", "up in transition")
                    res_2$node_A = xx
                    res_2$node_B = yy
                    res_2$celltype_A = edges$x_name[cnt]
                    res_2$celltype_B = edges$y_name[cnt]
                    res_2$edge_type = edges$edge_type[cnt]
                    res_2 = res_2[order(res_2$high_in_which),]
                } else {res_2 = NULL}
            }
        }
        
        if(sum(pd_sub$group == "group_3") >= 20 & sum(pd_sub$group == "group_4") >= 20){
            res_3 = FindMarkers(obj_sub, ident.1 = "group_3", ident.2 = "group_4", logfc.threshold = 0, min.pct = 0, features = gene_use)
            if(nrow(res_3) != 0){
                res_3 = res_3 %>% mutate(gene_ID = rownames(res_3)) %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>%
                    filter(p_val_adj < 0.05)
                if(nrow(res_3) != 0){
                    res_3$comparing = "group_3 vs. group_4"
                    res_3$high_in_which = if_else(res_3$avg_logFC > 0, "down in late transition", "up in late transition")
                    res_3$node_A = xx
                    res_3$node_B = yy
                    res_3$celltype_A = edges$x_name[cnt]
                    res_3$celltype_B = edges$y_name[cnt]
                    res_3$edge_type = edges$edge_type[cnt]
                    res_3 = res_3[order(res_3$high_in_which),]
                } else {res_3 = NULL}
            }
        }
        
        if(!is.null(res_1) | !is.null(res_2) | !is.null(res_3)){
            saveRDS(rbind(res_1, res_2, res_3), paste0(work_path, system_i, "_", xx, "_", yy, ".rds"))
            res_all = rbind(res_all, rbind(res_1, res_2, res_3))
        }
        
    }
    
}

saveRDS(res_all, paste0(work_path, "Neuroectoderm", ".keyGene.rds"))
write.csv(res_all, paste0(work_path, "Neuroectoderm", ".keyGene.csv"))





######################################################
### Summarize results after finishing the analysis ###
######################################################



system_list = c("Gastrulation",
                "Endothelium",             ### En
                "Epithelial_cells",        ### Ep
                "Eye",                     ### Ey
                "Gut",                     ### G
                "Notochord",               ### No
                "PNS_glia",                ### PG
                "PNS_neurons",             ### PN
                "Renal",                   ### R
                "Lateral_plate_mesoderm",  ### L
                "Blood",                   ### B
                "Neuroectoderm",
                "Brain_spinal_cord",       ### BS
                "Mesoderm")                ### M

result = NULL

for(i in 1:length(system_list)){
    system_i = system_list[i]
    print(paste0(i,"/",system_i))
    
    res = readRDS(paste0(work_path, system_i, ".keyGene.rds"))
    res$id = paste0("r", rep(1:nrow(res)))
    res_sub = res %>% group_by(high_in_which, node_A, node_B)
    res = res[res$id %in% res_sub$id,]
    res$p_val = res$id = NULL
    res$system = system_i
    
    if(system_i == "Neuroectoderm"){
        res$system = "Brain_spinal_cord"
    }
    
    result = rbind(result, res)
}

result_out = result[,c("system","gene_ID","gene_short_name",
                       "avg_logFC","pct.1","pct.2","p_val_adj","comparing",
                       "high_in_which","node_A","node_B","celltype_A","celltype_B","edge_type")]


write.csv(result_out, paste0(work_path, "All", ".keyGene.csv"))



### filtering top 10 significant Genes only

result = NULL

for(i in 1:length(system_list)){
    system_i = system_list[i]
    print(paste0(i,"/",system_i))
    
    res = readRDS(paste0(work_path, system_i, ".keyGene.rds"))
    
    res_x = res %>% group_by(node_A, node_B, celltype_A, celltype_B, edge_type, high_in_which) %>% tally()
    res_y = NULL
    for(j in 1:nrow(res_x)){
        res_j = res %>% filter(high_in_which == res_x$high_in_which[j],
                               node_A == res_x$node_A[j],
                               node_B == res_x$node_B[j]) %>% slice_min(order_by = p_val_adj, n = 10) %>% pull(gene_short_name)
        res_y = c(res_y, paste(as.vector(res_j), collapse = ", "))
    }
    res_x$n = NULL
    res_x$system = system_i
    res_x$keyGene = res_y
    
    if(system_i == "Neuroectoderm"){
        res_x$system = "Brain_spinal_cord"
    }
    
    result = rbind(result, res_x)
}

write.csv(as.data.frame(result), paste0(work_path, "All", ".keyGene.simple.csv"))




