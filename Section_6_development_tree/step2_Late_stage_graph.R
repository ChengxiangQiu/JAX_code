
#####################################
### Section - 6, Development tree ###
#####################################

#########################################################################################################
### First, we manually split all the cell types, from the organogenesis & fetal development, into 12 systems, 
### to perform dimension reducting using Scanpy, followed by identifying the kNNs across cells using annoy in Python

### For Brain_spinal_cord, we further split patterned neuroectoderm ("Neuroectoderm") to perform embedding.

### Python Dimension_reduction_subsystem.py

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell_graph.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")

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

patterned_neuroectoderm = c("Anterior floor plate",
                            "Diencephalon",
                            "Floorplate and p3 domain",
                            "Hypothalamus",
                            "Midbrain",
                            "Posterior roof plate",
                            "Telencephalon",
                            "Anterior roof plate",
                            "Dorsal telencephalon",
                            "Hindbrain",
                            "Hypothalamus (Sim1+)",
                            "Midbrain-hindbrain boundary",
                            "Spinal cord/r7/r8")

for(kk in 1:length(system_list)){
    
    system_i = system_list[kk]
    print(system_i)
    
    ### After you running the "Python Dimension_reduction_subsystem.py", you will get this files
    pd = read.csv(paste0(work_path, system_i, "_adata_scale.obs.csv"), as.is=T, row.names = 1)
    rownames(pd) = as.vector(pd$cell_id)
    pd = pd %>% left_join(pd_all[,c("cell_id","celltype_new","system","meta_group")], by = "cell_id") %>% as.data.frame()
    
    if(system_i == "Neuroectoderm"){
        pd$system = "Neuroectoderm"
    }
    
    ### create or read MNN pairs between individual cells
    
    nn_matrix = read.csv(paste0(work_path, system_i, "_adata_scale.kNN_15.csv"), as.is=T, header=F)
    nn_matrix = as.matrix(nn_matrix)
    nn_matrix = nn_matrix + 1 ### python and R using different start index
    
    ### extracting MNN pairs
    ### only retaining those edges which are considered twice (A -> B, B -> A)
    
    x = data.frame(i = rep(1:nrow(nn_matrix), ncol(nn_matrix)),
                   j = c(nn_matrix), stringsAsFactors = FALSE)
    
    dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                               j = as.numeric(as.vector(x$j)),
                               x = 1)
    
    dat_t = t(dat) + dat
    x = data.frame(summary(dat_t))
    x = x[x$x == 2 & x$i > x$j,]
    x$x = NULL   
    
    ### x saves the MNN pairs
    saveRDS(x, paste0(work_path, system_i, ".MNN.rds"))
    
    y = data.frame(i = 1:nrow(pd),
                   j = 1:nrow(pd),
                   meta_group = as.vector(pd$meta_group), stringsAsFactors = FALSE)
    
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
        print(i)
        for(j in 1:(i-1)){
            group = rbind(group, data.frame(system = system_i,
                                            x = colnames(obs_x)[j],
                                            y = rownames(obs_x)[i], stringsAsFactors = F))
        }
    }
    
    group$edge_num = obs_y
    
    group = group %>% 
        left_join(nodes %>% rename(x = meta_group) %>% select(x, celltype_new, celltype_num), by = "x") %>% rename(x_name = celltype_new, x_size = celltype_num) %>%
        left_join(nodes %>% rename(y = meta_group) %>% select(y, celltype_new, celltype_num), by = "y") %>% rename(y_name = celltype_new, y_size = celltype_num)
    
    group$min_size = if_else(group$x_size < group$y_size, group$x_size, group$y_size)
    group$edge_num_norm = group$edge_num/log2(15*group$min_size)
    group$min_size = NULL
    
    saveRDS(group, paste0(work_path, system_i, ".edges_new.rds"))
    
    ### output MNN pairs for manually reviewing
    
    edges = group
    
    edges_2 = edges
    edges_2$x = as.vector(edges$y); edges_2$y = as.vector(edges$x)
    edges_2$x_size = as.vector(edges$y_size); edges_2$y_size = as.vector(edges$x_size)
    edges_2$x_name = as.vector(edges$y_name); edges_2$y_name = as.vector(edges$x_name)
    
    dat = rbind(edges, edges_2) %>% as.data.frame() %>%
        filter(edge_num != 0) %>% rename(MNN_pairs = edge_num, MNN_pairs_normalized = edge_num_norm) %>%
        group_by(x) %>% arrange(desc(MNN_pairs), .by_group = T) %>%
        as.data.frame()
    
    if(system_i == "Brain_spinal_cord"){
        tmp = read.table(paste0(work_path, "Neuroectoderm", ".MNN_pairs.txt"),header=T,as.is=T,sep="\t")
        dat = dat[!dat$x %in% c(tmp$x, tmp$y) | !dat$y %in% c(tmp$x, tmp$y), ]
    }
    
    write.table(dat, paste0(work_path, system_i, ".MNN_pairs.txt"), row.names=F, sep="\t", quote=F)

}




