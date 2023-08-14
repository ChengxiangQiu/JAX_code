
#####################################
### Section - 6, Development tree ###
#####################################

###################################
### Summary the edges and nodes ###
###################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

nodes = read.table(paste0(work_path, "nodes.txt"), header=T, as.is=T, sep="\t")

### now we merged edges which have been manually reviewed.

### edges_1 includes edges from pre-gastrulation and gastrulation stages
edges_1 = read.table(paste0(work_path, "edges_1.txt"), header=F, as.is=T, sep="\t")

### edges_2 includes edges from organogenesis & fetal development
edges_2 = read.table(paste0(work_path, "edges_2.txt"), header=F, as.is=T, sep="\t")

### edges_3 includes edges which are manually added to connect blood and PNS-neuron
edges_3 = read.table(paste0(work_path, "edges_3.txt"), header=F, as.is=T, sep="\t")

edges = rbind(edges_1, edges_2, edges_3)
names(edges) = c("system", "x", "y", "x_name", "y_name", "edge_type")

length((unique(c(edges$x, edges$y))))

write.table(edges, paste0(work_path, "edges.txt"), row.names=F, sep="\t", quote=F)

### To better visualize the result, we took out the spatial continuity edges, and also collapse reundant nodes
edges_sub = edges[edges$edge_type != "Spatial continuity",]
length((unique(c(edges_sub$x, edges_sub$y))))

edges_sub = rbind(edges_sub, edges[edges$x %in% c("BS_M37", "BS_M39") | edges$y %in% c("BS_M37", "BS_M39"),])
length((unique(c(edges_sub$x, edges_sub$y))))

write.table(edges_sub, paste0(work_path, "edges_sub.txt"), row.names=F, sep="\t", quote=F)

### removing redundant nodes
edges_sub$x_y = paste0(edges_sub$x, ":", edges_sub$y)
edges_x_1 = edges_sub[edges_sub$x_name == edges_sub$y_name & edges_sub$system == "Pre_gastrulation",]
edges_x_2 = edges_sub[edges_sub$x_name == edges_sub$y_name & edges_sub$system == "Gastrulation_E8.5b",]

edges_x_3 = edges_sub[edges_sub$x %in% as.vector(edges_x_1$y),]
edges_x_4 = edges_sub[edges_sub$x %in% as.vector(edges_x_2$y),]

edges_x_3_ = edges_x_3 %>% left_join(edges_x_1 %>% select(x,y) %>% rename(new_x = x, x=y), by = "x")
edges_x_3$x = as.vector(edges_x_3_$new_x)

edges_x_4_ = edges_x_4 %>% left_join(edges_x_2 %>% select(x,y) %>% rename(new_x = x, x=y), by = "x")
edges_x_4$x = as.vector(edges_x_4_$new_x)

edges_x_5 = edges_sub[!edges_sub$x_y %in% c(edges_x_1$x_y, edges_x_2$x_y, edges_x_3$x_y, edges_x_4$x_y),]

edges_x = rbind(edges_x_3, edges_x_4, edges_x_5)
print(edges_x[edges_x$x_name == edges_x$y_name,])
edges_x = edges_x[edges_x$x_name != edges_x$y_name,]
edges_x$x_y_name = paste0(edges_x$x_name, ":", edges_x$y_name)
x_table = table(edges_x$x_y_name)
tmp = edges_x[edges_x$x_y_name %in% names(x_table)[x_table != 1],]
print(tmp[order(tmp$x_name),])

redundant_edges = c("En_M5:En_M1", "Ga_M5:Ga_M6", "L_M7:L_M3", "En_M7:En_M5", "Ga_M23:En_M5", "BS_M20:BS_M2", "Ga_M17:En_M7")
edges_x = edges_x[!edges_x$x_y %in% redundant_edges,]
print(length(unique(c(edges_x$x, edges_x$y))))
print(length(unique(c(edges_x$x_name, edges_x$y_name))))

write.table(edges_x, paste0(work_path, "edges_sub.txt"), row.names=F, sep="\t", quote=F)

nodes_sub = nodes[nodes$meta_group %in% c(edges_x$x, edges_x$y),]
write.table(nodes_sub, paste0(work_path, "nodes_sub.txt"), row.names=F, sep="\t", quote=F)


##############################################################
### making Histogram for accepted edges and rejected edges ###
##############################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

dat = read.table(paste0(work_path, "edges_MNNs.txt"), header=T, sep="\t", as.is=T)
dat = dat[dat$MNN_pairs_normalized > 1,]

dat_1 = dat[dat$Comments %in% c("Developmental progression", "Spatial continuity"),]
dat_2 = dat[dat$Comments %in% c("x","X"),]

dat_uniq = NULL
x_uniq = NULL
for(i in 1:nrow(dat_2)){
    tmp = paste0(dat_2$x[i], ":", dat_2$y[i])
    tmp_r = paste0(dat_2$y[i], ":", dat_2$x[i])
    if(tmp %in% x_uniq | tmp_r %in% x_uniq){
        next
    } else {
        dat_uniq = rbind(dat_uniq, dat_2[i,])
        x_uniq = c(x_uniq, tmp)
    }
}

dat_1$group = "Accepted"
dat_uniq$group = "Rejected"
df = rbind(dat_1, dat_uniq)
df$log2_MNN_pairs_normalized = log2(df$MNN_pairs_normalized)

### Fig. S21a

p <- df %>%
    ggplot( aes(x=log2_MNN_pairs_normalized, fill=group)) +
    geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
    scale_fill_manual(values=c("#f85633", "#0058d6")) +
    theme_ipsum() +
    labs(fill="")



