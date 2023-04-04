
###################################
### Section - 1, Basci analysis ###
###################################

####################################
### Making 2D UMAP visualization ###
####################################

source("JAX_help_code.R")
source("JAX_color_code.R")

pd = readRDS("df_cell.rds")
### n = 11,441,407 cells retained

### Making 2D UMAP visualization for the global embedding

### Highlight cells with their major trajectories (Fig. 1f)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = major_trajectory), size=0.3) +
    scale_color_manual(values = major_trajectory_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Global_embedding_2D_UMAP_major_trajectory.png", width = 10, height = 10, dpi = 300)

### Highlight cells with their day timepoints (Fig. 1f)
### Of note, we need to downsample cells from each timepint to a similar number (i.e. 100,000)

pd_1 = pd %>% filter(pd$day %in% c("E8.75", "E17.25")) %>% as.data.frame()
pd_2 = pd %>% filter(!pd$day %in% c("E8.75", "E17.25")) %>% group_by(day) %>% sample_n(100000) %>% as.data.frame()
pd_sub = rbind(pd_1, pd_2)
p = pd_sub[sample(1:nrow(pd_sub))] %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.3) +
    scale_color_manual(values=day_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Global_embedding_2D_UMAP_day.png", width = 10, height = 10, dpi = 300)


###################
### Cell number ###
###################

### Making bar plot for cell number from individual timepoints (Fig. 1c)

pd_cell_num_1 = pd %>% group_by(day) %>% tally() %>% rename(cell_num = n) %>% as.data.frame()
pd_cell_num_1$day = factor(pd_cell_num_1$day, levels = rev(names(day_color_plate)))
pd_cell_num_1$log2_cell_num = log2(pd_cell_num_1$cell_num)

p1 = pd_cell_num_1 %>%
    ggplot(aes(day, cell_num, fill = day)) + 
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_manual(values = day_color_plate) + 
    scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    geom_text(aes(label = scales::comma(cell_num)), 
              hjust = -0.1,
              position = position_dodge(width = 1),
              inherit.aes = TRUE,
              size = 3) +
    labs(x = "", y = "Cell number") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
pdf("Cell_number_timepoints.pdf", 3.5, 8)
print(p1)
dev.off()

### Making bar plot for cell number from individual somite counts (Fig. 1c)

pd_cell_num_2 = pd[!is.na(pd$somite_count),] %>% group_by(somite_count) %>% tally() %>% rename(cell_num = n) %>% as.data.frame()
pd_cell_num_2$somite_count = factor(pd_cell_num_2$somite_count, levels = rev(names(somite_color_plate)))

p2 = pd_cell_num_2 %>%
    ggplot(aes(somite_count, cell_num, fill = somite_count)) + 
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_manual(values = somite_color_plate) + 
    scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    geom_text(aes(label = scales::comma(cell_num)), 
              hjust = -0.1,
              position = position_dodge(width = 1),
              inherit.aes = TRUE,
              size = 3) +
    labs(x = "", y = "Cell number") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
pdf("Cell_number_somite_counts.pdf", 3.5, 8)
print(p2)
dev.off()


###############################################################################
### Making 3D UMAP for individual major_trajectories (Supplementary Fig. 4) ###
###############################################################################

major_trajectory_list = names(major_trajectory_color_plate)

for(i in major_trajectory_list){
    print(i)
    
    if(sum(pd$global_celltype == i) > 300000){
        fig = pd %>%
            filter(global_celltype == i) %>%
            sample_n(300000) %>%
            plot_ly(x = ~sub_UMAP_3d_1, y = ~sub_UMAP_3d_2, z = ~sub_UMAP_3d_3, size=I(30), color = ~sub_celltype) %>% 
            layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                                yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                                zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2)))
    } else {
        fig = pd %>%
            filter(global_celltype == i) %>%
            plot_ly(x = ~sub_UMAP_3d_1, y = ~sub_UMAP_3d_2, z = ~sub_UMAP_3d_3, size=I(30), color = ~sub_celltype) %>% 
            layout(scene = list(xaxis=list(title = list(text ='UMAP_1', font = t1), tickfont = t2),
                                yaxis=list(title = list(text ='UMAP_2', font = t1), tickfont = t2),
                                zaxis=list(title = list(text ='UMAP_3', font = t1), tickfont = t2)))
    }
    
    saveWidget(fig, paste0(i, "_celltype_update.html"), selfcontained = FALSE, libdir = "tmp")
    
}



