
##################################
### Section - 5, Neuroectoderm ###
##################################

#################################################
### 2D UMAP of subclustering on early neurons ###
#################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "Neurons"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Fig. 5e

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = neurons_sub_clustering), size=0.3) +
    theme_void() +
    scale_color_manual(values=neuroectoderm_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)

### Extended Data Fig. 10c

day_list = names(day_color_plate)
pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.3) +
    theme_void() +
    scale_color_manual(values=neuron_day_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".day.2D_UMAP.png"), width = 6, height = 6, dpi = 300)


####################################################
### 2D UMAP of Intermediate neuronal progenitors ###
####################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "INP"; print(example_i)

pd = readRDS(paste0(work_path, example_i, "_adata_scale.obs.rds"))

### Extended Data Fig. 10a

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.3) +
    theme_void() +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".2D_UMAP.png"), width = 6, height = 6, dpi = 300)

day_list = names(day_color_plate)
pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = day), size=0.3) +
    theme_void() +
    scale_color_manual(values=neuron_day_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, example_i, ".day.2D_UMAP.png"), width = 6, height = 6, dpi = 300)


####################################################################################################
### Composition of embryos from each 6-hr bin by intermediate neuronal progenitor and CNS neuron ###
####################################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "df_cell.rds"))
day_list = names(day_color_plate)
pd$day = factor(pd$day, levels = day_list[day_list %in% pd$day])

pd_1 = pd[pd$major_trajectory == "CNS_neurons",]
pd_2 = pd[pd$major_trajectory == "Intermediate_neuronal_progenitors",]

x1 = pd_1 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x2 = pd_2 %>% group_by(day) %>% tally() %>%
    left_join(pd %>% group_by(day) %>% tally() %>% rename(total_n = n), by = "day") %>%
    mutate(frac = 100*n/total_n)

x = x1 %>% select(day, frac) %>% rename(direct_frac = frac) %>% left_join(x2 %>% select(day, frac) %>% rename(indirect_frac = frac), by = "day")
x$indirect_frac[is.na(x$indirect_frac)] = 0
x = data.frame(day = rep(x$day, 2),
               frac = c(x$direct_frac, x$indirect_frac),
               major_trajectory = rep(c("CNS_neurons","Intermediate_neuronal_progenitors"), each = nrow(x)))

### Fig. 4d

p = x %>% 
    ggplot(aes(x=day, y=frac, fill = day)) + 
    geom_bar(stat='identity') + facet_grid(rows = vars(major_trajectory)) + 
    labs(x='',y='% of cells') +
    scale_fill_manual(values=day_color_plate) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(color="black"))


#############################################################
### For each interneuron, what are the top expressed TFs? ### 
#############################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

dat = readRDS(paste0(work_path, "Neurons_heatmap_dat.rds"))

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "Neurons_heatmap.pdf"), 8, 5)
heatmap.2(as.matrix(t(dat)), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()








