
#####################################
### Section - 2, Posterior embryo ###
#####################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "pd_somites.rds"))
### n = 104,671 nuclei

celltype_color_plate = c("#54c15f", "#c34fb7", "#91b737", "#7b63d0", "#c9a63c",
                         "#7081ca", "#478734", "#da3f78", "#56c09e", "#cf4a35",
                         "#999999", "#4eacd7", "#e18f4f", "#be75b4", "#a0b46c",
                         "#a0445d", "#37845f", "#df7c82", "#72732b", "#a06432")
names(celltype_color_plate) = x

### Supplementary Figure 7b
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = anno), size=0.6) +
    theme_void() +
    scale_color_manual(values=color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "somites.anno.png"), width = 6, height = 6, dpi = 300)

somite_color_plate = c("#440154", "#482475", "#414487", "#355f8d",
                       "#2a788e", "#21918c", "#22a884", "#44bf70",
                       "#7ad151", "#bddf26", "#fde725")
names(somite_color_plate) = paste0(c(8,9,10,11,12,13,14,16,17,20,21), " somites")

### Supplementary Figure 7c
p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.8) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.6) +
    theme_void() +
    scale_color_manual(values=somite_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "somites.somite_count.png"), width = 6, height = 6, dpi = 300)

### Supplementary Figure 7a
pd$embryo_id = factor(pd$embryo_id, levels = rev(names(table(pd$embryo_id))))
p1 = pd %>%
    group_by(embryo_id, somite_count) %>% tally() %>% rename(cell_num = n) %>%
    ggplot(aes(embryo_id, cell_num, fill = somite_count)) + 
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_manual(values = somite_color_plate) + 
    geom_text(aes(label = scales::comma(cell_num)), 
              hjust = -0.1,
              position = position_dodge(width = 1),
              inherit.aes = TRUE,
              size = 5) +
    labs(x = "", y = "Cell number") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
pdf(paste0(work_path, "cell_num.pdf"), 5, 8)
print(p1)
dev.off()


#######################
### Focusing on NMP ###
#######################


pd_NMP = read.csv(paste0(work_path, "adata_somites_NMP.obs.csv"), row.names=1, as.is=T)

somite_color_plate = c("#440154", "#482475", "#414487", "#355f8d",
                       "#2a788e", "#21918c", "#22a884", "#44bf70",
                       "#7ad151", "#bddf26", "#fde725")
names(somite_color_plate) = paste0(c(8,9,10,11,12,13,14,16,17,20,21), " somites")

pd_NMP$somite_count = factor(pd_NMP$somite_count, levels = paste0(c(8,9,10,11,12,13,14,16,17,20,21), " somites"))
p = pd_NMP %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=1) +
    theme_void() +
    scale_color_manual(values=somite_color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "NMP.somite_count.png"), width = 6, height = 6, dpi = 300)


p = pd_NMP %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = anno), size=1) +
    theme_void() +
    scale_color_manual(values=color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "NMP.anno.png"), width = 6, height = 6, dpi = 300)

pd = data.frame(pData(cds))
pd$somite_count = factor(pd$somite_count, levels = paste0(c(8,9,10,11,12,13,14,16,17,20,21), " somites"))
pd$Cdx1 = as.vector(exprs(cds)["ENSMUSG00000024619",])
pd$Hoxa10 = as.vector(exprs(cds)["ENSMUSG00000000938",])
pd$T = as.vector(exprs(cds)["ENSMUSG00000062327",])
pd$Meis1 = as.vector(exprs(cds)["ENSMUSG00000020160",])

df = pd %>% filter(Cdx1 != 0) %>% group_by(somite_count, .drop = FALSE) %>% tally() %>%
    left_join(pd %>% group_by(somite_count) %>% tally() %>% rename(total_n = n)) %>%
    mutate(percent = 100*n/total_n)

p1 <-ggplot(data=df, aes(x=somite_count, y=percent, fill = somite_count)) +
    geom_bar(stat="identity") + labs(x="",y="% of cells expressed Cdx1") +
    scale_fill_manual(values=somite_color_plate) + theme_classic(base_size = 10) + theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90), axis.text.y = element_text(color="black")) 

df = pd %>% filter(Hoxa10 != 0) %>% group_by(somite_count, .drop = FALSE) %>% tally() %>%
    left_join(pd %>% group_by(somite_count) %>% tally() %>% rename(total_n = n)) %>%
    mutate(percent = 100*n/total_n)

p2 <-ggplot(data=df, aes(x=somite_count, y=percent, fill = somite_count)) +
    geom_bar(stat="identity") + labs(x="",y="% of cells expressed Hoxa10") +
    scale_fill_manual(values=somite_color_plate) + theme_classic(base_size = 10) + theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90), axis.text.y = element_text(color="black")) 

library(gridExtra) 
pdf(paste0(work_path, "NMP_Cdx1_Hoxa10.pdf"), 4, 6)
grid.arrange(p1, p2, nrow=2, ncol=1) 
dev.off()



