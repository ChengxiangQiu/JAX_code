
#####################################
### Section - 2, Posterior embryo ###
#####################################

#### Supplementary Figure 9

### First, we downsampled the dataset to ~1M cells using geosketch (Hie et al. 2019) 
### and performed k-means clustering to ensure that each cluster contained roughly 500 cells. 
### Second, we aggregated UMI counts for cells within each cluster to generate 2,289 meta-cells, 
### and normalized the UMI counts for each meta-cell followed by log2-transformation. 
### Third, we performed Pearson correlation between each protein-coding gene and Npm1, 
### and selected genes with correlation coefficients > 0.6 (738 genes, ~3% of the 
### total protein coding genes). A gene set enrichment analysis suggests that the module 
### is associated with RNP complexes (corrected p-value = 1.4e-105), cytoplasmic translation 
### (corrected p-value = 2.8e-90), and ribosomal proteins (corrected p-value = 7.4e-71). 
### Finally, we summed the normalized UMI counts of these genes to calculate a Npm1 signature 
### for individual cells. 

### First, performing geosketch to downsample the original dataset to ~1M cells
### python run_geosketch.py

### After running geosktech, performing dimension reduction on the subset using monocle/3

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd_all = readRDS(paste0(work_path, "df_cell.rds"))

pd = read.csv(paste0(work_path, "adata_scale.obs.csv"))
pd_geosketch = read.csv(paste0(work_path, "adata_scale_geosketch_downsample.csv"))
pd = pd[as.vector(pd_geosketch$V1),]

pd_sub = pd_all[pd_all$cell_id %in% as.vector(pd$cell_id),]
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding')) & mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

gene_count = doExtractData(pd_sub, mouse_gene_sub)

cds = new_cell_data_set(gene_count,
                        cell_metadata = pd_sub,
                        gene_metadata = mouse_gene[rownames(gene_count), c("gene_ID", "gene_short_name")])
cds = preprocess_cds(cds)
cds = reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", umap.min_dist = 0.1, umap.n_neighbors = 50)
print(dim(cds))
saveRDS(cds, paste0(work_path, "cds_geosketch.rds"))


### Performing kmeans clustering to generate meta-cells

exp = exprs(cds)
cell_size = 500

ks = floor(ncol(cds)/cell_size) + 1
emb = as.matrix(reducedDims(cds)$PCA[,c(1:30)])
meta_group = paste0("M", kmeans(emb, centers = ks, nstart = 25)$cluster)
saveRDS(meta_group, paste0(work_path, "meta_group_", cell_size, ".rds"))

meta_group = readRDS(paste0(work_path, "meta_group_", cell_size, ".rds"))

start_time = Sys.time()

gene_count = NULL
for(i in 1:length(table(meta_group))){
    print(paste0("M",i))
    gene_count = cbind(gene_count, 
                       Matrix::rowSums(exp[, meta_group == paste0("M",i), drop=FALSE]))
}
colnames(gene_count) = paste0("M", 1:length(table(meta_group)))

end_time = Sys.time()

print(end_time - start_time)

saveRDS(gene_count, paste0(work_path, "meta_group_", cell_size, ".gene_count.rds"))

### Identifying genes which are correlated with Npm1

cds <- new_cell_data_set(gene_count,
                         gene_metadata = mouse_gene[rownames(gene_count),])

gene_count_norm = t(t(gene_count)/as.vector(cds$Size_Factor))
gene_count_norm_log2 = log2(gene_count_norm + 1)
gene_count_norm_log2_trans = t(gene_count_norm_log2)

gene_count_norm_log2_trans_sum = apply(gene_count_norm_log2_trans, 2, sum)

cor_res = cor(as.matrix(gene_count_norm_log2_trans[,c("ENSMUSG00000057113"), drop = F]),
              as.matrix(gene_count_norm_log2_trans[,gene_count_norm_log2_trans_sum != 0]))
cor_res = data.frame(gene_ID = colnames(cor_res),
                     corr = as.vector(cor_res), stringsAsFactors = F)
cor_res = cor_res %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")]) %>% as.data.frame()

saveRDS(cor_res, paste0(work_path, "Npm1_corr_genes.rds"))

### calculating Npm1 signature

cds = readRDS(paste0(work_path, "cds_geosketch.rds"))
cor_res = readRDS(paste0(work_path, "Npm1_corr_genes.rds"))
gene_count = exprs(cds)
gene_count_sub = gene_count[rownames(gene_count) %in% as.vector(cor_res$gene_ID[cor_res$corr > 0.6]),]
gene_count_norm = t(t(gene_count_sub)/as.vector(cds$Size_Factor))

pd = data.frame(pData(cds))
pd$score = Matrix::colSums(gene_count_norm)

sample_sheet = read.table(paste0(work_path, "sample_sheet.txt"), sep="\t", header = T, as.is=T)

pd = pd %>% left_join(sample_sheet[,c("embryo_id", "Plug_Date", "Date_Harvested", "litter_ID", "shipment_ID")]) %>% as.data.frame()
pd$shipment_ID = unlist(lapply(as.vector(pd$shipment_ID), function(x) strsplit(x,"[_]")[[1]][1]))
pd$shipment_ID = factor(pd$shipment_ID, levels = paste0("Shipment", 4:1))

x = as.vector(pd$sequencing_batch)
x[pd$sequencing_batch == "run_23"] = "run_23_A"
x[pd$sequencing_batch == "run_27"] = "run_23_B"
pd$sequencing_batch = factor(x, levels = rev(paste0("run_", c(4, 13:22, "23_A", "23_B", 24:26))))
table(pd$sequencing_batch)

x = as.vector(pd$litter_ID)
x[pd$litter_ID == "12.5-L5"] = 'E12.5_L5'
x[pd$litter_ID == "E10.0-L1"] = 'E10.0_L1'
x[pd$litter_ID == "E10.0-L4"] = 'E10.0_L4'
x[pd$litter_ID == "E12.5-L3"] = 'E12.5_L3'
y = data.frame(id = names(table(x)), stringsAsFactors = F)
y$tmp1 = unlist(lapply(as.vector(y$id), function(x) strsplit(x,"[_]")[[1]][1])) 
day_list_x = day_list
day_list_x[day_list == "E14.333"] = "E14.5"
y$tmp1 = factor(y$tmp1, levels = day_list_x)
y$tmp2 = unlist(lapply(as.vector(y$id), function(x) strsplit(x,"[_]")[[1]][2])) 
y = y %>% arrange(tmp1, tmp2)
pd$litter_ID = factor(x, levels = rev(as.vector(y$id)))
table(pd$litter_ID)

y = data.frame(id = names(table(pd$Date_Harvested)), stringsAsFactors = F)
y$tmp1 = unlist(lapply(as.vector(y$id), function(x) strsplit(x,"[/]")[[1]][1])) 
y$tmp2 = unlist(lapply(as.vector(y$id), function(x) strsplit(x,"[/]")[[1]][2])) 
y$tmp3 = unlist(lapply(as.vector(y$id), function(x) strsplit(x,"[/]")[[1]][3])) 
y = y %>% arrange(tmp3, tmp1, tmp2)
pd$Date_Harvested = factor(as.vector(pd$Date_Harvested), levels = rev(as.vector(y$id)))
table(pd$Date_Harvested)

x = sort(table(pd$major_trajectory), decreasing = T)
selected_major_trajectory = names(x)[1:10]
pd$major_trajectory = factor(pd$major_trajectory, levels = names(x))
table(pd$major_trajectory)

p = ggplot(pd, aes(sequencing_batch, score, fill = sequencing_batch)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 13) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_batch.pdf"), 4, 6)
print(p)
dev.off()

p = ggplot(pd %>% filter(major_trajectory %in% selected_major_trajectory), aes(sequencing_batch, score, fill = sequencing_batch)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(cols = vars(major_trajectory), scale = "free") +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 8) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_batch_x.pdf"), 8, 4)
print(p)
dev.off()


p = ggplot(pd, aes(Date_Harvested, score, fill = Date_Harvested)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 13) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_Date_Harvested.pdf"), 4, 6)
print(p)
dev.off()

p = ggplot(pd %>% filter(major_trajectory %in% selected_major_trajectory), aes(Date_Harvested, score, fill = Date_Harvested)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(cols = vars(major_trajectory), scale = "free") +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 8) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_Date_Harvested_x.pdf"), 8, 4)
print(p)
dev.off()


p = ggplot(pd, aes(litter_ID, score, fill = litter_ID)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 13) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_litter_ID.pdf"), 4, 6)
print(p)
dev.off()

p = ggplot(pd %>% filter(major_trajectory %in% selected_major_trajectory), aes(litter_ID, score, fill = litter_ID)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(cols = vars(major_trajectory), scale = "free") +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 8) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_litter_ID_x.pdf"), 8, 4)
print(p)
dev.off()


p = ggplot(pd, aes(shipment_ID, score, fill = shipment_ID)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 13) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_shipment_ID.pdf"), 4, 6)
print(p)
dev.off()

p = ggplot(pd %>% filter(major_trajectory %in% selected_major_trajectory), aes(shipment_ID, score, fill = shipment_ID)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(cols = vars(major_trajectory), scale = "free") +
    labs(x = "", y = "Npm1 signature") +
    theme_classic(base_size = 8) +
    theme(legend.position="none") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylim(0, 450)
pdf(paste0("Npm1_signature_shipment_ID_x.pdf"), 8, 4)
print(p)
dev.off()





