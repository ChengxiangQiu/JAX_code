
###################################
### Section - 1, Basic analysis ###
###################################

################################################
### Data preprocessing and doublets removing ###
################################################

### Please note that we have provided all of the processed data, after preprocessing and removing doublets. 
### This script is only necessary if you want to reprocess the raw data or process your own data.

### Data from each individual sci-RNA-seq3 experiment was processed independently. 
### For each experiment, read alignment and gene count matrix generation was performed using the pipeline that we developed for sci-RNA-seq3 (https://github.com/JunyueC/sci-RNA-seq3_pipeline)

### The data generated in this study can be downloaded in raw and processed forms from the NCBI Gene Expression Omnibus under accession number GSE186069 and GSE228590. 


##################################################################################
### Step-1, read the output of pipeline and roughly removing low quality cells ###
##################################################################################

### After running the pipeline, you will obtain a profile named sci_summary.RData. 
### The file run_16_RT-sci-samplesheet.csv contains the RT barcode information for individual cells. 
### In this example, we used the run_16 data, but this script can be used for any experiment.

work_path = "./"

load("sci_summary.RData")
### It will load three profiles, gene_count, df_cell, and df_gene

rownames(df_gene)    = unlist(lapply(rownames(df_gene), function(x) strsplit(x,"[.]")[[1]][1]))
rownames(gene_count) = unlist(lapply(rownames(gene_count), function(x) strsplit(x,"[.]")[[1]][1]))
df_gene$gene_id      = unlist(lapply(as.vector(df_gene$gene_id), function(x) strsplit(x,"[.]")[[1]][1]))

### read RT_barcode to extract RT_group for each individual cell
RT_barcode = read.csv(paste0(work_path, "run_16_RT-sci-samplesheet.csv"), as.is=T, header=T)
names(RT_barcode) = c("RT_barcode_well", "RT_barcode_sequence", "RT_group", "genome")

df_cell$RT_barcode_sequence = str_sub(unlist(lapply(as.vector(df_cell$sample), function(x) strsplit(x,"[.]")[[1]][2])), -10, -1)
df_cell = df_cell %>%
    left_join(RT_barcode[,c("RT_barcode_sequence", "RT_group")], by = "RT_barcode_sequence")

day = rep(NA, nrow(df_cell))
day[df_cell$RT_group %in% c("E11.75_L4-01")] = "E11"
day[df_cell$RT_group %in% c("E11.5_L2-01")] = "E11.25"
day[df_cell$RT_group %in% c("E11.5_L5-06")] = "E11.5"
day[df_cell$RT_group %in% c("E11.75_L7-09")] = "E11.75"
df_cell$day = as.vector(day)
df_cell$batch = "run_16"

df_cell$UMI_count = Matrix::colSums(gene_count)
gene_count_copy = gene_count
gene_count_copy@x[gene_count_copy@x > 0] = 1
df_cell$gene_count = Matrix::colSums(gene_count_copy)

mouse_gene = read.table("mouse.v12.geneID.txt", header=T, sep="\t", as.is=T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene = mouse_gene[as.vector(df_gene$gene_id),]
df_gene$chr = as.vector(mouse_gene$chr)
names(df_gene) = c("gene_id", "gene_type", "gene_short_name", "chr")
gene_keep = df_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y"))

keep = df_cell$UMI_count >= 200 & 
    df_cell$gene_count >= 100 & 
    df_cell$unmatched_rate < 0.4 &
    df_cell$RT_group != "Barnyard"
df_cell = df_cell[keep,]
df_gene = df_gene[gene_keep,]
gene_count = gene_count[gene_keep, keep]
rownames(gene_count) = as.vector(df_gene$gene_id)
colnames(gene_count) = as.vector(df_cell$sample)
rownames(df_cell) = as.vector(df_cell$sample)

print(dim(gene_count))

saveRDS(gene_count, paste0(work_path, "gene_count.rds"))
saveRDS(df_cell, paste0(work_path, "df_cell.rds"))

print(nrow(df_cell))
print(median(df_cell$UMI_count))
print(median(df_cell$gene_count))
print(sum(is.na(df_cell$RT_group)))


######################################################
### Step-2, identifying doublets based on Scrublet ###
######################################################

### We performed three steps with the aim of exhaustively detecting and removing potential doublets. Of note, all these analyses were performed separately on data from each experiment.
### Please see more details in our Methods section

################################################################################
### Step-1, Cells with doublet scores over 0.2 were annotated as detected doublets. 

### We first randomly split the dataset into multiple subsets (six for most of the experiments) in order to reduce the time and memory requirements

n = 6
df_cell$split_batch = sample(c(1:n), nrow(df_cell), replace = T)

gene_count = readRDS(paste0(work_path, "gene_count.rds"))
df_cell = readRDS(paste0(work_path, "df_cell.rds"))
df_cell = df_cell[colnames(gene_count),]

for(i in 1:n){
    gene_count_i = gene_count[,df_cell$split_batch == i]
    df_cell_i = df_cell[df_cell$split_batch ==i,]
    
    print(sum(colnames(gene_count_i) != rownames(df_cell_i)))
    
    writeMM(t(gene_count_i), paste0(work_path, "gene_count_", i, ".mtx"))
    write.csv(df_cell_i, paste0(work_path, "df_cell_", i, ".csv"))
}

### python Run_Scrublet.py

df = NULL
for(i in 1:n){
    print(i)
    df_i = read.csv(paste0(work_path, "df_cell_", i, ".csv"), header=T, row.names=1, as.is=T)
    
    doublet_scores_observed_cells_i = read.csv(paste0("doublet_scores_observed_cells_", i, ".csv"), header=F)
    df_i$doublet_score = as.vector(doublet_scores_observed_cells_i$V1)
    df = rbind(df, df_i)
}

df = df[rownames(df_cell),]
print(sum(rownames(df) != rownames(df_cell)))

df_cell$doublet_score = as.vector(df$doublet_score)
df_cell$detected_doublets = df_cell$doublet_score > 0.2

################################################################################
### Step-2, Subclusters with a detected doublet ratio (by Scrublet) over 15% were annotated as doublet-derived subclusters.

### python Detecting_doublets_by_subclustering.py

global = read.csv(paste0(work_path, "doublet_cluster/global.csv"), header=T)
main_cluster_list = sort(as.vector(unique(global$louvain)))

res = NULL

for(i in 1:length(main_cluster_list)){
    print(paste0(i, "/", length(main_cluster_list)))
    dat = read.csv(paste0(work_path, "doublet_cluster/adata.obs.louvain_", (i-1), ".csv"), header=T)
    print(nrow(dat))
    dat$louvain = as.vector(paste0("cluster_", dat$louvain))
    dat = dat %>%
        left_join(df_cell[,c("sample", "detected_doublets", "doublet_score")], by = "sample")
    
    tmp2 = dat %>%
        group_by(louvain) %>%
        tally() %>%
        dplyr::rename(n_sum = n)
    
    tmp1 = dat %>%
        filter(detected_doublets == "TRUE") %>%
        group_by(louvain) %>%
        tally() %>%
        left_join(tmp2, by = "louvain") %>%
        mutate(frac = n/n_sum) %>%
        filter(frac > 0.15)
    
    dat$doublet_cluster = dat$louvain %in% as.vector(tmp1$louvain) 
    
    dat[,!colnames(dat) %in% c("umap_1", "umap_2")]
    dat$main_louvain = (i-1)
    
    res = rbind(res, dat)
}

rownames(res) = as.vector(res$sample)
res = res[rownames(df_cell),]
df_cell$doublet_cluster = res$doublet_cluster

saveRDS(df_cell, paste0(work_path, "df_cell.rds"))

################################################################################
### Step-3, we found that the above Scrublet and iterative clustering based approach has difficulty identifying doublets in clusters derived from rare cell types (e.g. clusters comprising less than 1% of the total cell population), so we applied a third step to further detect and remove doublets.
### In this step, subclusters that showed low expression of target cell cluster-specific markers and enriched expression of non-target cell cluster-specific markers were identified as doublet-driven clusters. 

### Of note, this step was mainly done by Monocle/3-alpha
### http://cole-trapnell-lab.github.io/monocle-release/monocle3/

suppressMessages(library(dplyr))
suppressMessages(library(reticulate))
suppressMessages(library(monocle))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
suppressMessages(library(VGAM))
suppressMessages(library(ggrastr))
library(reticulate)
import("louvain")
print(packageVersion('monocle'))
### monocle 2.99

work_path = "./"

count = readRDS(paste0(work_path, "gene_count.rds"))
pd = readRDS(paste0(work_path, "df_cell.rds"))
fd = read.csv(paste0(work_path, "df_gene_all.csv"), header=T, row.names=1, as.is=T)

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
fd = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & 
            (!fd$chr %in% c('chrX', 'chrY')),]

### cells labeled as doublets (by Scrublet) or from doublet-derived subclusters were filtered out 
pd = pd[!(pd$detected_doublets | pd$doublet_cluster),]

count = count[rownames(fd), rownames(pd)]

### Genes expressed in less than 10 cells and cells expressing less than 100 genes were further filtered out
min.features = 100
min.cells = 10

# filter genes on the number of cells expressing
if (min.cells > 0) {
    num.cells <- Matrix::rowSums(x = count > 0)
    count <- count[which(x = num.cells >= min.cells), ]
}

# filer cells on the number of genes expressing
if (min.features > 0) {
    nfeatures <- Matrix::colSums(x = count > 0)
    count <- count[, which(x = nfeatures >= min.features)]
}

sum(!colnames(count) %in% rownames(pd))
sum(!rownames(count) %in% rownames(fd))

pd = pd[colnames(count),]
fd = fd[rownames(count),]

pData = new("AnnotatedDataFrame",data=pd) 
fData = new("AnnotatedDataFrame",data=fd)

cds <- newCellDataSet(count, 
                      phenoData = pData,
                      featureData =fData,
                      expressionFamily = negbinomial.size())

saveRDS(cds, paste0(work_path, "monocle_cds_alpha_1.rds"))

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table = dispersionTable(cds)
disp_table = disp_table %>% 
    mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 5000)$gene_id)

cds = setOrderingFilter(cds, top_subset_genes)
cds <- preprocessCDS(cds,  method = 'PCA',
                     norm_method = 'log',
                     num_dim = 50,
                     verbose = T)

cds <- reduceDimension(cds, 
                       max_components = 2,
                       reduction_method = 'UMAP',
                       metric="cosine",
                       min_dist = 0.1,
                       n_neighbors = 50,
                       verbose = T)

cds <- clusterCells(cds,
                    method = 'louvain',
                    res = 1e-6,
                    louvain_iter = 1,
                    verbose = T)

cds <- partitionCells(cds)

saveRDS(cds, paste0(work_path, "monocle_cds_alpha_2.rds"))

jpeg("louvain_component.jpeg", width=10, height=8, units = 'in', res=300)
plot_cell_clusters(cds,
                   color_by = 'louvain_component',
                   cell_size = 0.1,
                   show_group_id = T) +
    theme(legend.text=element_text(size=6)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
dev.off()

saveRDS(list(count = exprs(cds),
             pd = data.frame(pData(cds)),
             fd = data.frame(fData(cds))), paste0(work_path, "dat_2.rds"))

### We took the cell clusters identified by Monocle/3, downsampled each cluster to 2,500 cells, and computed differentially expressed genes across cell clusters with the top_markers function of Monocle/3 (reference_cells=1000).
### We selected a gene set combining the top ten gene markers for each cell cluster (filtering out genes with fraction_expressing < 0.1 and then ordering by pseudo_R2).

library(monocle3) ### Of note, this is done by Monocle/3, rather than Monocle/3-alpha

work_path = "./"

dat = readRDS(paste0(work_path, "dat_2.rds"))
exp = dat[['count']]
pd = dat[['pd']]
fd = dat[['fd']]

pd$myCluster = paste0('cluster_', pd$louvain_component)

myCluster_table = table(pd$myCluster)
small_cluster = names(myCluster_table)[myCluster_table < 2500]

pd_sub_1 = pd %>%
    filter(!myCluster %in% small_cluster) %>%
    group_by(myCluster) %>%
    sample_n(2500) %>%
    as.data.frame()
pd_sub_2 = pd %>%
    filter(myCluster %in% small_cluster) %>%
    as.data.frame()
pd_sub = rbind(pd_sub_1, pd_sub_2)
rownames(pd_sub) = as.vector(pd_sub$sample)
exp_sub = exp[,rownames(pd_sub)]

cds <- new_cell_data_set(exp_sub,
                         cell_metadata = pd_sub,
                         gene_metadata = fd)

cds <- preprocess_cds(cds, num_dim = 50)
cds = reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

marker_test_res <- top_markers(cds, 
                               group_cells_by="myCluster", 
                               reference_cells=1000,
                               cores = 8)

markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(10, pseudo_R2)

saveRDS(markers, file=paste0(work_path, "marker_top10.rds"))

### Cells from each main cell cluster were subjected to dimensionality reduction by PCA (10 components) on the selected set of top cluster-specific gene markers.
### Each cell cluster was further reduced to 2D using UMAP (max_components = 2, n_neighbors = 50, min_dist = 0.1, metric = 'cosine'). 
### The cells within each cluster were further subclustered using the Louvain algorithm implemented in Monocle/3 (res = 1e-04 for most clustering analysis). 

### mkdir doublet_cluster_2

suppressMessages(library(dplyr))
suppressMessages(library(reticulate))
suppressMessages(library(monocle))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
suppressMessages(library(VGAM))
suppressMessages(library(ggrastr))
library(reticulate)
import("louvain")
print(packageVersion('monocle'))
### monocle 2.99

work_path = "./"

dat = readRDS(paste0(work_path, "dat_2.rds"))
exp = dat[['count']]
pd = dat[['pd']]
fd = dat[['fd']]

pd$myCluster = paste0('cluster_', pd$louvain_component)
cluster_list = as.vector(unique(pd$myCluster))
cluster_list = paste0("cluster_", 1:length(cluster_list))

markers = readRDS(paste0(work_path, "/marker_top10.rds"))

for(i in 1:length(cluster_list)){
    
    cluster_i = cluster_list[i]
    print(cluster_i)
    
    index = pd$myCluster == cluster_i
    
    pData = new("AnnotatedDataFrame",data=pd[index,]) 
    fData = new("AnnotatedDataFrame",data=fd)
    
    cds_subset <- newCellDataSet(exp[,index], 
                                 phenoData = pData,
                                 featureData =fData,
                                 expressionFamily = negbinomial.size())
    
    print(dim(cds_subset))
    
    DelayedArray:::set_verbose_block_processing(TRUE)
    options(DelayedArray.block.size=1000e6)
    cds_subset <- estimateSizeFactors(cds_subset)
    cds_subset <- estimateDispersions(cds_subset)
    
    cds_subset = setOrderingFilter(cds_subset, as.vector(markers$gene_id))
    cds_subset <- preprocessCDS(cds_subset,  
                                method = 'PCA',
                                norm_method = 'log',
                                num_dim = 10,
                                verbose = T)
    
    cds_subset <- reduceDimension(cds_subset, 
                                  max_components = 2,
                                  reduction_method = 'UMAP',
                                  metric="cosine",
                                  min_dist = 0.1,
                                  n_neighbors = 50,
                                  verbose = T)
    
    for(j in 1:length(cluster_list)){
        cluster_j = cluster_list[j]; print(cluster_j)
        markers_sub = markers %>% filter(cell_group == cluster_j)
        jpeg(paste0(work_path, 'doublet_cluster_2/', cluster_i, "_", cluster_j, '.jpeg'), width=10, height=8, units = 'in', res=300)
        print(plot_cell_clusters(cds_subset,
                                 markers = as.character(markers_sub$gene_short_name),
                                 cell_size = 0.5))
        dev.off()
    }
    
    p = plot_cell_clusters(cds_subset,
                           color_by = 'doublet_score',
                           cell_size = 0.8) + scale_color_viridis() 
    
    jpeg(paste0(work_path, 'doublet_cluster_2/', cluster_i, '_doublet_score.jpeg'), width=10, height=10, units = 'in', res=300)
    print(p)
    dev.off()
    
    saveRDS(cds_subset, paste0(work_path, 'doublet_cluster_2/', cluster_i, '.rds'))
    
    ### testing different clustering parameters
    if(nrow(pData(cds_subset)) > 50000){
        resolution_list = c(1e-3, 5e-4, 1e-4)
    } else if(nrow(pData(cds_subset)) > 10000){
        resolution_list = c(5e-3, 1e-3, 1e-4)
    } else if (nrow(pData(cds_subset)) < 500){
        resolution_list = c(0.5, 0.2, 0.1)
    } else {
        resolution_list = c(1e-2, 5e-3, 1e-3)
    }
    
    for(xx in 1:length(resolution_list)){
        
        res_para = resolution_list[xx]
        
        cds_subset <- clusterCells(cds_subset,
                                   method = 'louvain',
                                   res = res_para,
                                   louvain_iter = 1,
                                   verbose = T)
        
        if(xx == 1){
            cds_subset$cluster_1 = cds_subset$Cluster
        } else if(xx == 2){
            cds_subset$cluster_2 = cds_subset$Cluster
        } else {
            cds_subset$cluster_3 = cds_subset$Cluster
        }
        
        saveRDS(data.frame(pData(cds_subset)), paste0(work_path, 'doublet_cluster_2/',  cluster_i, '_pd.rds'))
        
        p = plot_cell_clusters(cds_subset,
                               color_by = 'Cluster',
                               cell_size = 1,
                               show_group_id = T) + theme(legend.position="none")
        
        jpeg(paste0(work_path, 'doublet_cluster_2/', cluster_i, "_res_", xx, '.jpeg'), width=10, height=10, units = 'in', res=300)
        print(p)
        dev.off()
        
        pd = data.frame(pData(cds_subset))
        pd$Cluster = factor(paste0("cluster_", pd$Cluster))
        p = ggplot(pd, aes(factor(Cluster), doublet_score)) +
            geom_boxplot() +
            coord_flip() 
        jpeg(paste0(work_path, 'doublet_cluster_2/', cluster_i, "_boxplot_res_", xx, '.jpeg'), width=10, height=15, units = 'in', res=300)
        print(p)
        dev.off()
    }
    
}

### Subclusters showing low expression of the differentially expressed genes identified in step 5 and enriched expression of specific markers of a different cluster were annotated as doublet-derived subclusters and filtered out. This procedure eliminated 0.5-13.2% of the cells in each experiment.

pd_all = NULL

for(i in 1:length(cluster_list)){
    
    print(i)
    
    pd_tmp = readRDS(paste0(work_path, "doublet_cluster_2/cluster_", i, "_pd.rds"))
    
    pd_tmp$cluster_1 = as.vector(paste0("cluster_", pd_tmp$cluster_1))
    pd_tmp$cluster_2 = as.vector(paste0("cluster_", pd_tmp$cluster_2))
    pd_tmp$cluster_3 = as.vector(paste0("cluster_", pd_tmp$cluster_3))
    
    if(i == 1){
        pd_tmp$doublet_DEG = pd_tmp$cluster_1 %in% paste0("cluster_", c(21))
    }
    
###    ...
### manually checking each individual cluster to identify doublets-driven subclusters
    
    pd_tmp = pd_tmp[,c("sample", "doublet_DEG")]
    pd_all = rbind(pd_all, pd_tmp)
    
}

saveRDS(pd_all, paste0(work_path, "doublet_cluster_2/res_doubelt_DEG.rds"))


##################################################################################
### Step-3, read the output of pipeline and roughly removing low quality cells ###
##################################################################################

count = readRDS(paste0(work_path, "gene_count.rds"))
pd = readRDS(paste0(work_path, "df_cell.rds"))
fd = read.csv(paste0(work_path, "df_gene_all.csv"), header=T, row.names=1, as.is=T)

res_doublet_DEG = readRDS(paste0("doublet_cluster_2/res_doubelt_DEG.rds"))
res_doublet_DEG = res_doublet_DEG[res_doublet_DEG$doublet_DEG,]
pd$doublet_DEG = rownames(pd) %in% as.vector(res_doublet_DEG$sample)
pd$log2_umi = log2(pd$UMI_count)
pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)

### calculate MT_pct and Ribo_pct per cell
gene = fd
MT_gene = as.vector(gene[grep("^mt-",gene$gene_short_name),]$gene_id)
Rpl_gene = as.vector(gene[grep("^Rpl",gene$gene_short_name),]$gene_id)
Mrpl_gene = as.vector(gene[grep("^Mrpl",gene$gene_short_name),]$gene_id)
Rps_gene = as.vector(gene[grep("^Rps",gene$gene_short_name),]$gene_id)
Mrps_gene = as.vector(gene[grep("^Mrps",gene$gene_short_name),]$gene_id)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd$MT_pct = 100 * Matrix::colSums(count[gene$gene_id %in% MT_gene, ])/Matrix::colSums(count)
pd$RIBO_pct = 100 * Matrix::colSums(count[gene$gene_id %in% RIBO_gene, ])/Matrix::colSums(count)

### removing doublet cells that are detected by the above three strategies.
pd = pd[!(pd$detected_doublets | pd$doublet_cluster | pd$doublet_DEG),]

### first we removed cells with exon% > 85%
### then we identified the peak of the histgram and then set the log2_UMI cutoff as 
### [mean(x) - sd(x), mean(x) + 2*sd(x)]
### for a handful of experiments, we manually set the low boundary of cutoff

x_tmp = pd$log2_umi[pd$EXON_pct <= 85]
x1 = mean(x_tmp) - sd(x_tmp)
x2 = mean(x_tmp) + 2*sd(x_tmp)

pd = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <= 85,]

### save the final dataset with cells after removing doublets and low-quality cells
saveRDS(pd, paste0(work_path, "pd.rds"))














