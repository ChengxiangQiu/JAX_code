
###########################################
### Profiles that are used for analysis ###
###########################################

### All those data are provided from
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax/download/other/
### Please let me know if any questions (cxqiu@uw.edu)

### Section_1_basic_analysis
### mouse.v12.geneID.txt                                  ### gene annotation based on GENCODE.M12
### df_gene_all.csv                                       ### gene annotation in a simple format, including 49,585 genes on chr1-19, chrM, chrX, and chrY
### run_16_RT-sci-samplesheet.csv                         ### RT barcode information for individual cells in run_16
### df_cell.rds                                           ### meta information of individual cells
### embryo_cds.rds                                        ### pseudobulk dataset of individual embryos, in monocle/v3 format
### cell_num_prediction.rds                               ### predicted cell number for individual embryos

### Section_2_posterior_embryo
### posterior_embryo_adata_scale.obs.csv                  ### meta information of individual cells for posterior embryo subset
### posterior_embryo_adata_scale.NMP_Mesoderm.obs.csv     ### meta information of individual cells for posterior embryo subset (NMP)
### posterior_embryo_adata_scale.Notochord.obs.csv        ### meta information of individual cells for posterior embryo subset (Notochord)
### posterior_embryo_adata_scale.Gut.obs.csv              ### meta information of individual cells for posterior embryo subset (Gut)
### pd_somites.rds                                        ### meta information of individual cells for the somitogenesis validation dataset
### adata_somites_NMP.obs.csv                             ### meta information of individual cells for the NMP in the validation dataset

### Section_3_kidney_mesenchyme
### Renal_adata_scale.obs.rds                             ### meta information of renal subset
### Renal_CDI_adata_scale.obs.rds                         ### meta information of a subset of renal cells
### LPM_adata_scale.obs.rds                               ### meta information of lateral plate & intermediate mesoderm
### Mosta_file_list.txt                                   ### a list of file names that downloaded from Mosta database

### Section_4_eye
### Eye_adata_scale.obs.rds                               ### meta information of eye subset
### Eye_RGC_adata_scale.obs.rds                           ### meta information of RGCs subset
### Eye_RGC_heatmap_dat.rds                               ### data used for plot the heatmap in Fig.4e
### Eye_early_adata_scale.obs.rds                         ### meta information of early eay subset (<=E12.5)
### Eye_iris_adata_scale.obs.rds                          ### meta information of iris related cells

### Section_5_neuroectoderm
### Neuroectoderm_backbone_adata_scale.obs.rds            ### meta information of neuroectoderm backbone cells (<E13.0)
### Neuroectoderm_derivative_adata_scale.obs.rds          ### meta information of neuroectoderm and derivatives (<E13.0) 
### Neuroectoderm_derivative_adata_scale.PCs.rds          ### PC features of neuroectoderm and derivatives (<E13.0) 
### Neurons_adata_scale.obs.rds                           ### meta information of subclustering result of early neurons (Fig.5e)
### INP_adata_scale.obs.rds                               ### meta information of subclustering result of intermediate neuronal progenitors major cell cluster
### Neurons_heatmap_dat.rds                               ### data used for plot the heatmap in Fig.5f
### Astrocytes_adata_scale.obs.rds                        ### meta information of astrocytes (Fig.S13, <E13.0)
### Neurons_interneurons_top_TFs.csv                      ### top key TFs identified for each spinal interneurons (Fig.S14b)
### Neurons_interneurons_overlap_TFs.csv                  ### top TFs that highly expressed in the progenitors of each spinal interneurons (Fig.S14c)

### Section_6_development_tree

### Section_7_birth_series





################################################
### Packges that are needed for the analysis ###
################################################

library(monocle3)
library(Seurat)
library(lattice)
library(FNN)
library(dplyr)
library(ggplot2)
library(viridis)
library(Matrix)
library(stringr)
library(plotly)
library(htmlwidgets)
library(gridExtra) 
library(reshape2)
library(gplots)
library(RColorBrewer)
library(rdist)
library(ggridges)

mouse_gene = read.table("mouse.v12.geneID.txt", header = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

################################################################################
### Function: estimating gender based on ratio of gene expression on X and Y ###
################################################################################

estimateSex <- function(count){
    Y_gene = mouse_gene$gene_ID[mouse_gene$chr == "chrY" & mouse_gene$gene_short_name != "Erdr1"]
    
    ### Xist
    res = data.frame(Y_exp = Matrix::colSums(count[rownames(count) %in% Y_gene,]),
                     X_exp = as.vector(count["ENSMUSG00000086503",]))
    
    estimated_sex = rep("Unsure", nrow(res))
    estimated_sex[res$X_exp > res$Y_exp] = "Female"
    estimated_sex[res$X_exp < res$Y_exp] = "Male"
    res$estimated_sex = as.vector(estimated_sex)
    
    return(res)
}


#########################################################################
### Function: transition data object between monocle/v3 and seurat/v3 ###
#########################################################################

doObjectTransform <- function(x, transform_to = NULL){
    if (transform_to == "seurat"){
        count = exprs(x)
        meta.data = data.frame(pData(x))
        obj = CreateSeuratObject(count, meta.data = meta.data)
        return(obj)
    } else {
        count = GetAssayData(object = x, slot = "counts")
        meta.data = data.frame(x[[]])
        mouse_gene_sub = mouse_gene[rownames(count),]
        rownames(mouse_gene_sub) = rownames(count)
        
        cds <- new_cell_data_set(count,
                                 cell_metadata = meta.data,
                                 gene_metadata = mouse_gene_sub)
        return(cds)
    }
}


##############################################################
### Function: extract subset of cells to perform embedding ###
##############################################################

### We provided UMI count matrix for individual experiments, please download those data if you going to apply this function.
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax/download/mtx/

doExtractData <- function(df_cell, df_gene){
    experiment_list = paste0("run_", c(4, 13, 14, 15, 16, "17_sub1", "17_sub2", 18,
                                       19, 20, 21, 22, 23, 24, 25, 26, 27))
    gene_count = NULL
    for(i in experiment_list){
        print(paste0(i, "/", length(experiment_list)))
        
        gene_count_tmp = Matrix::readMM(paste0("gene_count.", i, "mtx.gz"))
        df_gene_tmp = read.csv("gene_annotation.csv.gz", row.names=1, as.is=T)
        df_cell_tmp = read.csv(paste0("cell_annotation.", i, ".csv.gz"), row.names=1, as.is=T)
        rownames(gene_count_tmp) = as.vector(df_gene_tmp$gene_ID)
        colnames(gene_count_tmp) = as.vector(df_gene_tmp$cell_id)
        
        gene_count = cbind(gene_count, 
                           gene_count_tmp[rownames(gene_count_tmp) %in% as.vector(df_gene$gene_ID), 
                                          colnames(gene_count_tmp) %in% as.vector(df_cell$cell_id), drop=FALSE])
    }
    
    gene_count = gene_count[,as.vector(df_cell$cell_id)]
    return(gene_count)
}




