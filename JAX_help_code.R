
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
### posterior_embryo_gene_count.rds                       ### UMI count matrix of subset of cells for posterior embryo development
### posterior_embryo_adata_scale.NMP_Mesoderm.obs.csv     ### meta information of individual cells for posterior embryo subset (NMP)
### posterior_embryo_adata_scale.Notochord.obs.csv        ### meta information of individual cells for posterior embryo subset (Notochord)
### posterior_embryo_adata_scale.Gut.obs.csv              ### meta information of individual cells for posterior embryo subset (Gut)
### pd_somites.rds                                        ### meta information of individual cells for the somitogenesis validation dataset
### adata_somites_NMP.obs.csv                             ### meta information of individual cells for the NMP in the validation dataset

### Section_3_kidney_mesenchyme
### Renal_adata_scale.obs.rds                             ### meta information of renal subset
### Renal_CDI_adata_scale.obs.rds.                        ### meta information of a subset of renal cells
### 

################################################
### Packges that are needed for the analysis ###
################################################

library(monocle3)
library(Seurat)
require(lattice)
library(FNN)
library(dplyr)
library(ggplot2)
library(viridis)
library(Matrix)
library(stringr)
library(plotly)
library(htmlwidgets)
library(gridExtra) 

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
