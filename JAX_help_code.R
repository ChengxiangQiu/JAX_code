
###########################################
### Profiles that are used for analysis ###
###########################################

### 1. df_cell.csv                     ### meta information of individual cells
### 2. mouse.v12.geneID.txt            ### gene annotation based on GENCODE.M12
### 3. embryo_cds.rds                  ### pseudobulk dataset of individual embryos, in monocle/v3 format

### 4. posterior_embryo_gene_count.rds ### gene count of posterior embryos (somites 0-32)
### 5. pd_somites.rds                  ### meta information of individual cells for the somitogenesis validation dataset
### 6. adata_somites_NMP.obs.csv       ### meta information of individual cells for the NMP in the validation dataset

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
