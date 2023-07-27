
#####################################
### Section - 2, Posterior embryo ###
#####################################

#######################################################################################
### Which genes are correlated with each PC identified from NMPs, Notochord, or Gut ###
#######################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

celltype_list = c("NMP_Mesoderm", "Notochord", "Gut")
example_i = "posterior_embryo"

pd = read.csv(paste0(work_path, example_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count = doExtractData(pd, mouse_gene_sub)

for(i in celltype_list){
    print(i)
    
    pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)
    res = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA.rds"))
    emb = res[["cell.embeddings"]]
    pd_x = pd_x[rownames(emb),]
    
    gene_count_x = gene_count[,rownames(emb)]
    obj_x = CreateSeuratObject(gene_count_x, meta.data = pd_x)
        
    obj_x = NormalizeData(obj_x, normalization.method = "LogNormalize", scale.factor = 10000)
    obj_x = FindVariableFeatures(obj_x, selection.method = "vst", nfeatures = 2500)
    obj_x = ScaleData(obj_x, verbose = FALSE)
    scale_dat = GetAssayData(obj_x, slot = "scale.data")
    
    if(nrow(emb) > 20000){
        keep = sample(1:nrow(emb), 20000)
    } else {
        keep = 1:nrow(emb)
    }
    
    emb_sub = emb[keep,]
    scale_dat_sub = scale_dat[,keep]
    
    corr_matrix = matrix(NA, 2500, 3)
    pval_matrix = matrix(NA, 2500, 3)
    for(j in 1:3){
        for(k in 1:2500){
            print(k)
            fit = cor.test(emb_sub[,j], scale_dat_sub[k,])
            corr_matrix[k,j] = fit$estimate
            pval_matrix[k,j] = fit$p.value
        }
    }
    
    mouse_gene_sub = mouse_gene[rownames(scale_dat),]
    
    df = data.frame(corr = c(corr_matrix[,1], corr_matrix[,2], corr_matrix[,3]),
                    pval = c(pval_matrix[,1], pval_matrix[,2], pval_matrix[,3]),
                    PC = rep(c("PC_1", "PC_2", "PC_3"), each = 2500),
                    gene_short_name = rep(as.vector(mouse_gene_sub$gene_short_name), 3), stringsAsFactors = F)
    
    saveRDS(df, paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.rds"))
    
    ### further filter significant ones
    names(df) = c("corr", "pval", "PC", "gene_short_name")
    
    df_sub = df %>% filter(!is.na(pval), PC == "PC_1") %>%
        mutate(qval = p.adjust(pval, method = "fdr")) 
    x1 = mean(df_sub$corr) - sd(df_sub$corr)
    x2 = mean(df_sub$corr) + sd(df_sub$corr)
    df_1 = df_sub %>% filter(qval < 0.05) %>%
        filter(corr < x1 | corr > x2) %>%
        arrange(corr) %>%
        select(-pval)
    
    df_sub = df %>% filter(!is.na(pval), PC == "PC_2") %>%
        mutate(qval = p.adjust(pval, method = "fdr")) 
    x1 = mean(df_sub$corr) - sd(df_sub$corr)
    x2 = mean(df_sub$corr) + sd(df_sub$corr)
    df_2 = df_sub %>% filter(qval < 0.05) %>%
        filter(corr < x1 | corr > x2) %>%
        arrange(corr) %>%
        select(-pval)
    
    df_sub = df %>% filter(!is.na(pval), PC == "PC_3") %>%
        mutate(qval = p.adjust(pval, method = "fdr")) 
    x1 = mean(df_sub$corr) - sd(df_sub$corr)
    x2 = mean(df_sub$corr) + sd(df_sub$corr)
    df_3 = df_sub %>% filter(qval < 0.05) %>%
        filter(corr < x1 | corr > x2) %>%
        arrange(corr) %>%
        select(-pval)
    
    df_out = rbind(df_1, df_2, df_3)
    
    write.csv(df_out, paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.csv"))
    
}



############################################################################
### Identify differential expressed genes between early vs. late somites ###
############################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "posterior_embryo"

i = "NMP_Mesoderm"

pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count_x = doExtractData(pd_x, mouse_gene_sub)

mouse_gene %>% filter(gene_short_name %in% c("T","Meis1"))

pd_x$Meis1_exp = as.vector(gene_count_x["ENSMUSG00000020160",])
pd_x$T_exp = as.vector(gene_count_x["ENSMUSG00000062327",])
pd_x$NMP = if_else(pd_x$T_exp >= 5 & pd_x$Meis1_exp == 0, "yes", "no")
table(pd_x$NMP)
table(pd_x$group[pd_x$NMP=="yes"])

NMP_early_late = rep("no", nrow(pd_x))
NMP_early_late[pd_x$NMP=="yes" & pd_x$group == "E85"] = "early"
NMP_early_late[pd_x$NMP=="yes" & pd_x$group != "E85"] = "late"
pd_x$NMP_early_late = NMP_early_late
p = ggplot() +
    geom_point(data = pd_x, aes(x = UMAP_2d_1, y = UMAP_2d_2, color = NMP_early_late), size=0.35) +
    geom_point(data = pd_x[pd_x$NMP == "yes",], aes(x = UMAP_2d_1, y = UMAP_2d_2, color = NMP_early_late), size=0.5) +
    theme_void() +
    scale_color_manual(values=c("early"="blue", "late"="red", "no"="grey80")) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "Early_late_NMPs.png"), width = 4, height = 3, dpi = 300)


### n = 8859 cells for NMP
pd_sub = pd_x[pd_x$NMP == "yes",]
gene_count_sub = gene_count[,rownames(pd_sub)]
obj = CreateSeuratObject(gene_count_sub, meta.data = pd_sub)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)
Idents(obj) = if_else(obj$group == "E85", "early", "late")

res = FindMarkers(obj, ident.1 = "early", ident.2 = "late", features = genes_include, logfc.threshold = 0)

res %>% mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>%
    select(-p_val) %>%
    mutate(highly_expressed = if_else(avg_logFC > 0, "early", "late")) %>%
    write.csv(paste0(work_path, example_i, "_adata_scale.", i, ".DEGs_early_late.csv"))

### filtering significant ones
res %>% mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene %>% select(gene_ID, gene_short_name), by = "gene_ID") %>%
    filter(p_val_adj < 0.05, abs(avg_logFC) > 0.25) %>%
    select(-p_val) %>%
    mutate(highly_expressed = if_else(avg_logFC > 0, "early", "late")) %>%
    write.csv(paste0(work_path, example_i, "_adata_scale.", i, ".DEGs_early_late_sig.csv"))


##########################################################
### Identifying overlaps between different comparisons ###
##########################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

example_i = "posterior_embryo"

i = "NMP_Mesoderm"
df_all_1 = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.rds"))
names(df_all_1) = c("NMP_corr", "NMP_pval", "PC", "gene_short_name")
df_all_1 = df_all_1[!is.na(df_all_1$NMP_pval),]
df_sig_1 = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.csv"))

i = "Notochord"
df_all_2 = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.rds"))
names(df_all_2) = c("Notochord_corr", "Notochord_pval", "PC", "gene_short_name")
df_all_2 = df_all_2[!is.na(df_all_2$Notochord_pval),]
df_sig_2 = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.csv"))

i = "Gut"
df_all_3 = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.rds"))
names(df_all_3) = c("Gut_corr", "Gut_pval", "PC", "gene_short_name")
df_all_3 = df_all_3[!is.na(df_all_3$Gut_pval),]
df_sig_3 = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".PCA_corr.csv"))

df_all_4 = read.csv(paste0(work_path, "posterior_embryo_adata_scale.NMP_Mesoderm.DEGs_early_late.csv"), row.names=1)
df_sig_4 = read.csv(paste0(work_path, "posterior_embryo_adata_scale.NMP_Mesoderm.DEGs_early_late_sig.csv"), row.names=1)
NMP_early = df_sig_4[df_sig_4$highly_expressed == "early",]
NMP_late = df_sig_4[df_sig_4$highly_expressed == "late",]


###########################################################
### comparing PC1 of Notochord and PC1 of Gut (Fig. 2l) ###
###########################################################

gene_overlap = intersect(df_sig_2 %>% filter(PC == "PC_1") %>% pull(gene_short_name),
                         df_sig_3 %>% filter(PC == "PC_1") %>% pull(gene_short_name))

df = df_all_2 %>% filter(PC == "PC_1") %>% mutate(Notochord_fdr = p.adjust(Notochord_pval, method = "fdr")) %>%
    left_join(df_all_3 %>% filter(PC == "PC_1") %>% mutate(Gut_fdr = p.adjust(Gut_pval, method = "fdr")) %>% select(Gut_corr, Gut_fdr, gene_short_name), by = "gene_short_name") %>%
    filter(!is.na(Gut_corr)) %>%
    mutate(if_sig = if_else(gene_short_name %in% gene_overlap, "yes", "no")) 

df_sub = df[df$if_sig == "yes",]
df_sub$Notochord_PC = df_sub$Gut_PC = df_sub$PC
df_sub = df_sub[,c("gene_short_name", "Notochord_PC", "Notochord_corr", "Notochord_fdr",
                   "Gut_PC", "Gut_corr", "Gut_fdr")]
write.csv(df_sub, paste0(work_path, "Corr_Notochord_PC1_Gut_PC1.csv"))

### four quadrant
x1 = sum(Notochord_PC1_pos %in% Gut_PC1_pos)
x2 = sum(Notochord_PC1_pos %in% Gut_PC1_neg)
x3 = sum(Notochord_PC1_neg %in% Gut_PC1_pos)
x4 = sum(Notochord_PC1_neg %in% Gut_PC1_neg)
chisq.test(c(x1,x2,x3,x4))$p.value
### p-value = 3.013274e-29

p1 = df %>%
    ggplot(aes(Notochord_corr, Gut_corr, color = if_sig)) + geom_point() + 
    labs(x="Corr_Notochord_PC1", y="Corr_Gut_PC1", title="") +
    theme_minimal(base_size = 15) +
    theme(legend.position="none") +
    scale_color_manual(values=c("yes" = "blue", "no" = "grey80")) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-0.9,0.9) + ylim(-0.9,0.9) + geom_hline(yintercept=0) + geom_vline(xintercept=0)  +
    geom_text(data=subset(df, (Notochord_corr >= 0.3 & Gut_corr >= 0.3) | (Notochord_corr <= (-0.5) & Gut_corr <= (-0.5)) & if_sig == "yes"),
              aes(Notochord_corr, Gut_corr, label=gene_short_name), hjust = 0, nudge_x = 0.01, color = "red", size = 2)

pdf(paste0(work_path, "Corr_Notochord_PC1_Gut_PC1.pdf"), 5, 5)
print(p1)
dev.off()

### line plot 

example_i = "posterior_embryo"
pd = read.csv(paste0(work_path, example_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count = doExtractData(pd, mouse_gene_sub)

wnt_genes = c("Ptk7","Rspo3","Axin2","Fzd10","Lef1","Nkd1","Rnf43","Tle4","Wnt3a")

i = "Notochord"
pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)
res = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA.rds"))
emb = res[["cell.embeddings"]]
pd_x = pd_x[rownames(emb),]
gene_count_x = gene_count[,rownames(pd_x)]

mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% wnt_genes,]

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(t(as.matrix(gene_count_x))),
                gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                PC_1 = rep(emb[,1], nrow(gene_count_x)), stringsAsFactors = F)

pdf(paste0(work_path, "Notochord_marker_genes.pdf"), 5, 3)
df %>%
    ggplot(aes(PC_1, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 
dev.off()


i = "Gut"
pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)
res = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA.rds"))
emb = res[["cell.embeddings"]]
pd_x = pd_x[rownames(emb),]
gene_count_x = gene_count[,rownames(pd_x)]

mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% wnt_genes,]

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(t(as.matrix(gene_count_x))),
                gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                PC_1 = rep(emb[,1], nrow(gene_count_x)), stringsAsFactors = F)

pdf(paste0(work_path, "Gut_marker_genes.pdf"), 5, 3)
df %>%
    ggplot(aes(PC_1, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 
dev.off()



##############################################################
### comparing Early vs. Late NMPs and PC2 of Gut (Fig. 2m) ###
##############################################################

gene_overlap = intersect(df_sig_4 %>% pull(gene_short_name),
                         df_sig_3 %>% filter(PC == "PC_2") %>% pull(gene_short_name))

df = df_all_4 %>% select(avg_logFC, pct.1, pct.2, p_val_adj, highly_expressed, gene_short_name) %>% 
    left_join(df_all_3 %>% filter(PC == "PC_2") %>% mutate(Gut_fdr = p.adjust(Gut_pval, method = "fdr")) %>% select(PC, Gut_corr, Gut_fdr, gene_short_name), by = "gene_short_name") %>%
    filter(!is.na(Gut_corr)) %>%
    mutate(if_sig = if_else(gene_short_name %in% gene_overlap, "yes", "no")) 

df_sub = df[df$if_sig == "yes",]
df_sub = df_sub[,c("gene_short_name", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "highly_expressed", 
                   "PC", "Gut_corr", "Gut_fdr")]
write.csv(df_sub, paste0(work_path, "NMP_early_late_Gut_PC2_overlap.csv"))

NMP_early = df_sub %>% filter(avg_logFC > 0) %>% pull(gene_short_name)
NMP_late = df_sub %>% filter(avg_logFC < 0) %>% pull(gene_short_name)
Gut_PC2_pos = df_sub %>% filter(Gut_corr > 0) %>% pull(gene_short_name)
Gut_PC2_neg = df_sub %>% filter(Gut_corr < 0) %>% pull(gene_short_name)

### four quadrant
x1 = sum(NMP_late %in% Gut_PC2_pos)
x2 = sum(NMP_late %in% Gut_PC2_neg)
x3 = sum(NMP_early %in% Gut_PC2_pos)
x4 = sum(NMP_early %in% Gut_PC2_neg)
chisq.test(c(x1,x2,x3,x4))$p.value
### p-value = 1.617402e-16

overlap_gene_early = intersect(NMP_early, Gut_PC2_pos)

p2 = df %>%
    ggplot(aes(avg_logFC, Gut_corr, color = if_sig)) + geom_point() + 
    labs(x="avg_logFC_NMPs", y="Corr_Gut_PC2", title="") +
    theme_minimal(base_size = 15) +
    theme(legend.position="none") +
    scale_color_manual(values=c("yes" = "blue", "no" = "grey80")) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-1.5,1.5) + ylim(-0.7,0.7) + geom_hline(yintercept=0) + geom_vline(xintercept=0)  +
    geom_text(data=subset(df, (avg_logFC >= 0.5 & Gut_corr > 0.2) | (avg_logFC <= (-0.5) & Gut_corr < (-0.2)) & if_sig == "yes"),
              aes(avg_logFC, Gut_corr, label=gene_short_name), hjust = 0, nudge_x = 0.01, color = "red", size = 2)

pdf(paste0(work_path, "NMP_early_late_Gut_PC2_overlap.pdf"), 5, 5)
print(p2)
dev.off()

### line plot 
example_i = "posterior_embryo"
pd = read.csv(paste0(work_path, example_i, "_adata_scale.obs.csv"), header=T, row.names=1, as.is=T)
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]
gene_count = doExtractData(pd, mouse_gene_sub)

target_genes = c("Hspd1","Npm1","Eno1","Ncl","Lin28a","Hsp90aa1")

i = "NMP_Mesoderm"
pd_x = read.csv(paste0(example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)
res = readRDS(paste0(example_i, "_adata_scale.", i, ".PCA.rds"))
emb = res[["cell.embeddings"]]
pd_x = pd_x[rownames(emb),]
gene_count_x = gene_count[,rownames(pd_x)]

mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

gene_count_early = gene_count_x[,colnames(gene_count_x) %in% colnames(obj)[obj$group=="E85"]]
gene_count_late = gene_count_x[,colnames(gene_count_x) %in% colnames(obj)[obj$group!="E85"]]

df = rbind(data.frame(exp = c(t(as.matrix(gene_count_early))),
                gene = rep(rownames(gene_count_early), each = ncol(gene_count_early)), 
                group = "early",stringsAsFactors = F),
           data.frame(exp = c(t(as.matrix(gene_count_late))),
                      gene = rep(rownames(gene_count_late), each = ncol(gene_count_late)),
                      group = "late", stringsAsFactors = F))

pdf(paste0(work_path, "NMPs_marker_genes.pdf"), 4, 3)
df %>%
    ggplot( aes(gene, exp, fill = group)) + 
    geom_boxplot(outlier.shape = NA) + 
    labs(x="", y="exp") +
    theme_classic(base_size = 5) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
dev.off()



i = "Gut"
pd_x = read.csv(paste0(work_path, example_i, "_adata_scale.", i, ".obs.csv"), header=T, row.names=1, as.is=T)
res = readRDS(paste0(work_path, example_i, "_adata_scale.", i, ".PCA.rds"))
emb = res[["cell.embeddings"]]
pd_x = pd_x[rownames(emb),]
gene_count_x = gene_count[,rownames(pd_x)]

mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% target_genes,]

gene_count_x = t(t(gene_count_x) / colSums(gene_count_x)) * 100000
gene_count_x = gene_count_x[rownames(mouse_gene_sub),]
rownames(gene_count_x) = as.vector(mouse_gene_sub$gene_short_name)
gene_count_x@x = log(gene_count_x@x + 1)

df = data.frame(exp = c(t(as.matrix(gene_count_x))),
                gene = rep(rownames(gene_count_x), each = ncol(gene_count_x)),
                PC_2 = rep(emb[,2], nrow(gene_count_x)), stringsAsFactors = F)

pdf(paste0(work_path, "Gut_marker_genes.pdf"), 5, 3)
df %>%
    ggplot(aes(PC_2, exp, color = gene)) + geom_smooth(method = loess, se = FALSE) +
    labs(x="", y="", title="") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_color_brewer(palette = "Set1") 
dev.off()

