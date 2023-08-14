
####################################
### Section - 7, Key TFs & genes ###
####################################

### making plot to show which TF or gene are appeared for different edges

################
### key TFs ####
################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

dat = read.csv(paste0(work_path, "All.keyTF.csv"), header=T, as.is=T)

df = dat %>% select(node_A, node_B, gene_short_name) %>% unique() %>%
    group_by(node_A, node_B) %>% tally()
print(paste0(mean(df$n), "+/-", sd(df$n)))
### 39.64 +/- 43.44

print(paste(quantile(df$n, 0.25), quantile(df$n, 0.5), quantile(df$n, 0.75)))
### 12 28 51

df_1 = dat %>% filter(comparing %in% c("group_1 vs. group_2","group_3 vs. group_4")) %>% group_by(gene_short_name) %>% tally() %>% pull(gene_short_name)
df_2 = dat %>% filter(comparing %in% c("group_2 vs. group_3")) %>% group_by(gene_short_name) %>% tally() %>% pull(gene_short_name)
df_3 = dat %>% group_by(gene_short_name) %>% tally() %>% pull(gene_short_name)
sum(!df_1 %in% df_2)/length(df_3)
### 5%

df_1 = dat %>% filter(edge_type == "Developmental progression", comparing == "group_1 vs. group_2") %>% 
    select(node_A, node_B, gene_short_name) %>% unique() %>%
    group_by(gene_short_name) %>% tally() 
print(head(df_1[order(df_1$n, decreasing = T),], 20))
p1 = df_1 %>% 
    ggplot(aes(n)) + geom_histogram(binwidth = 0.5) +
    labs(x="# of edges that have been involved", y="# of key TFs", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

df_3 = dat %>% 
    select(node_A, node_B, gene_short_name) %>% unique() %>%
    group_by(gene_short_name) %>% tally()
print(head(df_3[order(df_3$n, decreasing = T),], 20))
p3 = df_3 %>%
    ggplot(aes(n)) + geom_histogram(binwidth = 0.6) +
    labs(x="# of edges that have been involved", y="# of key TFs", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

### Fig. S22c

pdf(paste0(work_path, "Hist_TF.pdf"), 6, 3)
grid.arrange(p1, p3, nrow=1, ncol=2) 
dev.off()

##################
### key Genes ####
##################


dat = read.csv(paste0(work_path, "All.keyGene.csv"), header=T, as.is=T)

df = dat %>% select(node_A, node_B, gene_short_name) %>% unique() %>%
    group_by(node_A, node_B) %>% tally()
print(paste0(mean(df$n), "+/-", sd(df$n)))
### 293.24 +/- 358.04

print(paste(quantile(df$n, 0.25), quantile(df$n, 0.5), quantile(df$n, 0.75)))
### 76 171 389


df_1 = dat %>% filter(comparing %in% c("group_1 vs. group_2","group_3 vs. group_4")) %>% group_by(gene_short_name) %>% tally() %>% pull(gene_short_name)
df_2 = dat %>% filter(comparing %in% c("group_2 vs. group_3")) %>% group_by(gene_short_name) %>% tally() %>% pull(gene_short_name)
df_3 = dat %>% group_by(gene_short_name) %>% tally() %>% pull(gene_short_name)
sum(!df_1 %in% df_2)/length(df_3)
### 7%


df_1 = dat %>% filter(edge_type == "Developmental progression", comparing == "group_1 vs. group_2") %>% 
    select(node_A, node_B, gene_short_name) %>% unique() %>%
    group_by(gene_short_name) %>% tally() 
print(head(df_1[order(df_1$n, decreasing = T),], 20))
df_1[order(df_1$n, decreasing = T),] %>% filter(n > 10) %>% write.csv("~/Dropbox/tmp/Fig.S16.d_1.csv")
p1 = df_1 %>% 
    ggplot(aes(n)) + geom_histogram(binwidth = 0.5) +
    labs(x="# of edges that have been involved", y="# of key Genes", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

df_3 = dat %>% 
    select(node_A, node_B, gene_short_name) %>% unique() %>%
    group_by(gene_short_name) %>% tally()
print(head(df_3[order(df_3$n, decreasing = T),], 20))
df_3[order(df_3$n, decreasing = T),] %>% filter(n > 10) %>% write.csv("~/Dropbox/tmp/Fig.S16.d_3.csv")
p3 = df_3 %>%
    ggplot(aes(n)) + geom_histogram(binwidth = 0.6) +
    labs(x="# of edges that have been involved", y="# of key Genes", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

### Fig. S22d

pdf(paste0(work_path, "Hist_Gene.pdf"), 6, 3)
grid.arrange(p1, p3, nrow=1, ncol=2) 
dev.off()


