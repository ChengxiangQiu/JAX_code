
#################################
### Section - 8, Birth series ###
#################################

###########################################
### Plotting cell number across samples ###
###########################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

df = readRDS(paste0(work_path, "pd_birth.rds"))
day = rep(NA, nrow(df))
day[df$RT_group == "E18.75_L1-01"] = "10"
day[df$RT_group == "E18.75_L2-01"] = "1"
day[df$RT_group == "E18.75_L2-02"] = "2"
day[df$RT_group == "E18.75_L2-03"] = "3"
day[df$RT_group == "E18.75_L2-04"] = "6"
day[df$RT_group == "E18.75_L2-05"] = "7"
day[df$RT_group == "E18.75_L2-06"] = "8"
day[df$RT_group == "E18.75_L2-07"] = "4"
day[df$RT_group == "E18.75_L2-08"] = "5"
day[df$RT_group == "E18.75_L2-09"] = "9"
day[df$RT_group == "P0_L1-01"    ] = "11"
day[df$RT_group == "P0_L1-04"    ] = "12"
df$tmp = as.vector(day)

pd_cell_num_1 = df %>% group_by(tmp, day) %>% tally() %>% rename(cell_num = n) %>% as.data.frame()
pd_cell_num_1$tmp = factor(pd_cell_num_1$tmp, levels = rev(1:12))

### Extended Data Fig. 12a

p = pd_cell_num_1 %>%
    ggplot(aes(tmp, cell_num, fill = day)) + 
    geom_bar(stat="identity") +
    coord_flip() +
    scale_fill_manual(values = birth_color_plate) + 
    geom_text(aes(label = scales::comma(cell_num)), 
              hjust = -0.1,
              position = position_dodge(width = 1),
              inherit.aes = TRUE,
              size = 5) +
    labs(x = "", y = "Cells profiled") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))


##################################################
### 2D UMAP of cells from birth series dataset ###
##################################################

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

pd = readRDS(paste0(work_path, "pd_birth.rds"))

### Extended Data Fig. 12b

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.5) +
    geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = major_trajectory), size=0.3) +
    scale_color_manual(values=major_trajectory_color_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, "birth_anno.png"), width = 10, height = 10, dpi = 300)


#########################################################################
### 2D UMAP of subclustering results for adipocytes and lung & airway ###
#########################################################################

##########################
###  1 - Hepatocytes

source("JAX_help_code.R")
source("JAX_color_code.R")
work_path = "./"

trajectory_i = "Hepatocytes"
df = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
df$day = paste0("Csection_", as.vector(df$day), "m")

### Fig. 6e (1st row)

birth_color_plate = c(birth_color_plate, "other" = "grey80")

for(i in paste0("Csection_", c(0,20,40,60,80), "m")){
    df$tmp = if_else(df$day == i, i, "other")
    try(df %>%
            ggplot() +
            geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=1.5, alpha = 0.3) +
            geom_point(data = subset(df, tmp == i),
                       aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=1.5) +
            theme_void() +
            scale_color_manual(values=birth_color_plate) +
            theme(legend.position="none") + 
            ggsave(paste0(work_path, trajectory_i, "_", i, ".png"), width = 4, height = 4, dpi = 300))
    df$tmp = NULL
}


##########################
###  2 - Adipocytes

trajectory_i = "Adipocytes"
df = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
df$day = paste0("Csection_", as.vector(df$day), "m")

### Extended Data Fig. 12c

color_plate = c("Brown adipocyte cells" = "#cb6751",
                "Adipocyte progenitor cells" = "#7aa457",
                "Adipocyte cells (Cyp2e1+)" = "#9e6ebd")

p = ggplot() +
    geom_point(data = df, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2, color="black") +
    geom_point(data = df, aes(x = UMAP_2d_1, y = UMAP_2d_2, color = anno_subclustering), size=1) +
    theme_void() +
    scale_color_manual(values=color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, trajectory_i, "_anno.png"), width = 4, height = 4, dpi = 300)

### Fig. 6e (2nd row)

for(i in paste0("Csection_", c(0,20,40,60,80), "m")){
    df$tmp = if_else(df$day == i, i, "other")
    try(df %>%
            ggplot() +
            geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=1.5, alpha = 0.3) +
            geom_point(data = subset(df, tmp == i),
                       aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=1.5) +
            theme_void() +
            scale_color_manual(values=birth_color_plate) +
            theme(legend.position="none") + 
            ggsave(paste0(work_path, trajectory_i, "_", i, ".png"), width = 4, height = 4, dpi = 300))
    df$tmp = NULL
}


##########################
###  3 - Lung & airway

trajectory_i = "Lung_and_airway"
df = readRDS(paste0(work_path, "Birth_series_", trajectory_i, "_Csections.obs.rds"))
df$day = paste0("Csection_", as.vector(df$day), "m")

### Extended Data Fig. 12d

color_plate = c("Airway club cells" = "#6dd9b4",
                "Alveolar Type 1 cells" = "#008cff",
                "Alveolar Type 2 cells" = "#dab300",
                "Lung cells (Eln+)" = "#185e3e",
                "Airway goblet cells" = "#663fc6")

p = ggplot() +
    geom_point(data = df, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1.2, color="black") +
    geom_point(data = df, aes(x = UMAP_2d_1, y = UMAP_2d_2, color = anno_subclustering), size=1) +
    theme_void() +
    scale_color_manual(values=color_plate) +
    theme(legend.position="none") + 
    ggsave(paste0(work_path, trajectory_i, "_anno.png"), width = 4, height = 4, dpi = 300)

### Fig. 6e (3rd row)

for(i in paste0("Csection_", c(0,20,40,60,80), "m")){
    df$tmp = if_else(df$day == i, i, "other")
    try(df %>%
            ggplot() +
            geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=1.5, alpha = 0.3) +
            geom_point(data = subset(df, tmp == i),
                       aes(x = UMAP_2d_1, y = UMAP_2d_2, color = tmp), size=1.5) +
            theme_void() +
            scale_color_manual(values=birth_color_plate) +
            theme(legend.position="none") + 
            ggsave(paste0(work_path, trajectory_i, "_", i, ".png"), width = 4, height = 4, dpi = 300))
    df$tmp = NULL
}

















