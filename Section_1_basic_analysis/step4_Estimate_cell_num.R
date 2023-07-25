
###################################
### Section - 1, Basic analysis ###
###################################

######################################################################
### Estimating absolute number of cells from individual timepoints ###
######################################################################

### Cell number estimated by qPCR experiment (million)
cell_num = c("E8.5" = 0.21,
             "E9.5" = 0.94,
             "E10.5" = 10.10,
             "E11.5" = 22.98,
             "E12.5" = 45.03,
             "E13.5" = 60.59,
             "E14.5" = 131.00,
             "E15.5" = 216.79,
             "E16.5" = 353.17,
             "E17.5" = 515.85,
             "E18.5" = 584.78,
             "E19.5" = 671.50)

df = data.frame(x = as.numeric(gsub("E", "", as.vector(names(cell_num)))),
                log2_y = log2(cell_num * 1000000))

### fit polynomial regression with degree 5
fit3 = lm(log2_y~poly(x,3,raw=TRUE), data=df)
print(summary(fit3)$adj.r.squared) ###  0.9858479

day_list = c("E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75", "E10.0", "E10.25", 
             "E10.5", "E10.75", "E11.0", "E11.25", "E11.5", "E11.75", "E12.0", 
             "E12.25", "E12.5", "E12.75", "E13.0", "E13.25", "E13.5", "E13.75", 
             "E14.0", "E14.25", "E14.333", "E14.75", "E15.0", "E15.25", "E15.5", 
             "E15.75", "E16.0", "E16.25", "E16.5", "E16.75", "E17.0", "E17.25", 
             "E17.5", "E17.75", "E18.0", "E18.25", "E18.5", "E18.75", "E19.5")

x_axis = as.numeric(gsub("E","",day_list))

plot(df$x, df$log2_y, pch=19, xlab='x', ylab='log2_y')
lines(x_axis, predict(fit3, data.frame(x=x_axis)), col='purple')

cell_num_pred = round(2^predict(fit3, data.frame(x=x_axis)))

df_x = data.frame(day = x_axis,
                  cell_num_pred_log2 = predict(fit3, data.frame(x=x_axis)),
                  cell_num_pred = round(2^predict(fit3, data.frame(x=x_axis))))

day_x = day = paste0("E", df_x$day)
day[day_x == "E9"] = "E9.0"
day[day_x == "E10"] = "E10.0"
day[day_x == "E11"] = "E11.0"
day[day_x == "E12"] = "E12.0"
day[day_x == "E13"] = "E13.0"
day[day_x == "E14"] = "E14.0"
day[day_x == "E15"] = "E15.0"
day[day_x == "E16"] = "E16.0"
day[day_x == "E17"] = "E17.0"
day[day_x == "E18"] = "E18.0"
day[day_x == "E19.5"] = "P0"
df_x$day = as.vector(day)

### summary of fit3 curve
a3 = 0.011369
a2 = -0.583861
a1 = 10.397036
a4 = -35.469755

### the function of curve
p_function = function(x){
    y = a3*x^3 + a2*x^2 + a1*x^1 + a4
    return(y)
}

### derivative of the curve, which is increasing time (in log2 scale) at a given timepoint
d_function = function(x){
    y = 3*a3*x^2 + 2*a2*x + a1
    return(y)
}

### timepoint
x_axis = as.numeric(gsub("E","",day_list))

### the doubling time
doubling_time = 24*2/(2^d_function(x_axis))

df_x$doubling_time = doubling_time
df_x$x_axis = x_axis

write.csv(df_x, "cell_num_prediction.csv")


################################################################################################
### Plotting the cell composition of each major trajectories as a function of time (Fig. 1e) ###
################################################################################################

source("JAX_help_code.R")
source("JAX_color_code.R")

pd = readRDS("df_cell.rds")
### n = 11,441,407 cells

x = as.vector(pd$day)
x[pd$day == "E8.0-E8.5"] = "E8.5"
pd$day = as.vector(x)

cell_num = read.csv("cell_num_prediction.csv")

df = pd %>%
    group_by(day, major_trajectory) %>%
    tally() %>%
    dplyr::rename(cell_num = n)
df_sub = pd %>%
    group_by(day) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)
df = df %>%
    left_join(df_sub, by = "day") %>%
    mutate(percentage = cell_num/cell_num_total) %>%
    left_join(cell_num, by = "day") %>%
    mutate(cell_num_pred_log2 = cell_num_pred_log2 * percentage)
df$day = factor(df$day, levels = names(day_color_plate))
df$major_trajectory = factor(df$major_trajectory, levels = names(major_trajectory_color_plate))

cell_num_x = c("E8.5" = 0.21,
               "E9.5" = 0.94,
               "E10.5" = 10.10,
               "E11.5" = 22.98,
               "E12.5" = 45.03,
               "E13.5" = 60.59,
               "E14.333" = 131.00,
               "E15.5" = 216.79,
               "E16.5" = 353.17,
               "E17.5" = 515.85,
               "E18.5" = 584.78,
               "P0" = 671.50)

df_y = data.frame(day = as.vector(names(cell_num_x)),
                  log2_y = log2(cell_num_x * 1000000))

df_y$day = factor(df_y$day, levels = names(day_color_plate))

cell_num$day = factor(cell_num$day, levels = names(day_color_plate))

p = ggplot() +
    geom_bar(data = df, aes(x = day, y = cell_num_pred_log2, group = major_trajectory, fill = major_trajectory), stat="identity", width = 1) +
    geom_point(data=df_y, aes(x=day, y=log2_y), shape = 21, colour = "black", fill = "white", size = 2, stroke = 1.5, alpha = 0.8) +
    labs(x = "", y = "") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = major_trajectory_color_plate) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1), axis.text.y = element_text(color="black"))

pdf("Cell_composition_over_time.pdf", 7, 5)
print(p)
dev.off()







