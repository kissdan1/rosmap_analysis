library(tidyverse)
library(broom)
library(cowplot)
library(ggseg3d)
library(ggseg)
library(patchwork)
library(dplyr)
library(ggpubr)


#Data from rosmap 
rosmap_mgps = read_csv('/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmapMGPs.csv')

map_mri_thickness = read_csv('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_l_map_share.csv', na = 'NULL')
map_mri_thickness_2 = read_csv('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_r_map_share.csv', na = 'NULL')
map_mri_thickness = left_join(map_mri_thickness, map_mri_thickness_2)

ros_mri_thickness = read_csv('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_l_ros_share.csv', na = 'NULL')
ros_mri_thickness_2 = read_csv('/external/rprshnas01/external_data/rosmap/neuroimaging/data_transfer_june6_2022/freesurfer/v6/mri_cortical_thick_v6d_r_ros_share.csv', na = 'NULL')
ros_mri_thickness = left_join(ros_mri_thickness, ros_mri_thickness_2)

rosmap_mri_thickness = bind_rows(map_mri_thickness, ros_mri_thickness)
rosmap_mri_thickness = rosmap_mri_thickness %>% mutate(projid = as.numeric(projid))

joined_df = left_join(rosmap_mgps, rosmap_mri_thickness)

joined_df = joined_df %>% mutate(visit = as.numeric(visit), mri_from_death_years = age_death - (age_bl + visit))

#Reformat cortical areas to be compatible with ggseg 
cortical_areas = names(ros_mri_thickness)[-(38:40)]
cortical_areas = cortical_areas[-(1:2)]
cortical_areas = cortical_areas[-35]
cortical_areas_scaled = as.list(paste0("scale(", cortical_areas, ")"))
cortical_areas = substr(cortical_areas, start = 0, stop = (str_length(cortical_areas)-2))
cortical_areas = gsub(x = cortical_areas, pattern = "r_", replacement = 'rh_')
cortical_areas = gsub(x = cortical_areas, pattern = "l_", replacement = 'lh_')

#Which cortical areas correlate most with gpath? 
#Create a lm for each cortical area and compute statistics in list cortical_v_gpath_summaries
cortical_v_gpath = lapply((paste(cortical_areas_scaled, "~ scale(gpath) + as.factor(msex) +  scale(pmi) + as.factor(study) + scale(age_death)")), as.formula)
cortical_v_gpath_results = lapply(cortical_v_gpath, function(x) lm(x, data = joined_df))  
cortical_v_gpath_summaries = lapply(cortical_v_gpath_results, summary) 
cortical_v_gpath_summaries = lapply(cortical_v_gpath_summaries, tidy)
names(cortical_v_gpath_summaries) = cortical_areas

#Which cortical areas correlate most with cogdx? 
#Create a lm for each cortical area and compute statistics in list cortical_v_cogdx_summaries
cortical_v_cogdx = lapply((paste(cortical_areas_scaled, "~ scale(cogdx) + as.factor(msex) +  scale(pmi) + as.factor(study) + scale(age_death)")), as.formula)
cortical_v_cogdx_results = lapply(cortical_v_cogdx, function(x) lm(x, data = joined_df))  
cortical_v_cogdx_summaries = lapply(cortical_v_cogdx_results, summary) 
cortical_v_cogdx_summaries = lapply(cortical_v_cogdx_summaries, tidy)
names(cortical_v_cogdx_summaries) = cortical_areas

#Do certain cortical areas correlate with certain cell types? 
#Create a lm for each cortical area and compute statistics in list cortical_v_cells_summaries
cells = joined_df[,3:21] %>% colnames()
cortical_v_cells = list()
cortical_v_cells_results = list()
cortical_v_cells_summaries = list()
for(cell in cells) {
  cortical_v_cells[[cell]] = lapply((paste(cortical_areas_scaled, (paste0("~ scale(", cell, ") + scale(pmi)")))), as.formula)
  names(cortical_v_cells[[cell]]) = cortical_areas
  cortical_v_cells_results[[cell]] = lapply(cortical_v_cells[[cell]], function(x) lm(x, data = joined_df))
  cortical_v_cells_summaries[[cell]] = lapply(cortical_v_cells_results[[cell]], summary)
  cortical_v_cells_summaries[[cell]] = lapply(cortical_v_cells_summaries[[cell]], tidy)
  names(cortical_v_cells_summaries[[cell]]) = cortical_areas
}

#Visualizing summary statistics 

#Plot beta values for gpath 
gpath_df = data.frame()
for( x in 1:length(cortical_v_gpath_summaries)) {
  gpath_df = rbind(gpath_df, cortical_v_gpath_summaries[[x]][2,])
}
gpath_df$cortical_area = names(cortical_v_gpath_summaries)
gpath_df = subset(gpath_df) #Can set pval threshold here
gpath_df = gpath_df[order(gpath_df$estimate),]
gpath_df$hemisphere = substr(gpath_df$cortical_area, start = 1, stop = 1)
#Plot gpath B estimate for each cortical area where p is significant
gpath_plot = gpath_df %>% ggplot(aes(x = reorder(cortical_area, estimate), y = estimate, fill = hemisphere)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Effect of gpath on Cortical Thickness") + ylab("gpath β Estimate") + xlab("Cortical Area")
plot(gpath_plot)

#Plot beta values for cogdx 
cogdx_df = data.frame()
for( x in 1:length(cortical_v_cogdx_summaries)) {
  cogdx_df = rbind(cogdx_df, cortical_v_cogdx_summaries[[x]][2,])
}
cogdx_df$cortical_area = names(cortical_v_cogdx_summaries)
cogdx_df = subset(cogdx_df) #Can set pval threshold here
cogdx_df = cogdx_df[order(cogdx_df$estimate),]
cogdx_df$hemisphere = substr(cogdx_df$cortical_area, start = 1, stop = 1)
#Plot cogdx B estimate for each cortical area where p is significant
cogdx_plot = cogdx_df %>% ggplot(aes(x = reorder(cortical_area, estimate), y = estimate, fill = hemisphere)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Effect of cogdx on Cortical Thickness") + ylab("cogdx β Estimate") + xlab("Cortical Area")
plot(cogdx_plot)

#Graph gpath using ggseg to project cortical areas onto brain maps 
ggseg_gpath = data.frame(matrix(nrow = 68, ncol = 0))
ggseg_gpath$label = gpath_df$cortical_area
ggseg_gpath$cortical_area = gpath_df$estimate

ggseg_gpath %>% ggseg(mapping=aes(fill=(cortical_area)), position = "stacked", colour = "black") + ggtitle("Effect of gpath on Cortical Thickness ")

#Graph cogdx using ggseg to project cortical areas onto brain maps 
ggseg_cogdx = data.frame(matrix(nrow = 68, ncol = 0))
ggseg_cogdx$label = cogdx_df$cortical_area
ggseg_cogdx$cortical_area = cogdx_df$estimate

ggseg_cogdx %>% ggseg(mapping=aes(fill=(cortical_area)), position = "stacked", colour = "black", atlas = "dk") + ggtitle("Effect of cogdx on Cortical Thickness ") + scale_fill_gradient2(low="red",high="green", mid = "white")

#Plot beta values for ALL cell types (multiple graphs)
#Make a list of data for each cell type vs. each cortical area 
cell_list = list()
ggseg_cell_list = list()
for(y in 1:length(cortical_v_cells_summaries)) {
  cell_df = data.frame(matrix(, nrow=68, ncol=0))
  ggseg_cell = data.frame(matrix(, nrow=68, ncol=0))
  estimates = c()
  pvals = c()
  for(x in 1:length(cortical_v_cells_summaries[[y]])) {
    cell_df$cortical_area = names(cortical_v_cells_summaries[[y]])
    estimates = append(estimates, cortical_v_cells_summaries[[y]][x][[1]]$estimate[2])
    pvals = append(pvals, cortical_v_cells_summaries[[y]][x][[1]]$p.value[2])
  }
  cell_df$estimate = estimates
  cell_df$p.value = pvals
  cell_df$hemisphere = substr(cell_df$cortical_area, start = 1, stop = 1)
  ggseg_cell$label = cell_df$cortical_area
  ggseg_cell$cortical_area = cell_df$estimate
  cell_list[[y]] = cell_df
  ggseg_cell_list[[y]] = ggseg_cell
}
names(cell_list) = names(cortical_v_cells_summaries)
names(ggseg_cell_list) = names(cortical_v_cells_summaries)

#Graph each cell type as ggplot bar plot
for(x in 1:length(cell_list)) {
  cell_plot = cell_list[[x]] %>% ggplot(aes(x = reorder(cortical_area, estimate), y = estimate, fill = hemisphere)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + ggtitle(paste0("Effect of variation in ", names(cell_list)[x], " Cell Proportion on Cortical Thickness")) + ylab("Cell Type Proportion β Estimate") + xlab("Cortical Area")
  plot(cell_plot)
}
#Graph each cell type as ggpseg brain map 
for(x in 1:length(cell_list)) {
  ggseg_cell_plot = ggseg_cell_list[[x]] %>% ggseg(mapping=aes(fill=(cortical_area)), position = "stacked", colour = "black") + ggtitle(paste0("Effect of ", names(ggseg_cell_list)[x], " Cell Proportion on Cortical Thickness")) + scale_fill_gradient2(low="red",high="green", mid = "white")
  plot(ggseg_cell_plot)
}
