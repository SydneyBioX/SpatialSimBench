library(dplyr)
library(stringr)
library(tidyr)
library(funkyheatmap)
library(kableExtra)
library(ggplot2)
library(forestplot)
library(SingleCellExperiment)
library(MASS)
library(readr)
library(dotwhisker)
library(broom)

final_sample_by_ds <- readRDS("final_sample_by_ds.rds")
final_feature_by_ds <- readRDS("final_feature_by_ds.rds")
final_cluster_by_ds <- readRDS("final_cluster_by_ds.rds")
final_spatial_by_ds <- readRDS("final_spatial_by_ds.rds")

# panel a

processing_boxplot <- function(dataset, label){
  
  rownames(dataset) <- c("scDesign3_gau(rf)", "scDesign3_nb(rf)","scDesign3_poi(rf)","SRTsim(rf)","scDesign2" ,"scDesign3_gau","scDesign3_nb","scDesign3_poi", "SPARsim","splatter", "SRTsim","symsim","zinbwave")
  
  colnames(dataset) <- c("Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5", "Dataset6", "Dataset7","Dataset8", "Dataset9", "Dataset10")
  
  df <- as_tibble(dataset, .name_repair = "universal")
  
  
  df <- df %>% 
    mutate(model = rownames(dataset)) %>%
    pivot_longer(-model, names_to = "dataset", values_to = "value")
  
  df$dataset <- factor(df$dataset)
  
  df$dataset <- factor(df$dataset, levels = c("Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5", "Dataset6", "Dataset7","Dataset8", "Dataset9", "Dataset10"))
  
  g <- ggplot(data = df, aes(x = model, y = value)) +
    geom_boxplot(outlier.shape = NA) +  # This will draw boxplots for all data without distinguishing by dataset
    geom_jitter(aes(color = dataset), width = 0.2, alpha = 0.5) +  # Color points by dataset
    theme_minimal() +
    labs(x = "Model", y = label, color = "Dataset") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),  # Remove panel background
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.border = element_rect(colour = "black", fill=NA, size=1)  # Add border
    )
  
  return(g)
}


bp_sample <- processing_boxplot(final_sample_by_ds, "zstat")
bp_feature <- processing_boxplot(final_feature_by_ds, "zstat")
bp_cluster <- processing_boxplot(final_cluster_by_ds, "diff")
bp_spatial <- processing_boxplot(final_spatial_by_ds, "diff")

bp_sample
bp_feature
bp_cluster
bp_spatial

ggplot2::ggsave("final_result/boxplot/bp_sample.pdf", bp_sample, width = 6, height = 4, dpi=300)
ggplot2::ggsave("final_result/boxplot/bp_feature.pdf", bp_feature, width = 6, height = 4, dpi=300)
ggplot2::ggsave("final_result/boxplot/bp_cluster.pdf", bp_cluster, width = 6, height = 4, dpi=300)
ggplot2::ggsave("final_result/boxplot/bp_spatial.pdf", bp_spatial, width = 6, height = 4, dpi=300)

# panel b

sample_cor <- read_csv("final_result/zstat/overall/sample_merged.csv")
feature_cor <- read_csv("final_result/zstat/overall/feature_merged.csv")
spatial_cor <- read_csv("final_result/zstat/overall/spatial_merged.csv")
cluster_cor <- read_csv("final_result/zstat/overall/cluster_merged.csv")

forest_plot_processing <- function(dataset){
  to_plot <-  dataset  %>%
    group_by(metric) %>%
    do(tidy(lm( value ~  model, data = .))) %>%
    rename(model = metric) %>%  
    relabel_predictors(c(
      'modelSRTsim(rf)' = "SRTsim(rf)",
      'modelscDesign3_gau(rf)' = "scDesign3_gau(rf)",
      'modelscDesign3_nb(rf)' = "scDesign3_nb(rf)",
      'modelscDesign3_poi(rf)' = "scDesign3_poi(rf)",
      modelSRTsim = "SRTsim",
      modelscDesign3_gau = "scDesign3_gau",
      modelscDesign3_nb = "scDesign3_nb",
      modelscDesign3_poi = "scDesign3_poi",
      '(Intercept)' = "scDesign2",
      modelsplatter = "splatter",
      modelSPARsim = "SPARsim",
      
      modelsymsim = "symsim",
      modelzinbwave = "zinbwave"
    ))
  
  g <- dwplot(to_plot,
              vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + theme_minimal()  + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1, fill=NA)) 
  
  return(to_plot)
}

fp_sample <- forest_plot_processing(sample_cor)
fp_feature <- forest_plot_processing(feature_cor)
fp_cluster <- forest_plot_processing(cluster_cor)
fp_spatial <- forest_plot_processing(spatial_cor)

ggplot2::ggsave("final_result/forestplot/fp_sample.pdf", fp_sample, width = 10, height = 6, dpi=300)
ggplot2::ggsave("final_result/forestplot/fp_feature.pdf", fp_feature, width = 10, height = 6, dpi=300)
ggplot2::ggsave("final_result/forestplot/fp_cluster.pdf", fp_cluster, width = 10, height = 6, dpi=300)
ggplot2::ggsave("final_result/forestplot/fp_spatial.pdf", fp_spatial, width = 10, height = 6, dpi=300)

