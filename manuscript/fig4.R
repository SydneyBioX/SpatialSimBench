library(dplyr)
library(magrittr)
library(janitor)
library(tidyr)
library(stringr)
library(funkyheatmap)
library(ggthemr)
library(readr)

column_info_detailed <- tribble(
  ~id,               ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "method",           "method",             NA,                         "text",       NA,          list(hjust = 0, width = 3),
  "adaptor",           "method",             "Adaptor",                         "text",       NA,          list(hjust = 0, width = 1),
  
  "overall_summary",     "summary",       "Overall",    "bar",  "overall",  list(legend=TRUE),
  "data_pro_summary",     "summary",       "data proporties",    "bar",  "palette1",  list(legend=TRUE),
  "spatial_task_summary",     "summary",       "Spatial Task",    "bar",  "palette2",  list(legend=TRUE),
  "scalability_summary",     "summary",       "Scalability",    "bar",  "palette3",  list(legend=TRUE),
  
  
  "spfracZero",     "data_pro_spot",       "fracZero",    "funkyrect",  "palette1",  list(legend=TRUE),
  "spLibSize",      "data_pro_spot",       "LibSize",    "funkyrect",  "palette1",  list(legend=FALSE),
  "spEffLibSize",  "data_pro_spot",       "EffLibSize",    "funkyrect",  "palette1",  list(legend=TRUE),
  "spTMM",          "data_pro_spot",       "TMM",    "funkyrect",  "palette1",  list(legend=TRUE),
  "spScaledVar",    "data_pro_spot",       "ScaledVar",    "funkyrect",  "palette1",  list(legend=FALSE),
  "spScaledMean",   "data_pro_spot",       "ScaledMean",    "funkyrect",  "palette1",  list(legend=TRUE),
  "splibsizeFracezero",  "data_pro_spot",       "libFraczero(higher-order)",    "funkyrect",  "palette1",  list(legend=FALSE),
  "spCorPearson",    "data_pro_spot",       "Pearson(higher-order)",    "funkyrect",  "palette1",  list(legend=TRUE),
  
  "ftFracZero",     "data_pro_gene",       "fracZero",    "funkyrect",  "palette1",  list(legend=TRUE),
  "ftScaledVar",      "data_pro_gene",       "ScaledVar",    "funkyrect",  "palette1",  list(legend=FALSE),
  "ftScaledMean",  "data_pro_gene",       "ScaledMean",    "funkyrect",  "palette1",  list(legend=TRUE),
  "ftPearson",          "data_pro_gene",       "Pearson(higher-order)",    "funkyrect",  "palette1",  list(legend=TRUE),
  "ftScaleMeanVal",    "data_pro_gene",       "MeanVal(higher-order)",    "funkyrect",  "palette1",  list(legend=FALSE),
  "ftMeanFraczero",   "data_pro_gene",       "MeanFracZero(higher-order)",    "funkyrect",  "palette1",  list(legend=TRUE),
  
  "TM",    "data_pro_spatial",       "TM",    "funkyrect",  "palette1",  list(legend=TRUE),
  "NWE",     "data_pro_spatial",       "NWE",    "funkyrect",  "palette1",  list(legend=TRUE),
  "CSM",      "data_pro_spatial",       "CSM",    "funkyrect",  "palette1",  list(legend=FALSE),
  "celltype_interaction",  "data_pro_spatial",       "celltype_interaction",    "funkyrect",  "palette1",  list(legend=TRUE),
  "L_stats",          "data_pro_spatial",       "L stats",    "funkyrect",  "palette1",  list(legend=TRUE),
  "nn_correlation",    "data_pro_spatial",       "nn_correlation",    "funkyrect",  "palette1",  list(legend=FALSE),
  "morans_I",   "data_pro_spatial",       "Moran's I",    "funkyrect",  "palette1",  list(legend=TRUE),
  
  "ARI",    "spatial_task_clustering",       "ARI",    "funkyrect",  "palette2",  list(legend=TRUE),
  "NMI",     "spatial_task_clustering",       "NMI",    "funkyrect",  "palette2",  list(legend=TRUE),
  "recall_svg",      "spatial_task_svg",       "recall",    "funkyrect",  "palette2",  list(legend=FALSE),
  "precision_svg",  "spatial_task_svg",       "precision",    "funkyrect",  "palette2",  list(legend=TRUE),
  "rmse_ctprob",          "spatial_task_den",       "RMSE",    "funkyrect",  "palette2",  list(legend=TRUE),
  "jsd_ctprob",    "spatial_task_den",       "JSD",    "funkyrect",  "palette2",  list(legend=FALSE),
  "mantel_cor",   "spatial_task_cor",       "Mantel stat",    "funkyrect",  "palette2",  list(legend=TRUE),
  "cosine_cor",   "spatial_task_cor",       "cosine similarity",    "funkyrect",  "palette2",  list(legend=TRUE),
  
  "pred_time_cells200_features200",    "scalability_time",       "200 x 200",    "rect",  "palette3",  list(legend=TRUE),
  "pred_time_cells200_features200",    "scalability_time",       NA,    "text",  "white6black4",  list(label="predstring_time_cells200_features200", overlay=TRUE, size = 3, scale = FALSE),
  "pred_time_cells3000_features500",     "scalability_time",       "3k x 500",    "rect",  "palette3",  list(legend=TRUE),
  "pred_time_cells3000_features500",     "scalability_time",       NA,    "text",  "white6black4",  list(label="predstring_time_cells3000_features500", overlay=TRUE, size = 3, scale = FALSE),
  "pred_time_cells5000_features1000",      "scalability_time",       "5k x 1k",    "rect",  "palette3",  list(legend=FALSE),
  "pred_time_cells5000_features1000",      "scalability_time",       NA,    "text",  "white6black4",  list(label="predstring_time_cells5000_features1000", overlay=TRUE, size = 3, scale = FALSE),
  
  "pred_memory_cells200_features200",  "scalability_memory",       "200 x 200",    "rect",  "palette3",  list(legend=TRUE),
  "pred_memory_cells200_features200",  "scalability_memory",       NA,    "text",  "white6black4",  list(label="predstring_memory_cells200_features200", overlay=TRUE, size = 3, scale = FALSE),
  "pred_memory_cells3000_features500",          "scalability_memory",       "3k x 500",    "rect",  "palette3",  list(legend=TRUE),
  "pred_memory_cells3000_features500",          "scalability_memory",       NA,    "text",  "white6black4",  list(label="predstring_memory_cells3000_features500", overlay=TRUE, size = 3, scale = FALSE),
  "pred_memory_cells5000_features1000",    "scalability_memory",       "5k x 1k",    "rect",  "palette3",  list(legend=FALSE),
  "pred_memory_cells5000_features1000",    "scalability_memory",      NA,    "text",  "white6black4",  list(label="predstring_memory_cells5000_features1000", overlay=TRUE, size = 3, scale = FALSE),
  
)


# Define column groups
column_groups <- tribble( # tribble_start
  ~Experiment, ~Category,  ~group,         ~palette,
  "Method", "",  "method",      "overall",
  "Summary", "Rank-based overall", "summary",       "overall",
  "Data proporties", "Spot-level", "data_pro_spot",       "palette1",
  "Data proporties", "Gene-level", "data_pro_gene",       "palette1",
  "Data proporties", "Spatial-level", "data_pro_spatial",       "palette1",
  
  "Spatial task", "Cluster", "spatial_task_clustering",       "palette2",
  "Spatial task", "SVG", "spatial_task_svg",       "palette2",
  "Spatial task", "Deconvolute", "spatial_task_den",       "palette2",
  "Spatial task", "Crosscorrelate", "spatial_task_cor",       "palette2",
  
  
  "Scalability", "Predicted time", "scalability_time",       "palette3",
  "Scalability", "Predicted memory", "scalability_memory",       "palette3",
) 


generateColorPalette <- function(baseColor, numColors, rev=FALSE) {
  # Create a gradient from the base color to a lighter version
  colors <- colorRampPalette(c(baseColor, "white"))(numColors)
  if (rev) {
    colors=rev(colors)
  }
  return(colors)
}

palettes <- tribble(
  ~palette,             ~colours,
  #c9c7b9 recapitulation gene expression
  "palette1",           generateColorPalette("#6ABAC4", 100),
  #e07a5f tcga prediction
  "palette2",           generateColorPalette("#FA8072", 100),
  #3d405b clinical prediction
  "palette3",           generateColorPalette("#FFE5B4", 100),
  #81b29a usability
  "palette4",           generateColorPalette("#c9c7b9", 100),
)

palettes_new <- dynbenchmark_data$palettes
palettes <- rbind(palettes, palettes_new[1,])
palettes <- rbind(palettes, palettes_new[7,])

data<- read_csv("output/data_sorted_final_NOTCHANGE.csv")

data$...1 <- NULL
data$adaptor <- c("Yes","No", "Yes", "No", "Yes","Yes","No","Yes","Yes","Yes","No","Yes","Yes")

rank_based_summary <- function(subset_data) {
  ranked_data <- apply(subset_data, 2, function(x) rank(x))  # Rank such that larger values get higher ranks
  rowMeans(ranked_data)
}

g <- funky_heatmap(
  data = data,
  column_info = column_info_detailed,
  column_groups = column_groups,
  palettes = palettes,
  # legends = legends_score,
  col_annot_offset = 3.2
)

g