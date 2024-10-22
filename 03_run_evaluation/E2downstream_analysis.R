library('SPARK')
library(readr)
library(dplyr)
library(philentropy)
library(Metrics)
library(spots)
library(spdep)
library(vegan)
library(psych)
library(lsa)
library(CellChat)
library(patchwork)

# spatial clustering

sim_models <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")
dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")

dataset_averages_df <- as.data.frame(matrix(nrow = length(dataset_names), ncol = length(sim_models)))
rownames(dataset_averages_df) <- dataset_names
colnames(dataset_averages_df) <- sim_models

for (model in sim_models) {
  for (dataset in dataset_names) {
    file_path <- paste0("output/overall/", model, "/", dataset, "/df_result.csv")
    if(file.exists(file_path)) {
      df <- readr::read_csv(file_path)
      metric_values <- df$Value[match(c("ARI", "NMI", "Precision", "Recall"), df$Metric)]
      print(metric_values)
      average_metric <- mean(metric_values, na.rm = TRUE)
      dataset_averages_df[dataset, model] <- average_metric
    }
  }
}

# SVG

generate_svg_sparkx <- function(sp_sce){
  
  sp_count <- counts(sp_sce)
  info <- cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(sp_sce@colData),split="x"),"[",1)),
                           y=as.numeric(sapply(strsplit(rownames(sp_sce@colData),split="x"),"[",2)))
  
  rownames(info) <- colnames(sp_count)
  location <- as.matrix(info)
  
  mt_idx <- grep("mt-",rownames(sp_count))
  if(length(mt_idx)!=0){
    sp_count    <- sp_count[-mt_idx,]
  }
  
  sparkX <- sparkx(sp_count,location,numCores=1,option="mixture")
  
  return(sparkX)
}

for(model in sim_model) {
  for(dataset in dataset_names) {
    message(model, " start ", dataset)
    file_path <- paste0("~/spatial_simulationV2/output/overall/", model, "/", dataset, "/final_sim_new.rds")
    sp_sce <- readRDS(file_path)
    processed_data <- generate_svg_sparkx(sp_sce)
    save_path <- paste0("~/spatial_simulationV2/output/overall/", model, "/", dataset, "/final_sim_svg.rds")
    
    saveRDS(processed_data, save_path)
    message(model, " end ", dataset)
  }
}

for(dataset in dataset_names) {
  file_path <- paste0("~/spatial_simulationV2/data/final/", dataset, ".rds")
  
  sp_sce <- readRDS(file_path)
  processed_data <- generate_svg_sparkx(sp_sce)
  save_path <- paste0("~/spatial_simulationV2/SVG/", dataset, "_svg.rds")
  saveRDS(processed_data, save_path)
}

calculate_precision <- function(real_data_path, compared_data_path) {
  real_data <- readRDS(real_data_path)
  compared_data <- readRDS(compared_data_path)
  filtered_real_data <- real_data$res_mtest %>% 
    filter(adjustedPval < 0.05)
  
  filtered_compared_data <- compared_data$res_mtest %>% 
    filter(adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fp <- length(setdiff(row.names(filtered_compared_data), row.names(filtered_real_data)))
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  
  return(precision)
}

calculate_recall <- function(real_data_path, compared_data_path) {
  real_data <- readRDS(real_data_path)
  compared_data <- readRDS(compared_data_path)
  filtered_real_data <- real_data$res_mtest %>% 
    filter(adjustedPval < 0.05)
  filtered_compared_data <- compared_data$res_mtest %>% 
    filter(adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fn <- length(setdiff(row.names(filtered_real_data), row.names(filtered_compared_data)))
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  
  return(recall)
}

dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")
sim_model <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")

precision_matrix <- matrix(NA, nrow = length(sim_model), ncol = length(dataset_names), 
                           dimnames = list(sim_model, dataset_names))

recall_matrix <- matrix(NA, nrow = length(sim_model), ncol = length(dataset_names), 
                        dimnames = list(sim_model, dataset_names))
for (model in sim_model) {
  for (dataset in dataset_names) {
    real_data_path <- paste0("~/spatial_simulationV2/SVG/", dataset, "_svg.rds")
    compared_data_path <- paste0("~/spatial_simulationV2/output/overall/", model, "/", dataset, "/final_sim_svg.rds")
    
    precision <- calculate_precision(real_data_path, compared_data_path)
    recall <- calculate_recall(real_data_path, compared_data_path)
    
    precision_matrix[model, dataset] <- precision
    recall_matrix[model, dataset] <- recall
  }
}

# ct proportion

generate_jds <- function(real, sim) {
  common_row_names <- intersect(rownames(real), rownames(sim))
  
  real_common <- real[common_row_names, , drop = FALSE]
  sim_common <- sim[common_row_names, , drop = FALSE]
  
  jsd_values <- sapply(1:nrow(real_common), function(i) {
    x.count <- rbind(as.vector(real_common[i, ]), as.vector(sim_common[i, ]))
    JSD(x.count, est.prob = "empirical")
  })
  
  average_jsd <- mean(jsd_values)
  
  return(average_jsd)
}

generate_rmse <- function(real, sim) {
  common_row_names <- intersect(rownames(real), rownames(sim))
  
  real_common <- real[common_row_names, , drop = FALSE]
  sim_common <- sim[common_row_names, , drop = FALSE]
  
  rmse_values <- sapply(1:nrow(real_common), function(i) {
    rmse(as.vector(real_common[i, ]), as.vector(sim_common[i, ]))
  })
  
  average_rmse <- mean(rmse_values)
  
  return(average_rmse)
}

models <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")
datasets <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN" ,"MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")

jsd_matrix <- matrix(nrow = length(models), ncol = length(datasets), dimnames = list(models, datasets))
rmse_matrix <- matrix(nrow = length(models), ncol = length(datasets), dimnames = list(models, datasets))


for (model in models) {
  for (dataset in datasets) {
    real_path <- paste0("~/spatial_simulationV2/data/final/", dataset, ".rds")
    sim_path <- paste0("~/spatial_simulationV2/output/overall/", model, "/", dataset, "/final_sim_new.rds")
    
    real_dataset <- readRDS(real_path)
    sim_dataset <- readRDS(sim_path)
    
    real <- real_dataset@metadata$celltype_prop
    sim <- sim_dataset@metadata$celltype_prop
    
    jsd_matrix[model, dataset] <- generate_jds(real, sim)
    rmse_matrix[model, dataset] <- generate_rmse(real, sim)
  }
}

# saveRDS(jsd_matrix, "spatialTask/ct_prob/jsd_matrix.rds")
# saveRDS(rmse_matrix, "spatialTask/ct_prob/rmse_matrix.rds")

# Moran's I

generate_moransI <- function(data){
  
  W <-  data.frame( colData(data)$row, colData(data)$col )
  W <- 1/as.matrix(dist(W))
  diag(W) <- 0
  res <- BivariateMoransI( X =  t( logcounts(data) ) , W =  W)
  return(res$Morans.I)
}

generate_cosine <- function(real, sim){
  
  real_new <- real[!is.na(real) & !is.na(sim)]
  sim_new <- sim[!is.na(real) & !is.na(sim)]
  similarity <- cosine(as.textmatrix(cbind(as.vector(real_new), as.vector(sim_new))))
  return(mean(similarity))
  
}

dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN" ,"MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")

for(dataset in dataset_names) {
  
  file_path <- paste0("~/spatial_simulationV2/data/final/", dataset, ".rds")
  sp_sce <- readRDS(file_path)
  processed_data <- generate_moransI(sp_sce)
  save_path <- paste0("~/spatial_simulationV2/spatialTask/morans/", dataset, "_morans.rds")
  saveRDS(processed_data, save_path)
  
}

sim_model <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")
dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN" ,"MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")

for(model in sim_model) {
  for(dataset in dataset_names) {
    message(model, " start ", dataset)
    file_path <- paste0("~/spatial_simulationV2/output/overall/", model, "/", dataset, "/final_sim_new.rds")
    sp_sce <- readRDS(file_path)
    processed_data <- generate_moransI(sp_sce)
    save_path <- paste0("~/spatial_simulationV2/output/overall/", model, "/", dataset, "/final_sim_morans.rds")
    saveRDS(processed_data, save_path)
    message(model, " end ", dataset)
  }
}

models <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")
datasets <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN" ,"MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")

cosine_matrix <- matrix(nrow = length(models), ncol = length(datasets), dimnames = list(models, datasets))

for (model in models) {
  for (dataset in datasets) {
    real_path <- paste0("~/spatial_simulationV2/spatialTask/morans/", dataset, "_morans.rds")
    sim_path <- paste0("~/spatial_simulationV2/output/domain/", model, "/", dataset, "/final_sim_morans.rds")
    
    real_morans  <- readRDS(real_path)
    sim_morans <- readRDS(sim_path)
    
    cosine_matrix[model, dataset] <-  generate_cosine(real_morans, sim_morans)
  }
}