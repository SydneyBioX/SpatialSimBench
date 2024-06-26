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
library(readr) 
library(reticulate)

dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")
sim_model <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")

# Spot-level

column_mapping <- list(
  spfracZero = list(c("real" = "realCount", "sim" = "simCount")),
  spLibSize = list(c("real" = "libReal", "sim" = "libSim")),
  spTMM = list(c("real" = "realTMM", "sim" = "simTMM")),
  effLibSize = list(c("real", "sim")),
  combinedMetric = list(c("real", "sim")),
  spCor_spearman_new = list(c("real" = "real", "sim" = "sim")),
  sampleCor_pearson_new = list(c("real" = "real", "sim" = "sim")),
  sp_meanVar = list(
    c("real" = "realVar", "sim" = "simVar", "name" = "sp_meanVar_Var"),
    c("real" = "scaled_realVar", "sim" = "scaled_simVar", "name" = "sp_meanVar_scaledVar"),
    c("real" = "realMean", "sim" = "simMean", "name" = "sp_meanVar_Mean"),
    c("real" = "scaled_realMean", "sim" = "scaled_simMean", "name" = "sp_meanVar_scaledMean")
  )
)

# Create a list of all unique metric names, including submetrics with descriptive names
all_metrics <- unlist(lapply(names(column_mapping), function(metric) {
  column_sets <- column_mapping[[metric]]
  if (length(column_sets) > 1) {
    sapply(column_sets, function(set) set["name"])
  } else {
    metric
  }
}))

zstat_matrix <- matrix(0, nrow = length(sim_model), ncol = length(all_metrics))
rownames(zstat_matrix) <- sim_model
colnames(zstat_matrix) <- all_metrics

zstat_list <- list()
app <- "overall"
for (sim_index in seq_along(sim_model)) {
  sim <- sim_model[sim_index]
  sim_list <- list()
  message(paste0(sim, " start"))
  metric_counter <- 1  # Counter for column index in zstat_matrix
  
  for (metric in names(column_mapping)) {
    metric_list <- list()
    # Check if the current metric is the special case
    message(paste0(sim, "-", metric))
    if (metric == "effLibSize") {
      zstats <- numeric()  #
      
      for (dataset in dataset_names) {
        dataset_name <- sub("\\.rds$", "", dataset)
        
        # Read libSize and TMMfactor data files
        libSize_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/spLibSize.csv")
        TMMfactor_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/spTMM.csv")
        
        libSize <- fread(libSize_path)
        TMMfactor <- fread(TMMfactor_path)
        
        # Perform the computation for effLibSize
        effLibSimNew <- libSize$libSim * TMMfactor$simTMM
        effLibRealNew <- libSize$libReal * TMMfactor$realTMM
        
        # Perform kde.test and store zstat
        zstat <- ks::kde.test(x1 = effLibRealNew, x2 = effLibSimNew)$zstat
        zstats <- c(zstats, zstat)
        
      }
      
      # Calculate the average zstat for effLibSize
      metric_list <- zstats
      
      zstat_matrix[sim_index, metric_counter] <- mean(zstats)
      metric_counter <- metric_counter + 1
    } else if (metric == "combinedMetric"){
      zstats <- numeric()  # Temporary storage for zstats for averaging
      
      for (dataset in dataset_names) {
        dataset_name <- sub("\\.rds$", "", dataset)
        
        # Read the datasets for 'spfracZero' and 'spLibSize'
        fracZero_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/spfracZero.csv")
        libSize_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/spLibSize.csv")
        
        fracZero <- fread(fracZero_path)
        libSize <- fread(libSize_path)
        
        # Construct two-dimensional dataframes
        real_df <- data.frame(libSize$libReal, fracZero$realCount)
        sim_df <- data.frame(libSize$libSim, fracZero$simCount)
        
        # Perform kde.test and store zstat
        zstat <- ks::kde.test(x1 = as.matrix(real_df), x2 = as.matrix(sim_df))$zstat
        zstats <- c(zstats, zstat)
      }
      
      metric_list <- zstats
      
      # Calculate the average zstat for the special metric 
      zstat_matrix[sim_index, which(all_metrics == "combinedMetric")] <- mean(zstats)
      metric_counter <- metric_counter + 1
    } else {
      
      for (column_set in column_mapping[[metric]]) {
        zstats <- numeric()  # Temporary storage for zstats for averaging
        
        
        for (dataset in dataset_names) {
          dataset_name <- sub("\\.rds$", "", dataset)
          file_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/", metric, ".csv")
          
          data <- fread(file_path)
          
          zstat <- ks::kde.test(x1 = data[[column_set["real"]]], x2 = data[[column_set["sim"]]])$zstat
          zstats <- c(zstats, zstat)
        }
        
        metric_list <- zstats
        
        zstat_matrix[sim_index, metric_counter] <- mean(zstats)
        metric_counter <- metric_counter + 1
      }
    }
    sim_list[[metric]] <- metric_list
  }
  zstat_list[[sim]] <- sim_list
}

zstat_list_preprocessed <- lapply(zstat_list, function(sim_model) {
  lapply(sim_model, function(metric) {
    # Check each value in the metric, if less than 0, set to 0
    sapply(metric, function(value) {
      if (value < 0) 0 else value
    })
  })
})

spot_zstat <- zstat_list_preprocessed

# gene-level

column_mapping <- list(
  fcFracZero = list(c("real" = "realCount", "sim" = "simCount")),
  genePearson_new = list(c("real" = "real", "sim" = "sim")),
  geneSpearman_new = list(c("real" = "real", "sim" = "sim")),
  
  feature_meanVar = list(
    c("real" = "realVar", "sim" = "simVar", "name" = "ft_Var"),
    c("real" = "scaled_realVar", "sim" = "scaled_simVar", "name" = "ft_scaledVar"),
    c("real" = "realMean", "sim" = "simMean", "name" = "ft_mean"),
    c("real" = "scaled_realMean", "sim" = "scaled_simMean", "name" = "ft_scaledMean")
  ),
  
  mean_val = list(c("real", "sim")),
  mean_val_scale = list(c("real", "sim")),
  mean_fraczero  = list(c("real", "sim"))
  
)

# Create a list of all unique metric names, including submetrics with descriptive names
all_metrics <- unlist(lapply(names(column_mapping), function(metric) {
  column_sets <- column_mapping[[metric]]
  if (length(column_sets) > 1) {
    sapply(column_sets, function(set) set["name"])
  } else {
    metric
  }
}))

zstat_matrix <- matrix(0, nrow = length(sim_model), ncol = length(all_metrics))
rownames(zstat_matrix) <- sim_model
colnames(zstat_matrix) <- all_metrics

zstat_list_gene <- list()
sim_list <- list()
app <- "overall"
for (sim_index in seq_along(sim_model)) {
  sim <- sim_model[sim_index]
  message(paste0(sim, " start"))
  metric_counter <- 1  
  
  for (metric in names(column_mapping)) {
    message(paste0(sim, "-", metric))
    if (metric == "mean_val") {
      zstats <- numeric() 
      
      for (dataset in dataset_names) {
        dataset_name <- sub("\\.rds$", "", dataset)
        
        # Read the datasets for 'spfracZero' and 'spLibSize'
        feature_meanVar_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/feature_meanVar.csv")
        feature_meanVar <- fread(feature_meanVar_path)
        
        # Construct two-dimensional dataframes
        real_df <- data.frame(feature_meanVar$realMean, feature_meanVar$realVar)
        sim_df <- data.frame(feature_meanVar$simMean, feature_meanVar$simVar)
        
        # Perform kde.test and store zstat
        zstat <- ks::kde.test(x1 = as.matrix(real_df), x2 = as.matrix(sim_df))$zstat
        zstats <- c(zstats, zstat)
      }
      
      metric_list <- zstats
      
      # Calculate the average zstat
      zstat_matrix[sim_index, which(all_metrics == "mean_val")] <- mean(zstats)
      # zstat_matrix[sim_index, metric_counter] <- mean(zstats)
      metric_counter <- metric_counter + 1
    } else if (metric == "mean_val_scale"){
      zstats <- numeric()  # Temporary storage for zstats for averaging
      
      for (dataset in dataset_names) {
        dataset_name <- sub("\\.rds$", "", dataset)
        
        # Read the datasets for 'spfracZero' and 'spLibSize'
        feature_meanVar_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/feature_meanVar.csv")
        
        feature_meanVar <- fread(feature_meanVar_path)
        
        
        # Construct two-dimensional dataframes
        real_df <- data.frame(feature_meanVar$scaled_realMean, feature_meanVar$scaled_realVar)
        sim_df <- data.frame(feature_meanVar$scaled_simMean, feature_meanVar$scaled_simVar)
        
        # Perform kde.test and store zstat
        zstat <- ks::kde.test(x1 = as.matrix(real_df), x2 = as.matrix(sim_df))$zstat
        zstats <- c(zstats, zstat)
      }
      metric_list <- zstats
      
      # Calculate the average zstat for the special metric
      zstat_matrix[sim_index, which(all_metrics == "mean_val_scale")] <- mean(zstats)
      metric_counter <- metric_counter + 1
    } else if (metric == "mean_fraczero"){
      
      zstats <- numeric()  # Temporary storage for zstats for averaging
      
      for (dataset in dataset_names) {
        dataset_name <- sub("\\.rds$", "", dataset)
        
        # Read the datasets for 'spfracZero' and 'spLibSize'
        feature_meanVar_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/feature_meanVar.csv")
        fcFracZero_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/fcFracZero.csv")
        
        feature_meanVar <- fread(feature_meanVar_path)
        fcFracZero <- fread(fcFracZero_path)
        
        # Construct two-dimensional dataframes
        real_df <- data.frame(feature_meanVar$realMean, fcFracZero$realCount)
        sim_df <- data.frame(feature_meanVar$simMean, fcFracZero$simCount)
        
        # Perform kde.test and store zstat
        zstat <- ks::kde.test(x1 = as.matrix(real_df), x2 = as.matrix(sim_df))$zstat
        zstats <- c(zstats, zstat)
      }
      metric_list <- zstats
      
      # Calculate the average zstat for the special metric
      zstat_matrix[sim_index, which(all_metrics == "mean_fraczero")] <- mean(zstats)
      metric_counter <- metric_counter + 1
      
      
    } else {
      
      
      for (column_set in column_mapping[[metric]]) {
        zstats <- numeric()  # Temporary storage for zstats for averaging
        
        for (dataset in dataset_names) {
          dataset_name <- sub("\\.rds$", "", dataset)
          file_path <- paste0("output/", app, "/", sim, "/", dataset_name, "/", metric, ".csv")
          
          data <- fread(file_path)
          
          # Perform kde.test and store zstat
          zstat <- ks::kde.test(x1 = data[[column_set["real"]]], x2 = data[[column_set["sim"]]])$zstat
          zstats <- c(zstats, zstat)
        }
        metric_list <- zstats
        
        # Calculate the average zstat for this specific metric set
        zstat_matrix[sim_index, metric_counter] <- mean(zstats)
        metric_counter <- metric_counter + 1
      }
      
    }
    sim_list[[metric]] <- metric_list
  }
  zstat_list_gene[[sim]] <- sim_list
}


zstat_list_preprocessed <- lapply(zstat_list_gene, function(sim_model) {
  lapply(sim_model, function(metric) {
    # Check each value in the metric, if less than 0, set to 0
    sapply(metric, function(value) {
      if (value < 0) 0 else value
    })
  })
})

gene_zstat <- zstat_list_preprocessed

# spatial level

use_condaenv("calibenv", required = TRUE)
source_python("S8spatial.py")

dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")
sim_model <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")

for (dataset_name in dataset_names) {
  message(paste0("start ",dataset_name))
  real_sce_path <- paste0("data/final/", dataset_name, ".rds")
  
  # Check if the real dataset file exists
  if (file.exists(real_sce_path)) {
    real_sce <- readRDS(real_sce_path)
    
    # Loop over each simulation model
    for (sim_model in sim_models) {
      # Path to the simulated dataset
      message(paste0("start ",dataset_name, "-", sim_model))
      sim_sce_path <- paste0("output/overall/", sim_model, "/", dataset_name, "/final_sim.rds")
      
      # Check if the simulated dataset file exists
      if (file.exists(sim_sce_path)) {
        # Read the simulated dataset
        sim_sce <- readRDS(sim_sce_path)
        
        df <- data.frame(x = real_sce@colData$row, 
                         y = real_sce@colData$col, 
                         real = real_sce@colData$spatial.cluster,
                         sim = sim_sce@colData$spatial.cluster)
        
        output_path <- paste0("output/overall/", sim_model, "/", dataset_name, "/spatial.png")
        
        final_spatial(df, output_path)
      }
    }
  }
}


metrics <- c("TM", "NWE", "CSM")
dataset_averages_df <- as.data.frame(matrix(nrow = length(dataset_names), ncol = length(sim_models)))
rownames(dataset_averages_df) <- dataset_names
colnames(dataset_averages_df) <- sim_models

for (model in sim_models) {
  for (dataset in dataset_names) {
    file_path <- paste0("output/overall/", model, "/", dataset, "/spatial_metrics.csv")
    if(file.exists(file_path)) {
      # Read the CSV file
      df <- read_csv(file_path)
      # Extract the values for the metrics and calculate their average
      metric_values <- df$Value[match(metrics, df$Metric)]
      average_metric <- mean(metric_values, na.rm = TRUE)
      dataset_averages_df[dataset, model] <- average_metric
    }
  }
}

final_spatial_by_ds <- t(dataset_averages_df)
# saveRDS(final_spatial_by_ds, "final_spatial_by_ds.rds")

# Define the Human_list if not already defined
Human_list <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE") 

# Multiple choices
approach <- c("overall")
sim_models <- c("scDesign2", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SPARsim", "splatter", "SRTsim", "symsim", "zinbwave", "scDesign3_gau_rf", "scDesign3_nb_rf", "scDesign3_poi_rf", "SRTsim_rf")
dataset_names <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")


# Loop through each approach, sim_model, and dataset_name
for (appr in approach) {
  for (sim_model in sim_models) {
    for (dataset_name in dataset_names) {
      # Check if dataset_name is in Human_list
      cat("Start for", appr, sim_model, dataset_name, "\n")
      
      if (dataset_name %in% Human_list) {
        sc_species <- "Homo sapiens"
      } else {
        sc_species <- "Mus musculus"
      }
      
      # Construct the file path based on current approach, sim_model, and dataset_name
      file_path <- paste0("output/", appr, "/", sim_model, "/", dataset_name, "/final_sim_new.rds")
      
      # Attempt to read the RDS file; proceed only if the file exists
      if (file.exists(file_path)) {
        sim_sce <- readRDS(file_path)
        
        # Proceed with your operations...
        sim_sce <- scater::logNormCounts(sim_sce)
        sampleID <- rep("sample1", ncol(sim_sce))
        feat_types <- c("celltype_interaction")
        sc_type <- "spatial_t"
        prob_matrix <- t(sim_sce@metadata$celltype_prop)
        
        ct_interaction <- scFeatures::scFeatures(assay(sim_sce, "logcounts"),
                                                 sample = sampleID,
                                                 spatialCoords = list(colData(sim_sce)$row, colData(sim_sce)$col),
                                                 feature_types = feat_types,
                                                 type = sc_type,
                                                 species = sc_species,
                                                 spotProbability = prob_matrix)
        
        scfeatures_sim_path <- paste0("output/", appr, "/", sim_model, "/", dataset_name, "/scfeatures_sim.rds")
        if (file.exists(scfeatures_sim_path)) {
          scfeatures_sim <- readRDS(scfeatures_sim_path)
        } else {
          scfeatures_sim <- list() 
        }
        
        scfeatures_sim$celltype_interaction <- ct_interaction$celltype_interaction
        saveRDS(scfeatures_sim, scfeatures_sim_path)
        
        cat("End with", appr, sim_model, dataset_name, "\n")
      }
    }
  }
}

