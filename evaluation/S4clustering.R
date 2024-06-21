cluster_metric <- function(real_seurat, sim_seurat, counts_real, counts_sim){
  
  real_seurat_cluster  <- data.frame(realCluster = Idents(real_seurat))
  sim_seurat_cluster  <- data.frame(simCluster = Idents(sim_seurat))
  realCluster <- real_seurat_cluster$realCluster
  simCluster <- sim_seurat_cluster$simCluster
  
  
  # ARI & NMI 
  ARI <- ARI(real_seurat_cluster$realCluster, sim_seurat_cluster$simCluster)
  NMI <- NMI(real_seurat_cluster$realCluster, sim_seurat_cluster$simCluster)
  
  # Confusion matrix
  # labels <- unique(realCluster)
  
  # precision <- numeric(length(labels))
  # recall <- numeric(length(labels))
  # fdr <- numeric(length(labels))
  
  # for (i in seq_along(labels)) {
  #  i <- 0
  #  label <- labels[i]
  
  #  TP <- sum(realCluster == label & simCluster == label)
  #  FP <- sum(realCluster != label & simCluster == label)
  #  FN <- sum(realCluster == label & simCluster != label)
  
  #  precision[i] <- TP / (TP + FP)
  #  recall[i] <- TP / (TP + FN)
  #  fdr[i] <- FP / (TP + FP)
  #}
  
  #macro_precision <- mean(precision, na.rm=TRUE)
  #macro_recall <- mean(recall, na.rm=TRUE)
  #macro_fdr <- mean(fdr, na.rm=TRUE)
  
  # update
  
  # Convert the clusters to integer type for processing
  realCluster_int <- as.integer(realCluster) - 1
  simCluster_int <- as.integer(simCluster) - 1
  
  # Create an empty confusion matrix
  confusion_matrix <- matrix(0, nrow=length(unique(realCluster_int)), ncol=length(unique(simCluster_int)))
  rownames(confusion_matrix) <- sort(unique(realCluster_int))
  colnames(confusion_matrix) <- sort(unique(simCluster_int))
  
  # Fill the confusion matrix
  for (i in unique(realCluster_int)) {
    for (j in unique(simCluster_int)) {
      confusion_matrix[i+1, j+1] <- sum(realCluster_int == i & simCluster_int == j)
    }
  }
  
  # Calculate precision, recall, and FDR
  precision <- numeric(length(unique(realCluster)))
  recall <- numeric(length(unique(realCluster)))
  fdr <- numeric(length(unique(realCluster)))
  
  for (i in 1:nrow(confusion_matrix)) {
    TP <- confusion_matrix[i, i]  # True Positives are the diagonal elements
    FP <- sum(confusion_matrix[-i, i])  # False Positives are the sum of column i excluding the diagonal
    FN <- sum(confusion_matrix[i, -i])  # False Negatives are the sum of row i excluding the diagonal
    
    precision[i] <- TP / (TP + FP)
    recall[i] <- TP / (TP + FN)
    fdr[i] <- FP / (TP + FP)
  }
  
  macro_precision <- mean(precision, na.rm=TRUE)
  macro_recall <- mean(recall, na.rm=TRUE)
  macro_fdr <- mean(fdr, na.rm=TRUE)
  
  
  rownames(counts_real) <- counts_real$Gene
  counts_real$Gene <- NULL
  rownames(counts_sim) <- counts_sim$Gene
  counts_sim$Gene <- NULL
  
  sil_width_real <- silhouette(as.numeric(as.character(realCluster)), dist(t(as.matrix(counts_real))))
  avg_sil_width_real <- mean(sil_width_real[, 3])
  
  sil_width_sim <- silhouette(as.numeric(as.character(simCluster)), dist(t(as.matrix(counts_sim))))
  avg_sil_width_sim <- mean(sil_width_sim[, 3])
  
  
  
  df_result <- data.frame(
    Metric = c("ARI","NMI","Precision", "Recall", "FDR", "silWidthReal", "silWidthSim"),
    Value = c(ARI,NMI, macro_precision, macro_recall, macro_fdr, avg_sil_width_real, avg_sil_width_sim)
  )
  return(df_result)
  
}

cluster_metric_BayerSpace <- function(cluster_df, counts_real, counts_sim){
  
  counts_real <- data.frame(Gene = rownames(counts_real),counts_real)
  counts_sim <- data.frame(Gene = rownames(counts_sim),counts_sim)
  
  
  # ARI & NMI 
  ARI <- ARI(cluster_df$real, cluster_df$sim)
  NMI <- NMI(cluster_df$real, cluster_df$sim)
  
  realCluster <- cluster_df$real
  simCluster <- cluster_df$sim
  
  # Confusion matrix
  labels <- unique(realCluster)
  
  precision <- numeric(length(labels))
  recall <- numeric(length(labels))
  fdr <- numeric(length(labels))
  
  for (i in seq_along(labels)) {
    label <- labels[i]
    
    TP <- sum(realCluster == label & simCluster == label)
    FP <- sum(realCluster != label & simCluster == label)
    FN <- sum(realCluster == label & simCluster != label)
    
    precision[i] <- TP / (TP + FP)
    recall[i] <- TP / (TP + FN)
    fdr[i] <- FP / (TP + FP)
  }
  
  
  macro_precision <- mean(precision, na.rm=TRUE)
  macro_recall <- mean(recall, na.rm=TRUE)
  macro_fdr <- mean(fdr, na.rm=TRUE)
  
  rownames(counts_real) <- counts_real$Gene
  counts_real$Gene <- NULL
  rownames(counts_sim) <- counts_sim$Gene
  counts_sim$Gene <- NULL
  
  sil_width_real <- silhouette(as.numeric(as.character(realCluster)), dist(t(as.matrix(counts_real))))
  avg_sil_width_real <- mean(sil_width_real[, 3])
  
  sil_width_sim <- silhouette(as.numeric(as.character(simCluster)), dist(t(as.matrix(counts_sim))))
  avg_sil_width_sim <- mean(sil_width_sim[, 3])
  
  df_result <- data.frame(
    Metric = c("ARI","NMI","Precision", "Recall", "FDR", "silWidthReal", "silWidthSim"),
    Value = c(ARI,NMI, macro_precision, macro_recall, macro_fdr, avg_sil_width_real, avg_sil_width_sim)
  )
  return(df_result)
  
}

cluster_loc_csv <- function(loc_list, cluster_df){
  
  loc_data <- as.data.frame(do.call(rbind, strsplit(loc_list, "x")), stringsAsFactors = FALSE)
  colnames(loc_data) <- c("x", "y")
  loc_data$x <- as.numeric(loc_data$x)
  loc_data$y <- as.numeric(loc_data$y)
  cluster_loc <- cbind(loc_data, cluster_df)
  
  return(cluster_loc)
}

reclassify_simsce <- function(real_sce, sim_sce){
  test <- data.frame(loc = colnames(counts(real_sce)), real = real_sce$spatial.cluster, sim = sim_sce$spatial.cluster)
  
  matrix_counts <- matrix(0, nrow = 4, ncol = 4, 
                          dimnames = list(paste0("real", 1:4), paste0("sim", 1:4)))
  
  for(r in 1:4) {
    for(s in 1:4) {
      subset_data <- subset(test, real == r & sim == s)
      matrix_counts[r, s] <- length(unique(subset_data$loc))
    }
  }
  
  reclassification <- numeric(4)
  
  for(i in 1:4) {
    max_value <- max(matrix_counts)
    if (max_value == 0) break 
    
    indices <- which(matrix_counts == max_value, arr.ind = TRUE)
    real_index <- indices[1, 1]
    sim_index <- indices[1, 2]
    
    reclassification[sim_index] <- real_index
    
    matrix_counts[real_index, ] <- -Inf
    matrix_counts[, sim_index] <- -Inf
  }
  
  
  cluster_map <- setNames(reclassification, seq_along(reclassification))
  sim_sce$spatial.cluster <- cluster_map[sim_sce$spatial.cluster]
  
  return(sim_sce)
  
}

find_cluster_indices <- function(cluster_column) {
  unique_clusters <- sort(unique(cluster_column))
  conditions <- list()
  
  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]
    indices <- which(cluster_column == cluster)
    range_name <- sprintf("cluster_%s_column_index", LETTERS[i])
    assign(range_name, c(min(indices):max(indices)), envir = .GlobalEnv)
    conditions[[range_name]] <- get(range_name)
  }
  
  return(conditions)
}
