reorder_matrix <- function(count_matrix){
  counts_df_sim <- data.frame(count_matrix)
  colnames(counts_df_sim) <- sub("^X", "", colnames(counts_df_sim))
  
  # optinal filter spot
  # counts_df_sim <- counts_df_sim[, colnames(counts_df_sim) %in% spots_cluster_2]
  
  other_colnames <- colnames(counts_df_sim)[colnames(counts_df_sim) != "Gene"]
  sorted_colnames <- other_colnames[order(as.numeric(gsub("x.*", "", other_colnames)))]
  new_order <- c(sorted_colnames)
  counts_df_sim <- counts_df_sim[, new_order]
  counts_df_sim <- counts_df_sim[order(rownames(counts_df_sim)), ]
  counts_df_sim$Gene <- rownames(counts_df_sim)
  counts_df_sim <- counts_df_sim[, c("Gene", colnames(counts_df_sim)[-ncol(counts_df_sim)])]
  
  return(counts_df_sim)
}

reorder_matrix_real <- function(real_sce){
  counts_matrix <- assay(real_sce, "counts")
  
  # optinal filter spot
  # counts_new_matrix <- counts_matrix[, colnames(counts_matrix) %in% spots_cluster_2]
  # counts_df_real <- as.data.frame(counts_new_matrix)
  
  counts_df_real <- data.frame(Gene = rownames(counts_df_real), counts_df_real)
  colnames(counts_df_real) <- sub("^X", "", colnames(counts_df_real))
  other_colnames <- colnames(counts_df_real)[colnames(counts_df_real) != "Gene"]
  sorted_colnames <- other_colnames[order(as.numeric(gsub("x.*", "", other_colnames)))]
  new_order <- c("Gene", sorted_colnames)
  counts_df_real <- counts_df_real[, new_order]
  counts_df_real <- counts_df_real[order(rownames(counts_df_real)), ]
  return(counts_df_real)
}

create_sce <- function(data) {
  expr_matrix <- as.matrix(data[,-1])
  rownames(expr_matrix) <- data$Gene
  sce_object <- SingleCellExperiment(assays = list(counts = expr_matrix))
  return(sce_object)
}

compute_logcounts <- function(sce_object) {
  lib_sizes <- colSums(assay(sce_object, "counts"))
  norm_counts <- t( t(assay(sce_object, "counts")) / lib_sizes ) * 1e6
  logcounts <- log1p(norm_counts)
  assay(sce_object, "logcounts") <- logcounts
  return(sce_object)
}

simbench_result <- function(sce_real, sce_sim){
  if (is.null(sce_sim$celltype)){
    sce_sim$celltype <- "celltype"
    sce_real$celltype <- "celltype"
  }
  
  if (is.null(sce_sim$logcounts)){
    sce_real <- compute_logcounts(sce_real)
    sce_sim <- compute_logcounts(sce_sim)
  }
  
  
  parameter_spatial_result <- SimBench::eval_parameter(real = sce_real, sim = sce_sim, type = "raw" , method = "samplemethod")
  distribution_celltype_spatial <- parameter_spatial_result$raw_value$celltype$raw_value
  fig_spa <- draw_parameter_plot(distribution_celltype_spatial) 
  return(fig_spa)
}

computePCA <- function(sce_obj) {
  if (!"logcounts" %in% assayNames(sce_obj)) {
    sce_obj <- compute_logcounts(sce_obj)
  }
  logcounts_data <- assay(sce_obj, "logcounts")
  pca_result <- prcomp(t(logcounts_data))
  reducedDims(sce_obj)[["PCA"]] <- pca_result$x
  return(sce_obj)
}

rebuilt_sce <- function(sce_obj){
  counts_single <- as.matrix(assay(sce_obj, "counts"))
  col_data <- data.frame(
    spatial1 = data.frame(sce_obj@colData)$spatial1,
    spatial2 = data.frame(sce_obj@colData)$spatial2,
    # celltype = data.frame(sce_obj@colData)$celltype,
    row.names = rownames(data.frame(sce_obj@colData))
  )
  
  col_data$loc <- paste0(col_data$spatial1, "x", col_data$spatial2)
  col_data$cell <- rownames(col_data)
  name_mapping <- setNames(col_data$loc, col_data$cell)
  colnames(counts_single) <- name_mapping[colnames(counts_single)]
  counts_single <- data.frame(counts_single)
  colnames(counts_single) <- sub("^x", "", colnames(counts_single))
  colnames(counts_single) <- sub("^X", "", colnames(counts_single))
  # optinal, add gene
  
  #col_data_sce <- data.frame(spatial1 = col_data$spatial1, spatial2 = col_data$spatial2, celltype = col_data$celltype)
  
  col_data_sce <- data.frame(spatial1 = col_data$spatial1, spatial2 = col_data$spatial2)
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts_single),
    colData = col_data_sce
  )
  
  return(sce)
}

rebuilt <- function(sce_obj){
  
  counts_data <- assays(sce_obj)$counts
  col_data_subset <- colData(sce_obj)[, c("row", "col", "spatial.cluster")]
  new_sce <- SingleCellExperiment(assays = list(counts = counts_data),
                                  colData = col_data_subset)
  
  return(new_sce)

}

extract_name <- function(path) {
  name_with_extension <- basename(path)  # Extracts the filename from the path
  name_without_extension <- sub("\\.rds$", "", name_with_extension)  # Removes the .rds extension
  return(name_without_extension)
}
