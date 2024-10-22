# SVG

library(nnSVG)
library(scran)
library(ggplot2)
library(SPARK)
library(MERINGUE)
library(Giotto)

run_sparkx <- function(sp_sce){
  
  sp_count <- counts(sp_sce)
  
  info <- cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(sp_sce@colData),split="x"),"[",1)), y=as.numeric(sapply(strsplit(rownames(sp_sce@colData),split="x"),"[",2)))
  rownames(info) <- colnames(sp_count)
  location <- as.matrix(info)
  mt_idx <- grep("mt-",rownames(sp_count))
  if(length(mt_idx)!=0){
    sp_count <- sp_count[-mt_idx,]
  }
  sparkX <- sparkx(sp_count,location,numCores=1,option="mixture")
  
  return(sum(sparkX$res_mtest <= 0.05))
}

run_nnSVG <- function(data){
  set.seed(123)
  spatial_coords <- data.frame(array_row = colData(data)$row, 
                               array_col = colData(data)$col)
  
  spatial_data <- SpatialExperiment::SpatialExperiment(
    assays = assays(data),
    rowData = rowData(data),
    colData = colData(data),
    spatialCoords = as.matrix(spatial_coords)
  )
  rowData(spatial_data)$gene_name <- rownames(spatial_data)
  spe <- filter_genes(spatial_data)
  spe <- computeLibraryFactors(spe)
  spe <- logNormCounts(spe)
  
  spe <- nnSVG(spe, n_threads = 1)
  
  return(sum(rowData(spe)$padj <= 0.05))
  
}

run_hvg <- function(data){
  
  data_seu <- Seurat::CreateSeuratObject( counts = logcounts(data))
  data_seu  <- Seurat::FindVariableFeatures( data_seu , selection.method = "vst" )
  hvg <- Seurat::HVFInfo(data_seu)
  hvg <- hvg[ order(hvg$variance.standardized, decreasing = T), ]
  
  return(hvg)
}

run_BOOSTGP <- function(data){
  sp_count <- counts(data)
  spatial_coords <- data.frame(x = colData(data)$row, 
                               y = colData(data)$col)
  rownames(spatial_coords) <- rownames(colData(data))
  res <- boost.gp(Y = Matrix::t(sp_count), loc = spatial_coords)
  return(sum(res$pval <= 0.05))
  
}

run_MERINGUE <- function(data) {
  set.seed(123)
  
  pos <- data.frame(
    x = data@colData$col,
    y = data@colData$row
  )
  
  cd <- assay(data, "counts")
  
  counts <- cleanCounts(counts = cd, 
                        min.reads = 100, 
                        min.lib.size = 100, 
                        plot = FALSE,
                        verbose = TRUE)
  
  pos <- pos[colnames(counts), ]
  mat <- normalizeCounts(counts = counts, 
                         log = FALSE,
                         verbose = TRUE)
  
  pcs.info <- prcomp(t(log10(as.matrix(mat) + 1)), center = TRUE)
  nPcs <- length(unique(data@colData$spatial.cluster))
  pcs <- pcs.info$x[, 1:nPcs]
  
  emb <- Rtsne::Rtsne(pcs,
                      is_distance = FALSE,
                      perplexity = 30,
                      num_threads = 1,
                      verbose = FALSE)$Y
  rownames(emb) <- rownames(pcs)
  
  k <- 30
  com <- getClusters(pcs, k, weight = TRUE)
  
  annot <- as.character(com)
  names(annot) <- names(com)
  
  dg <- getDifferentialGenes(as.matrix(mat), annot)
  dg.sig <- lapply(dg, function(x) {
    x <- x[x$p.adj < 0.05, ]
    x <- na.omit(x)
    x <- x[x$highest, ]
    rownames(x)
  })
  
  overall_gene_count <- sum(sapply(dg.sig, length))
  
  return(overall_gene_count)
}

run_giotto <- function( data, dir  ){
  
  my_giotto_object = createGiottoObject(raw_exprs = logcounts(data),
                                        spatial_locs =  data.frame( x = data$x_cord, y= data$y_cord  ))
  
  my_giotto_object@norm_expr <- my_giotto_object@raw_exprs
  
  
  my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k = 2)
  
  rank_spatialgenes = binSpect(my_giotto_object, bin_method = 'rank')
  
  return(sum(rank_spatialgenes$adj.p.value <= 0.05))
  
  
}

# Start running

outer_datasets <- c("HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", 
                    "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")


dataset_names <- c("SPARsim", "SRTsim", "SRTsim_rf", "scDesign3_poi_rf", "scDesign3_poi", 
                   "scDesign3_nb_rf", "scDesign3_nb", "scDesign3_gau_rf", "scDesign3_gau", 
                   "zinbwave", "symsim", "splatter", "scDesign2")

for (outer_dataset in outer_datasets) {
  

  message("Started processing outer dataset: ", outer_dataset)
  
  data <- readRDS(paste0("~/spatial_simulationV2/data/final/", outer_dataset, ".rds"))
  

  df_svg <- data.frame(Dataset = character(), sparkx = numeric(), nnSVG = numeric(),
                       hvg = numeric(), BOOSTGP = numeric(), stringsAsFactors = FALSE)
  
  test_sparkx <- tryCatch(run_sparkx(data), error = function(e) 0)
  test_nnSVG <- tryCatch(run_nnSVG(data), error = function(e) 0)
  test_hvg <- tryCatch(run_hvg(data), error = function(e) 0)
  
  if (is.data.frame(test_hvg)) {
    boundary_mean <- unname(summary(test_hvg$mean)[3])
    boundary_var <- unname(summary(test_hvg$var)[3])
    test_hvg_result <- sum(test_hvg$mean < boundary_mean & test_hvg$variance < boundary_var)
  } else {
    test_hvg_result <- 0
  }
  
  test_giotto <- tryCatch(run_giotto(data), error = function(e) 0)

  df_svg <- rbind(df_svg, data.frame(
    Dataset = "real", 
    sparkx = test_sparkx, 
    nnSVG = test_nnSVG, 
    hvg = test_hvg_result, 
    giotto = test_giotto
  ))
  
  for (dataset_name in dataset_names) {
    
    message("Started processing: ", dataset_name, " for outer dataset: ", outer_dataset)
    file_path <- paste0("~/spatial_simulationV2/output/overall/", dataset_name, "/", outer_dataset, "/final_sim.rds")
    
    data <- readRDS(file_path)
    
    test_sparkx <- tryCatch(run_sparkx(data), error = function(e) 0)
    test_nnSVG <- tryCatch(run_nnSVG(data), error = function(e) 0)
    test_hvg <- tryCatch(run_hvg(data), error = function(e) 0)
    
    if (is.data.frame(test_hvg)) {
      boundary_mean <- unname(summary(test_hvg$mean)[3])
      boundary_var <- unname(summary(test_hvg$var)[3])
      test_hvg_result <- sum(test_hvg$mean < boundary_mean & test_hvg$variance < boundary_var)
    } else {
      test_hvg_result <- 0
    }
    
    test_giotto <- tryCatch(run_giotto(data), error = function(e) 0)
    
    df_svg <- rbind(df_svg, data.frame(
      Dataset = dataset_name, 
      sparkx = test_sparkx, 
      nnSVG = test_nnSVG, 
      hvg = test_hvg_result, 
      giotto = test_giotto
    ))
    
    message("Finished processing: ", dataset_name, " for outer dataset: ", outer_dataset)
  }
  
  output_file <- paste0("~/spatial_simulationV2/spatialTask/SVG/revision/", outer_dataset, "_df_svg.csv")
  write.csv(df_svg, output_file, row.names = FALSE)
  
  message("Finished processing outer dataset: ", outer_dataset)
}

# spatial domain detection

library(BASS)
library(SingleCellExperiment)
library(SpatialPCA)
library(ggplot2)
library(PRECAST)
library(reticulate)
library(Seurat)
library(BayesSpace)

run_BASS <- function(sce_object, R) {
  count_matrix <- assay(sce_object, "counts")
  colnames(count_matrix) <- colnames(sce_object)
  
  xy <- data.frame(x = sce_object$col, y = sce_object$row)
  rownames(xy) <- colnames(sce_object)
  
  test_a <- list(count_matrix)
  test_b <- list(xy)
  
  # R <- length(unique(sce_object$manual_anno))
  
  BASS <- createBASSObject(test_a, test_b, C = 20, R = R, beta_method = "fix")
  BASS <- BASS.preprocess(BASS)
  BASS <- BASS.run(BASS)
  BASS <- BASS.postprocess(BASS, adjustLS = FALSE)
  
  df_output <- data.frame(slot1 = BASS@xy, slot2 = BASS@results$z)
  colnames(df_output)[ncol(df_output)] <- "spatial_cluster_BASS"
  rownames(df_output) <- paste(df_output$slot1.x, df_output$slot1.y, sep = "x")
  
  return(df_output)
}

run_SpatialPCA <- function(sce_object, R) {
  
  rawcount <- assay(sce_object, "counts")
  colnames(rawcount) <- colnames(sce_object)
  
  location <- as.matrix(data.frame(x = sce_object$col, y = sce_object$row))
  rownames(location) <- colnames(sce_object)
  
  ST <- CreateSpatialPCAObject(counts = rawcount, location = location, project = "SpatialPCA", 
                               gene.type = "spatial", sparkversion = "spark", 
                               gene.number = 3000, customGenelist = NULL, 
                               min.loctions = 20, min.features = 20)
  
  ST <- SpatialPCA_buildKernel(ST, kerneltype = "gaussian", bandwidthtype = "SJ")
  ST <- SpatialPCA_EstimateLoading(ST, fast = FALSE, SpatialPCnum = 20)
  ST <- SpatialPCA_SpatialPCs(ST, fast = FALSE)
  
  # R <- length(unique(sce_object$manual_anno))
  
  clusterlabel <- walktrap_clustering(R, ST@SpatialPCs, round(sqrt(dim(ST@location)[1])))
  clusterlabel_refine <- refine_cluster_10x(clusterlabel, ST@location, shape = "square")
  
  df_output <- data.frame(slot1 = ST@location, slot2 = clusterlabel_refine)
  colnames(df_output)[ncol(df_output)] <- "spatial_cluster_SpatialPCA"
  rownames(df_output) <- paste(df_output$slot1.x, df_output$slot1.y, sep = "x")
  
  return(df_output)
}


run_PRECAST_DESC <- function(sce_object, base, R) {
  
  rownames(sce_object) <- gsub("_", "-", rownames(sce_object))
  seurat_object <- as.Seurat(sce_object, counts = "counts", data = "logcounts")
  meta_data <- as.data.frame(colData(sce_object))
  rownames(meta_data) <- colnames(seurat_object)
  seurat_object <- AddMetaData(seurat_object, metadata = meta_data)
  
  preobj <- CreatePRECASTObject(seuList = list(seurat_object), selectGenesMethod = "HVGs", gene.number = 100)
  
  PRECASTObj <- AddAdjList(preobj, platform = base)
  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)
  
  # R <- length(unique(sce_object$manual_anno))
  PRECASTObj <- PRECAST(PRECASTObj, K = R)
  
  resList <- PRECASTObj@resList
  PRECASTObj <- SelectModel(PRECASTObj)
  
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  colnames(posList[[1]]) <- c("y", "x")
  
  seu_drsc <- DR.SC::DR.SC(PRECASTObj@seulist[[1]], K = R, verbose = T)
  
  df_output <- data.frame(
    slot1 = posList[[1]],
    spatial_cluster_PRECAST = factor(unlist(PRECASTObj@resList$cluster))
  )
  
  spatial_cluster_DESC <- seu_drsc$spatial.drsc.cluster
  matched_cluster_desc <- spatial_cluster_DESC[rownames(df_output)]
  df_output$spatial_cluster_DESC <- as.factor(matched_cluster_desc)
  rownames(df_output) <- paste(df_output$slot1.x, df_output$slot1.y, sep = "x")
  
  return(df_output)
}

run_leiden <- function(sce_object, R){
  
  seurat_object <- as.Seurat(sce_object, counts = "counts", data = "logcounts")
  meta_data <- as.data.frame(colData(sce_object))
  rownames(meta_data) <- colnames(seurat_object)
  seurat_object <- AddMetaData(seurat_object, metadata = meta_data)
  
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- FindNeighbors(seurat_object, dims = 1:10)  
  
  min_resolution <- 0.1
  max_resolution <- 1.0
  tolerance <- 1e-2
  
  target_clusters <- R
  current_resolution <- 1.0
  
  seurat_object <- FindClusters(seurat_object, algorithm = 4, resolution = current_resolution)
  current_clusters <- R
  
  while (abs(current_clusters - target_clusters) > tolerance) {
    if (current_clusters > target_clusters) {
      max_resolution <- current_resolution
    } else {
      min_resolution <- current_resolution
    }
    
    current_resolution <- (min_resolution + max_resolution) / 2
    seurat_object <- FindClusters(seurat_object, algorithm = 4, resolution = current_resolution)
    current_clusters <- length(unique(Idents(seurat_object)))
    
    cat("Current resolution:", current_resolution, " | Number of clusters:", current_clusters, "\n")
  }
  
  df_output <- data.frame(
    slot1.x = seurat_object@meta.data[, c("col")],  
    slot1.y = seurat_object@meta.data[, c("row")],  
    spatial_cluster_leiden = Idents(seurat_object)
  )
  rownames(df_output) <- paste(df_output$slot1.x, df_output$slot1.y, sep = "x")
  
  return(df_output)
}

run_BayesSpace <- function(sce_object, R, platform){
  
  sce_object <- spatialPreprocess(
    sce_object,
    platform = platform,
    n.PCs = 7,
    n.HVGs = 150,
    skip.PCA = FALSE,
    log.normalize = TRUE,
    assay.type = "logcounts"
  )
  
  single_sim_SRT <- BayesSpace::spatialCluster(sce=sce_object, 
                                               q=R, 
                                               d=7, 
                                               platform=platform, 
                                               nrep=50000, 
                                               gamma=2)
  
  df_output <- data.frame(
    slot1.x = colData(single_sim_SRT)$col,  
    slot1.y = colData(single_sim_SRT)$row,  
    spatial_cluster_BayesSpace = colData(single_sim_SRT)$spatial.cluster
  )
  
  rownames(df_output) <- paste(df_output$slot1.x, df_output$slot1.y, sep = "x")
  
  
  return(df_output)
}

library(mclust)
library(aricode)

dataset_filenames <- c(
  "real.rds"
)

dataset_filenames <- c(
  "scDesign2.rds", 
  "SPARsim.rds", "SRTsim.rds", "SRTsim_rf.rds", "scDesign3_poi_rf.rds", "scDesign3_poi.rds", "scDesign3_nb_rf.rds", "scDesign3_nb.rds",  "scDesign3_gau_rf.rds", "scDesign3_gau.rds", "zinbwave.rds", "symsim.rds"
)


base <- "Other_SRT" 
platform <- "ST"
real <- readRDS("~/spatial_simulationV2/data/final/EH8230.rds")
R <- length(unique(real$manual_anno))

base_dir <- "~/spatial_simulationV2/spatialTask/clustering"

df_ARI <- data.frame(Dataset = character(), BASS = numeric(), SpatialPCA = numeric(),
                     PRECAST = numeric(), DESC = numeric(), leiden = numeric(), 
                     BayesSpace = numeric(), stringsAsFactors = FALSE)

df_NMI <- data.frame(Dataset = character(), BASS = numeric(), SpatialPCA = numeric(),
                     PRECAST = numeric(), DESC = numeric(), leiden = numeric(), 
                     BayesSpace = numeric(), stringsAsFactors = FALSE)

for (filename in dataset_filenames) {
  
  set.seed(123)
  dataset_name <- sub("\\.rds$", "", filename)  # Remove .rds extension to get dataset name
  
  # Define the file path
  if (filename == "real.rds") {
    file_path <- "~/spatial_simulationV2/data/final/EH8230.rds"
  } else {
    file_path <- file.path("~/spatial_simulationV2", "output", "overall", dataset_name, "EH8230/final_sim.rds")
  }
  
  real <- readRDS("~/spatial_simulationV2/data/final/EH8230.rds")
  R <- length(unique(real$manual_anno))
  
  # Load the dataset
  dataset <- readRDS(file_path)
  selected_columns <- sample(ncol(dataset), 800)
  dataset <- dataset[, selected_columns]
  real <- real[, selected_columns]
  
  colnames(dataset) <- gsub("_", "-", colnames(dataset))
  assay(dataset, "logcounts") <- log1p(assay(dataset, "counts"))
  
  colnames(real) <- gsub("_", "-", colnames(real))
  assay(real, "logcounts") <- log1p(assay(real, "counts"))
  
  ari_BASS <- 0
  ari_SpatialPCA <- 0
  ari_PRECAST <- 0
  ari_DESC <- 0
  ari_leiden <- 0
  ari_BayesSpace <- 0
  
  nmi_BASS <- 0
  nmi_SpatialPCA <- 0
  nmi_PRECAST <- 0
  nmi_DESC <- 0
  nmi_leiden <- 0
  nmi_BayesSpace <- 0
  
  tryCatch({
    test_BASS <- run_BASS(dataset, R)
    ari_BASS <- adjustedRandIndex(real$manual_anno, test_BASS$spatial_cluster_BASS)
    nmi_BASS <- NMI(real$manual_anno, test_BASS$spatial_cluster_BASS)
  }, error = function(e) { message("Error in run_BASS for ", dataset_name) })
  
  tryCatch({
    test_SpatialPCA <- run_SpatialPCA(dataset, R)
    ari_SpatialPCA <- adjustedRandIndex(real$manual_anno, test_SpatialPCA$spatial_cluster_SpatialPCA)
    nmi_SpatialPCA <- NMI(real$manual_anno, test_SpatialPCA$spatial_cluster_SpatialPCA)
  }, error = function(e) { message("Error in run_SpatialPCA for ", dataset_name) })
  
  tryCatch({
    test_PRECAST_DESC <- run_PRECAST_DESC(dataset, base, R)
    ari_PRECAST <- adjustedRandIndex(real$manual_anno, test_PRECAST_DESC$spatial_cluster_PRECAST)
    ari_DESC <- adjustedRandIndex(real$manual_anno, test_PRECAST_DESC$spatial_cluster_DESC)
    nmi_PRECAST <- NMI(real$manual_anno, test_PRECAST_DESC$spatial_cluster_PRECAST)
    nmi_DESC <- NMI(real$manual_anno, test_PRECAST_DESC$spatial_cluster_DESC)
  }, error = function(e) { message("Error in run_PRECAST_DESC for ", dataset_name) })
  
  tryCatch({
    test_leiden <- run_leiden(dataset, R)
    ari_leiden <- adjustedRandIndex(real$manual_anno, test_leiden$spatial_cluster_leiden)
    nmi_leiden <- NMI(real$manual_anno, test_leiden$spatial_cluster_leiden)
  }, error = function(e) { message("Error in run_leiden for ", dataset_name) })
  
  tryCatch({
    test_BayesSpace <- run_BayesSpace(dataset, R, platform)
    ari_BayesSpace <- adjustedRandIndex(real$manual_anno, test_BayesSpace$spatial_cluster_BayesSpace)
    nmi_BayesSpace <- NMI(real$manual_anno, test_BayesSpace$spatial_cluster_BayesSpace)
  }, error = function(e) { message("Error in run_BayesSpace for ", dataset_name) })
  
  df_ARI <- rbind(df_ARI, data.frame(
    Dataset = dataset_name, 
    ARI_leiden = ari_leiden, 
    ARI_BASS = ari_BASS, 
    ARI_SpatialPCA = ari_SpatialPCA, 
    ARI_PRECAST = ari_PRECAST, 
    ARI_DE.SC = ari_DESC, 
    ARI_BayesSpace = ari_BayesSpace
  ))
  
  df_NMI <- rbind(df_NMI, data.frame(
    Dataset = dataset_name, 
    NMI_leiden = nmi_leiden,
    NMI_BASS = nmi_BASS, 
    NMI_SpatialPCA = nmi_SpatialPCA, 
    NMI_PRECAST = nmi_PRECAST, 
    NMI_DE.SC = nmi_DESC, 
    NMI_BayesSpace = nmi_BayesSpace
  ))
}

print(df_ARI)
print(df_NMI)

write.csv(df_NMI, "NMI_result.csv")
write.csv(df_ARI, "ARI_result.csv")




