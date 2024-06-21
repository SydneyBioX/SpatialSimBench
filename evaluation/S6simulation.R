sim_SRTsim <- function(real_sce, base){
  set.seed(1)
  real_count <- counts(real_sce)
  real_loc <- data.frame(x = colData(real_sce)$row,y = colData(real_sce)$col, region = colData(real_sce)$spatial.cluster)
  rownames(real_loc) <- rownames(colData(real_sce))
  
  simSRT<- createSRT(count_in=real_count,loc_in =real_loc)
  
  if (base == "domain"){
    message("SRTsim Simulation Start")
    simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")
  }else if (base == "tissue"){
    message("SRTsim Simulation Start")
    simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")
  }else{
    stop("wrong base parameter")
  }
  simSRT1 <- srtsim_count(simSRT1)
  
  counts_single <- as.matrix(simSRT1@simCounts)
  
  col_data <- data.frame(
    row = data.frame(simSRT1@simcolData)$x,
    col = data.frame(simSRT1@simcolData)$y,
    row.names = rownames(data.frame(simSRT1@simcolData))
  )
  
  col_data_sce <- data.frame(row = col_data$row, col = col_data$col)
  rownames(col_data_sce) <- rownames(col_data)
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts_single),
    colData = col_data_sce
  )
  
  message("SRTsim Simulation End")
  return(sce)
}

sim_scDesign3 <- function(real_sce, base, family, usebam_value){

  log_open("/home/xiaoqi/spatial_simulationV2/output/domain/scDesign3_nb/log_file.log")
  # this is for large dataset
  # real_sce <- real_sce[sample(rownames(real_sce), 200), ]
  set.seed(1)
  if (base == "tissue"){
    message("scDesign3 Simulation Start")
    real_sce@colData$cell_type <- "cell_type"
  }else if (base == "domain"){
    message("scDesign3 Simulation Start")
    real_sce@colData$cell_type <- as.character(real_sce@colData$spatial.cluster)
  }else{
    stop("wrong base parameter")
  }
  
  # message("scDesign3 Simulation Start")
  
  sce_simu <- scdesign3(
    sce = real_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = c("row", "col"), # spatial location
    other_covariates = NULL, 
    mu_formula = "s(row, col, bs = 'gp', k= 100)",
    sigma_formula = "1",
    family_use = family, # could be change another distribution
    n_cores = 1,
    usebam = usebam_value,
    corr_formula = "1",
    copula = "gaussian", 
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE,
    parallelization = "mcmapply" 
  )

  simu_sce <- SingleCellExperiment(list(counts = sce_simu$new_count), colData = sce_simu$new_covariate)
  logcounts(simu_sce) <- log1p(counts(simu_sce))
  
  message("scDesign3 Simulation End")
  log_close()
  return(simu_sce)
}

sim_scDesign2 <- function(real_sce, base){
  set.seed(1)
  if (base == "domain"){
    message("scDesign2 Simulation Start")
  }else{
    stop("ONLY domain base")
  }
  
  traincount <- as.matrix(counts(real_sce))
  real_sce$spatial.cluster <- as.character(real_sce$spatial.cluster)
  spatial_cluster_sel <- unique(real_sce$spatial.cluster)
  colnames(traincount) <- real_sce$spatial.cluster
  spatial_cluster_prop <- table(real_sce$spatial.cluster)
  copula_result <- fit_model_scDesign2(traincount, spatial_cluster_sel, sim_method = 'copula', ncores = 8)
  
  sim_count_copula <- simulate_count_scDesign2(copula_result, ncol(traincount),
                                               sim_method = 'copula',
                                               cell_type_prop = prop.table( spatial_cluster_prop))

  
  # simulated_result <- SingleCellExperiment(list(counts = sim_count_copula))
  # simulated_result$spatial.cluster <- colnames(simulated_result)
  
  # ordered_indices_sim <- colData(real_sce)$spatial.cluster
  # sce_sim <- simulated_result[, ordered_indices_sim]
  
  colnames(sim_count_copula) <- colnames(counts(real_sce))
  rownames(sim_count_copula) <- rownames(counts(real_sce))
  sce_sim <- real_sce
  counts(sce_sim) <- sim_count_copula
  
  colData(sce_sim) <- colData(sce_sim)[, c("col", "row")]
  sce_sim@metadata$celltype_prop <- NULL
  sce_sim@metadata$BayesSpace.data <- NULL
  
  message("scDesign2 Simulation End")
  return(sce_sim)
}

sim_splatter <- function(real_sce, base){
  set.seed(1)
  if (base == "domain"){
    message("splatter Simulation Start")
  }else{
    stop("ONLY domain base")
  }
  
  ordered_indices <- order(colData(real_sce)$spatial.cluster)
  real_sce_ordered <- real_sce[, ordered_indices]
  
  simulated_result <- NULL
  for ( thisSpatialCluster in (unique(real_sce_ordered$spatial.cluster)) ){
    print(thisSpatialCluster)
    
    res <- try({ 
      
      sce_thisSpatialCluster <- real_sce_ordered[ , real_sce_ordered$spatial.cluster == thisSpatialCluster]
      params <- splatter::splatEstimate(as.matrix(counts( sce_thisSpatialCluster)))
      sim_thisSpatialCluster <- splatter::splatSimulate(params)
      
      sim_thisSpatialCluster$spatial.cluster <- thisSpatialCluster
      colnames( sim_thisSpatialCluster) <-  paste0(thisSpatialCluster,  colnames(sim_thisSpatialCluster))
      names( rowData(sim_thisSpatialCluster)) <- paste  ( thisSpatialCluster,names(rowData(sim_thisSpatialCluster)))
      
      
      # combine the cell types 
      if (is.null( simulated_result)){
        simulated_result <-  sim_thisSpatialCluster
      }else{
        simulated_result <- SingleCellExperiment::cbind( simulated_result, sim_thisSpatialCluster)
      }
      
    })
    
  }

  
  colnames(simulated_result) <- colnames(real_sce_ordered)
  rownames(simulated_result) <- rownames(real_sce_ordered)
  
  simulated_result_order <- real_sce_ordered
  counts(simulated_result_order) <- counts(simulated_result)
  
  simulated_result_order <- simulated_result_order[,match(colnames(real_sce), colnames(simulated_result_order))]
  simulated_result_order <- simulated_result_order[match(rownames(real_sce), rownames(simulated_result_order)),]
  simulated_result_order@colData$spatial.cluster <- NULL
  message("splatter Simulation End")
  
  return(simulated_result_order)
}

sim_zinbwave <- function(real_sce, base){
  set.seed(1)
  if (base == "domain"){
    message("zinbwave Simulation Start")
  }else{
    stop("ONLY domain base")
  }
  
  ordered_indices <- order(colData(real_sce)$spatial.cluster)
  real_sce_ordered <- real_sce[, ordered_indices]
  # colData(real_sce_ordered)$spatial.cluster <- as.numeric(colData(real_sce_ordered)$spatial.cluster)
  
  multicoreParam <- MulticoreParam(workers = 8)
  
  X = model.matrix(~spatial.cluster, data=colData(real_sce_ordered))
  params <- zinbEstimate(as.matrix(counts(real_sce_ordered)), design.samples=X, BPPARAM = multicoreParam)
  
  simulated_result <- zinbSimulate(params)
  
  
  colnames(simulated_result) <- colnames(real_sce_ordered)
  rownames(simulated_result) <- rownames(real_sce_ordered)
  
  simulated_result_order <- real_sce_ordered
  counts(simulated_result_order) <- counts(simulated_result)
  
  simulated_result_order <- simulated_result_order[,match(colnames(real_sce), colnames(simulated_result_order))]
  simulated_result_order <- simulated_result_order[match(rownames(real_sce), rownames(simulated_result_order)),]
  simulated_result_order@colData$spatial.cluster <- NULL
  message("zinbwave Simulation End")
  
  return(simulated_result_order)
  
}

sim_symsim <- function(real_sce, base){
  set.seed(1)
  if (base == "domain"){
    message("symsim Simulation Start")
  }else{
    stop("ONLY domain base")
  }
  
  simulated_result <- NULL
  tech <-  "UMI"
  
  ordered_indices <- order(colData(real_sce)$spatial.cluster)
  real_sce_ordered <- real_sce[, ordered_indices]
  
  for (thisSpatialCluster in (unique(real_sce_ordered$spatial.cluster)) ){
    
    print(thisSpatialCluster)
    
    res <- try({ 
      # subset to one cell type 
      sce_thiscelltype <- real_sce_ordered[ , real_sce_ordered$spatial.cluster == thisSpatialCluster]
      
      #this is because if some genes are 0 , this will cause error in simulation 
      keep_feature <- rowSums(counts(sce_thiscelltype) > 0) > 0
      sce_thiscelltype_f <- sce_thiscelltype[keep_feature ,]
      
      best_matches_UMI <- BestMatchParams(tech = "UMI",
                                          counts = as.matrix( counts(sce_thiscelltype_f) ) ,
                                          plotfilename = 'best_params.umi.qqplot',
                                          n_optimal=1 ) 
      
      sim_thiscelltype <-  SimulateTrueCounts(ncells_total =  dim(sce_thiscelltype)[2] , 
                                              ngenes =  dim(sce_thiscelltype)[1] , 
                                              evf_type="one.population", 
                                              randseed = 1, 
                                              Sigma =  best_matches_UMI$Sigma[1], 
                                              gene_effects_sd = best_matches_UMI$gene_effects_sd[1],
                                              scale_s = best_matches_UMI$scale_s[1],
                                              gene_effect_prob = best_matches_UMI$gene_effect_prob[1],
                                              prop_hge = best_matches_UMI$prop_hge[1],
                                              mean_hge = best_matches_UMI$mean_hge[1]   )
      
      
      gene_len <- sample(gene_len_pool,  dim(sce_thiscelltype)[1] ,   replace = FALSE)
      
      sim_thiscelltype <- True2ObservedCounts(true_counts = sim_thiscelltype[[1]], 
                                              meta_cell = sim_thiscelltype[[3]], 
                                              protocol = tech, 
                                              alpha_mean = best_matches_UMI$alpha_mean[1], 
                                              alpha_sd = best_matches_UMI$alpha_sd[1] , 
                                              gene_len = gene_len ,
                                              depth_mean = best_matches_UMI$depth_mean[1],
                                              depth_sd = best_matches_UMI$depth_sd[1])
      
      # tidy up the names
      sim_thiscelltype <- SingleCellExperiment( list(counts = sim_thiscelltype$counts ) )
      sim_thiscelltype$spatial.cluster <- thisSpatialCluster
      
      # combine the cell types 
      if (is.null( simulated_result)   ){
        simulated_result <-  sim_thiscelltype
      }else{
        simulated_result  <- SingleCellExperiment::cbind( simulated_result , sim_thiscelltype  )
      }
      
    })
    
  }
  

  
  colnames(simulated_result) <- colnames(real_sce_ordered)
  rownames(simulated_result) <- rownames(real_sce_ordered)
  
  simulated_result_order <- real_sce_ordered
  counts(simulated_result_order) <- counts(simulated_result)
  
  simulated_result_order <- simulated_result_order[,match(colnames(real_sce), colnames(simulated_result_order))]
  simulated_result_order <- simulated_result_order[match(rownames(real_sce), rownames(simulated_result_order)),]
  simulated_result_order@colData$spatial.cluster <- NULL
  message("symsim Simulation End")
  
  return(simulated_result_order)
  
}

sim_SPARsim <- function(real_sce, base){
  set.seed(1)
  ordered_indices <- order(colData(real_sce)$spatial.cluster)
  real_sce_ordered <- real_sce[, ordered_indices]
  
  if (base == "domain"){
    message("SPARsim Simulation Start")
  }else{
    stop("ONLY domain base")
  }
  
  count_matrix <- data.frame(assay(real_sce_ordered)) 
  sce_scran <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(count_matrix)))
  sce_scran <- scran::computeSumFactors(sce_scran, sizes=seq(20, 100, 5), positive=F) 
  
  if (any(sce_scran$sizeFactor <= 0)) {
    threshold <- 1e-10
    sce_scran$sizeFactor[sce_scran$sizeFactor <= 0] <- threshold
  }
  
  count_matrix_norm <- scater::normalizeCounts(sce_scran, log = FALSE) 
  
  count_matrix_conditions <- find_cluster_indices(real_sce_ordered@colData$spatial.cluster)
  
  SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = count_matrix, 
                                                            norm_data = count_matrix_norm, 
                                                            conditions = count_matrix_conditions)
  sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)

  
  colnames(sim_result$count_matrix) <- gsub("\\.", "-", colnames(sim_result$count_matrix))
  
  simulated_result_order <- real_sce_ordered
  assays(simulated_result_order, withDimnames = FALSE) <- list(counts = sim_result$count_matrix)
  
  # counts(simulated_result_order) <- sim_result$count_matrix
  
  simulated_result_order <- simulated_result_order[,match(colnames(real_sce), colnames(simulated_result_order))]
  simulated_result_order <- simulated_result_order[match(rownames(real_sce), rownames(simulated_result_order)),]
  simulated_result_order@colData$spatial.cluster <- NULL
  
  message("SPARsim Simulation End")  
  
  return(simulated_result_order)
  
}