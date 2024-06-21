sample_meanVar_scaled <- function(sce_real, sce_sim){
  real_matrix_sample <- t(counts(sce_real))
  sim_matrix_sample <- t(counts(sce_sim))
  
  # Calculate variance and mean for each sample
  simVar <- apply(sim_matrix_sample, 1, var)
  realVar <- apply(real_matrix_sample, 1, var)
  
  simMean <- apply(sim_matrix_sample, 1, mean)
  realMean <- apply(real_matrix_sample, 1, mean)
  
  # Z-score standardize the variance and mean for each sample
  scaled_simVar <- (simVar - mean(simVar)) / sd(simVar)
  scaled_realVar <- (realVar - mean(realVar)) / sd(realVar)
  
  scaled_simMean <- (simMean - mean(simMean)) / sd(simMean)
  scaled_realMean <- (realMean - mean(realMean)) / sd(realMean)
  
  sample_var_spatial <- data.frame(simVar = simVar,
                                   scaled_simVar = scaled_simVar,
                                   realVar = realVar,
                                   scaled_realVar = scaled_realVar)
  
  sample_mean_spatial <- data.frame(simMean = simMean,
                                    scaled_simMean = scaled_simMean,
                                    realMean = realMean,
                                    scaled_realMean = scaled_realMean)
  
  sample_sp <- cbind(sample_var_spatial, sample_mean_spatial)
  
  return(sample_sp)
}



sampleTMM <- function(sce_real, sce_sim){
  
  real_counts <- assay(sce_real, "counts")
  sim_counts <- assay(sce_sim, "counts")
  real_dge <- edgeR::DGEList(counts = real_counts)
  sim_dge <- edgeR::DGEList(counts = sim_counts)
  real_dge <- edgeR::calcNormFactors(real_dge, method = "TMM")
  sim_dge <- edgeR::calcNormFactors(sim_dge, method = "TMM")
  real_TMM_factors <- real_dge$samples$norm.factors
  sim_TMM_factors <- sim_dge$samples$norm.factors
  TMMfactor <- cbind(data.frame(realTMM = real_TMM_factors), data.frame(simTMM = sim_TMM_factors))
  
  return(TMMfactor)
  
}



sampleLibSize <- function(sce_real, sce_sim){
  
  libReal <- log1p(rowSums(t(counts(sce_real))))
  libSim <- log1p(rowSums(t(counts(sce_sim))))
  libSize <- cbind(data.frame(libReal), data.frame(libSim))
  
  return(libSize)
  
}

compute_fraction_zeros <- function(sce_real, sce_sim) {
  
  fraction_zeros <- function(data) {
    sapply(data, function(col) sum(col == 0) / length(col))
  }
  
  real_counts <- data.frame(assay(sce_real, "counts"))
  colnames(real_counts) <- sub("^X", "", colnames(real_counts))
  
  sim_counts <- data.frame(assay(sce_sim, "counts"))
  colnames(sim_counts) <- sub("^X", "", colnames(sim_counts))
  
  
  
  countDf <- cbind(realCount = data.frame(fraction_zeros(real_counts)), simCount = data.frame(fraction_zeros(sim_counts)))
  names(countDf) <- c("realCount", "simCount")
  
  return(countDf)
  
}

corSample_pearson <- function(real_sce, sim_sce) {
  expression_matrix_real <- as.matrix(assay(real_sce,"counts"))
  correlation_long_real <- reshape2::melt(cor(expression_matrix_real, method = "pearson"))
  correlation_long_real$Dataset <- "real"
  
  expression_matrix_sim <- as.matrix(assay(sim_sce, "counts"))
  correlation_long_sim <- reshape2::melt(cor(expression_matrix_sim, method = "pearson"))
  correlation_long_sim$Dataset <- "sim"
  
  combined_data_pearson  <- rbind(correlation_long_real, correlation_long_sim)
  
  
  return(combined_data_pearson)
}

corSample_spearman <- function(real_sce, sim_sce) {
  expression_matrix_real <- as.matrix(assay(real_sce,"counts"))
  correlation_long_real <- reshape2::melt(cor(expression_matrix_real, method = "spearman"))
  correlation_long_real$Dataset <- "real"
  
  expression_matrix_sim <- as.matrix(assay(sim_sce, "counts"))
  correlation_long_sim <- reshape2::melt(cor(expression_matrix_sim, method = "spearman"))
  correlation_long_sim$Dataset <- "sim"
  
  combined_data_spearman  <- rbind(correlation_long_real, correlation_long_sim)
  
  
  return(combined_data_spearman)
}

plot_spatial_sample <- function(sample_meanVar, sample_distance, TMMfactor, libSize,fracZero,sampleCor_pearson, spCor_spearman){
  th <-   theme(text=element_text(size=12 ),
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(colour = "black", size=0.2, fill=NA) )
  
  sampleVar <- ggplot() +
    geom_density(data=sample_meanVar, aes(x=scaled_simVar, fill="simVar"), alpha=0.5) +
    geom_density(data=sample_meanVar, aes(x=scaled_realVar, fill="realVar"), alpha=0.5) +
    labs(title="spotScaleVar", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275" )) +
    scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + th
  
  sampleMean <- ggplot() +
    geom_density(data=sample_meanVar, aes(x=scaled_simMean, fill="simVar"), alpha=0.5) +
    geom_density(data=sample_meanVar, aes(x=scaled_realMean, fill="realVar"), alpha=0.5) +
    labs(title="spotScaleMean", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) +
    scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + th
  
  libSizePlot <- ggplot() +
    geom_density(data=libSize, aes(x=libReal, fill="real"), alpha=0.5) +
    geom_density(data=libSize, aes(x=libSim, fill="sim"), alpha=0.5) +
    labs(title="libSize", x="Value", y="Density") +
    scale_fill_manual(name="data", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="data", values=c(real="#b3202c", sim="#184275")) + th
  
  ymax <- max(
    max(density(TMMfactor$realTMM)$y),
    max(density(TMMfactor$simTMM)$y)
  )
  
  TMMPlot <- ggplot() +
    geom_density(data=TMMfactor, aes(x=realTMM, fill="realTMM"), alpha=0.5) +
    geom_density(data=TMMfactor, aes(x=simTMM, fill="simTMM"), alpha=0.5) +
    labs(title="TMM", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realTMM="#b3202c", simTMM="#184275")) +
    scale_colour_manual(name="Type", values=c(realTMM="#b3202c", simTMM="#184275")) +
    ylim(0, ymax) + th
  
  effLibSimNew <- libSize$libSim * TMMfactor$simTMM
  effLibRealNew <- libSize$libReal * TMMfactor$realTMM
  effLibSize <- cbind(effLibSim = data.frame(effLibSimNew), effLibReal = data.frame(effLibRealNew))
  
  effLib <- ggplot() +
    geom_density(data=effLibSize, aes(x=effLibRealNew, fill="realEffLib"), alpha=0.5) +
    geom_density(data=effLibSize, aes(x=effLibSimNew, fill="simEffLib"), alpha=0.5) +
    labs(title="effLibSize", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realEffLib="#b3202c", simEffLib="#184275")) +
    th  # I replaced 'th' with theme_minimal().
  
  fracZeroPlot <- ggplot() +
    geom_density(data=fracZero, aes(x=simCount, fill="simCount"), alpha=0.5) +
    geom_density(data=fracZero, aes(x=realCount, fill="realCount"), alpha=0.5) +
    labs(title="fracZeroSpot", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realCount="#b3202c", simCount="#184275")) +
    scale_colour_manual(name="Type", values=c(realCount="#b3202c", simCount="#184275")) + th
  
  merged_data <- cbind(libSize, fracZero)
  
  libsize_fraczero_sp <- ggplot(merged_data) +
    geom_point(aes(x=libReal, y=realCount, color="real"), alpha=0.5) +
    geom_point(aes(x=libSim, y=simCount, color="sim"), alpha=0.5) +
    labs(title="libsize_fraczero",x="Library Size", y="Fraction Zero") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + th
  
  density_plot_pearson <- ggplot(sampleCor_pearson, aes(x = value, fill = Dataset)) +
    geom_density(alpha = 0.5) +
    labs(title = "CorSamPear", x = "Pearson Correlation") +
    scale_fill_manual(values = c("real" = "#b3202c", "sim" = "#184275")) + th
  
  density_plot_spearman <- ggplot(spCor_spearman, aes(x = value, fill = Dataset)) +
    geom_density(alpha = 0.5) +
    labs(title = "CorSamSpea", x = "Spearman Correlation") +
    scale_fill_manual(values = c("real" = "#b3202c", "sim" = "#184275")) + th
  
  
  fig_spatial_sample <- list()
  fig_spatial_sample$libSizePlot <- libSizePlot
  fig_spatial_sample$TMMPlot <- TMMPlot
  fig_spatial_sample$effLib <- effLib
  fig_spatial_sample$sampleVar <- sampleVar
  fig_spatial_sample$sampleMean <- sampleMean
  fig_spatial_sample$fracZeroPlot <- fracZeroPlot
  fig_spatial_sample$libsize_fraczero <- libsize_fraczero_sp
  fig_spatial_sample$CorSamPear <- density_plot_pearson
  fig_spatial_sample$CorSamSpea <- density_plot_spearman
  
  
  fig_spatial_sample_final <- ggarrange( plotlist =  fig_spatial_sample ,  common.legend = T)
  return(fig_spatial_sample_final)
}

sampleTMMV1 <- function(sce_real, sce_sim){
  
  real_counts <- assay(sce_real, "counts")
  sim_counts <- assay(sce_sim, "counts")
  real_dge <- edgeR::DGEList(counts = real_counts)
  sim_dge <- edgeR::DGEList(counts = sim_counts)
  threshold <- 1e-10
  real_dge$samples$lib.size[real_dge$samples$lib.size <= 0] <- threshold
  sim_dge$samples$lib.size[sim_dge$samples$lib.size <= 0] <- threshold
  real_dge <- edgeR::calcNormFactors(real_dge, method = "TMM")
  sim_dge <- edgeR::calcNormFactors(sim_dge, method = "TMM")
  real_TMM_factors <- real_dge$samples$norm.factors
  sim_TMM_factors <- sim_dge$samples$norm.factors
  TMMfactor <- cbind(data.frame(realTMM = real_TMM_factors), data.frame(simTMM = sim_TMM_factors))
  
  return(TMMfactor)
  
}