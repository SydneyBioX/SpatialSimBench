feature_meanVar_scaled <- function(sce_real, sce_sim){
  real_matrix_sample <- counts(sce_real)
  sim_matrix_sample <- counts(sce_sim)
  
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

compute_fraction_zeros_gene <- function(sce_real, sce_sim) {
  
  fraction_zeros <- function(data) {
    sapply(data, function(col) sum(col == 0) / length(col))
  }
  
  real_counts <- data.frame(t(assay(sce_real, "counts")))
  sim_counts <- data.frame(t(assay(sce_sim, "counts")))
  countDf <- cbind(realCount = data.frame(fraction_zeros(real_counts)), simCount = data.frame(fraction_zeros(sim_counts)))
  names(countDf) <- c("realCount", "simCount")
  
  return(countDf)
  
}

corGene_pearson <- function(real_sce, sim_sce) {
  
  expression_matrix_real <- as.matrix(t(assay(real_sce,"counts")))
  correlation_long_real <- reshape2::melt(cor(expression_matrix_real, method = "pearson"))
  correlation_long_real$Dataset <- "real"
  
  expression_matrix_sim <- as.matrix(t(assay(sim_sce, "counts")))
  correlation_long_sim <- reshape2::melt(cor(expression_matrix_sim, method = "pearson"))
  correlation_long_sim$Dataset <- "sim"
  
  combined_data_pearson  <- rbind(correlation_long_real, correlation_long_sim)
  
  
  return(combined_data_pearson)
}

corGene_pearson_new <- function(real_sce, sim_sce) {
  # Convert to matrix and compute correlation directly
  expression_matrix_real <- as.matrix(t(assay(real_sce, "counts")))
  correlation_matrix_real <- cor(expression_matrix_real, method = "pearson")
  
  expression_matrix_sim <- as.matrix(t(assay(sim_sce, "counts")))
  correlation_matrix_sim <- cor(expression_matrix_sim, method = "pearson")
  
  # Convert to long format if needed
  correlation_long_real <- as.data.frame(as.table(correlation_matrix_real))
  names(correlation_long_real) <- c("Var1", "Var2", "value")
  correlation_long_real$Dataset <- "real"
  
  correlation_long_sim <- as.data.frame(as.table(correlation_matrix_sim))
  names(correlation_long_sim) <- c("Var1", "Var2", "value")
  correlation_long_sim$Dataset <- "sim"
  
  # Combine and return
  combined_data_pearson <- rbind(correlation_long_real, correlation_long_sim)
  
  return(combined_data_pearson)
}

corGene_spearman <- function(real_sce, sim_sce) {
  expression_matrix_real <- as.matrix(t(assay(real_sce, "counts")))
  correlation_long_real <- reshape2::melt(cor(expression_matrix_real, method = "spearman"))
  correlation_long_real$Dataset <- "real"
  
  expression_matrix_sim <- as.matrix(t(assay(sim_sce, "counts")))
  correlation_long_sim <- reshape2::melt(cor(expression_matrix_sim, method = "spearman"))
  correlation_long_sim$Dataset <- "sim"
  
  combined_data_spearman  <- rbind(correlation_long_real, correlation_long_sim)
  
  
  return(combined_data_spearman)
}

plot_spatial_gene <- function(fcFracZero, feature_meanVar, combined_data_pearson, combined_data_spearman){
  
  th <-   theme(text=element_text(size=12 ),
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(colour = "black", size=0.2, fill=NA) )
  
  fcFracZero_plot <- ggplot() +
    geom_density(data=fcFracZero, aes(x=simCount, fill="sim"), alpha=0.5) +
    geom_density(data=fcFracZero, aes(x=realCount, fill="real"), alpha=0.5) +
    labs(title="fracZeroGene", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + th
  
  featureVar_plot <- ggplot() +
    geom_density(data=feature_meanVar, aes(x=scaled_simVar, fill="simVar"), alpha=0.5) +
    geom_density(data=feature_meanVar, aes(x=scaled_realVar, fill="realVar"), alpha=0.5) +
    labs(title="featureScaleVar", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) +
    scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + th
  
  featureMean_plot <- ggplot() +
    geom_density(data=feature_meanVar, aes(x=scaled_simMean, fill="simVar"), alpha=0.5) +
    geom_density(data=feature_meanVar, aes(x=scaled_realMean, fill="realVar"), alpha=0.5) +
    labs(title="featureScaleMean", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275" )) +
    scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + th
  
  featureMeanVarScale_plot <- ggplot(feature_meanVar) +
    geom_point(aes(x=scaled_realMean, y=scaled_realVar, color="real"), alpha=0.5) +
    geom_point(aes(x=scaled_simMean, y=scaled_simVar, color="sim"), alpha=0.5) +
    labs(title="mean_variance(scale)",x="mean expression", y="var of gene expression") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + th
  
  featureMeanVar_plot <- ggplot(feature_meanVar) +
    geom_point(aes(x=realMean, y=realVar, color="real"), alpha=0.5) +
    geom_point(aes(x=simMean, y=simVar, color="sim"), alpha=0.5) +
    labs(title="mean_variance",x="mean expression", y="var of gene expression") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + th
  
  mean_fraczero <- cbind(fcFracZero, feature_meanVar)
  
  mean_fraczero_plot <-  ggplot(mean_fraczero) +
    geom_point(aes(x=realMean, y=realCount, color="real"), alpha=0.5) +
    geom_point(aes(x=simMean, y=simCount, color="sim"), alpha=0.5) +
    labs(title="mean_fraczero",x="mean expression", y="Fraction Zero") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + th
  
  density_plot_pearson <- ggplot(combined_data_pearson, aes(x = value, fill = Dataset)) +
    geom_density(alpha = 0.5) +
    labs(title = "CorGenePear", x = "Pearson Correlation") +
    scale_fill_manual(values = c("real" = "#b3202c", "sim" = "#184275")) + th
  
  density_plot_spearman <- ggplot(combined_data_spearman, aes(x = value, fill = Dataset)) +
    geom_density(alpha = 0.5) +
    labs(title = "CorGeneSpea", x = "Spearman Correlation") +
    scale_fill_manual(values = c("real" = "#b3202c", "sim" = "#184275")) + th
  
  
  
  fig_spatial_gene <- list()
  fig_spatial_gene$fcFracZero <- fcFracZero_plot
  fig_spatial_gene$featureVar <- featureVar_plot
  fig_spatial_gene$featureMean <- featureMean_plot
  fig_spatial_gene$featureMeanVarScale <- featureMeanVarScale_plot
  fig_spatial_gene$featureMeanVar <- featureMeanVar_plot
  fig_spatial_gene$meanFraczero <- mean_fraczero_plot
  fig_spatial_gene$featurePear <- density_plot_pearson
  fig_spatial_gene$featureSpear <- density_plot_spearman
  
  fig_spatial_gene_final <- ggarrange( plotlist =  fig_spatial_gene ,  common.legend = T)
  
  return(fig_spatial_gene_final)
}