eva_feature <- function(real_sce, sim_sce, path){

  message("feature Start")
  fcFracZero <- compute_fraction_zeros_gene(real_sce,sim_sce)
  feature_meanVar <- feature_meanVar_scaled(real_sce, sim_sce)
  genePearson <- corGene_pearson_new(real_sce,sim_sce)
  geneSpearman <- corGene_spearman(real_sce,sim_sce)
  gene_final <- plot_spatial_gene(fcFracZero,feature_meanVar,genePearson,geneSpearman)
  
  write.csv(fcFracZero, file = paste0(path, "/fcFracZero.csv"), row.names = FALSE)
  write.csv(feature_meanVar, file = paste0(path, "/feature_meanVar.csv"), row.names = FALSE)
  write.csv(genePearson, file = paste0(path, "/genePearson.csv"), row.names = FALSE)
  write.csv(geneSpearman, file = paste0(path, "/geneSpearman.csv"), row.names = FALSE)

  
  message("feature End")
  return(gene_final)
}

eva_sample <- function(real_sce, sim_sce, path){
  
  message("sample Start")
  sp_meanVar <- sample_meanVar_scaled(real_sce,sim_sce)
  spTMM <- sampleTMM(real_sce,sim_sce)
  spLibSize <- sampleLibSize(real_sce,sim_sce)
  spfracZero <- compute_fraction_zeros(real_sce,sim_sce)
  ft_meanVar <- feature_meanVar_scaled(real_sce,sim_sce)
  sampleCor_pearson <- corSample_pearson(real_sce, sim_sce)
  spCor_spearman <- corSample_spearman(real_sce, sim_sce)
  spot_final <- plot_spatial_sample(sp_meanVar, sp_distance, spTMM, spLibSize, spfracZero,sampleCor_pearson, spCor_spearman)
  
  write.csv(sp_meanVar, file = paste0(path, "/sp_meanVar.csv"), row.names = FALSE)
  write.csv(spTMM, file = paste0(path, "/spTMM.csv"), row.names = FALSE)
  write.csv(spLibSize, file = paste0(path, "/spLibSize.csv"), row.names = FALSE)
  write.csv(spfracZero, file = paste0(path, "/spfracZero.csv"), row.names = FALSE)
  write.csv(ft_meanVar, file = paste0(path, "/ft_meanVar.csv"), row.names = FALSE)
  write.csv(sampleCor_pearson, file = paste0(path, "/sampleCor_pearson.csv"), row.names = FALSE)
  write.csv(spCor_spearman, file = paste0(path, "/spCor_spearman.csv"), row.names = FALSE)
  
  message("sample End")
  
  return(spot_final)

}


eva_spatialCluster <- function(real_sce, sim_sce, path){
  # sim_sce already has spatial.cluster from bayerspace
  
  # message("Reclassify Start")
  # sim_sce <- reclassify_simsce(real_sce, sim_sce)
  # message("Reclassify End")
  
  cluster_SRT <- data.frame(real = real_sce$spatial.cluster, sim = sim_sce$spatial.cluster)
  df_result <- cluster_metric_BayerSpace(cluster_SRT, data.frame(assay(real_sce, "counts")) , data.frame(assay(sim_sce, "counts")))
  message("spatial cluster plot Start")
  plot <- ggplot(df_result, aes(x=Metric, y=abs(Value))) + 
    geom_bar(stat="identity", fill=c("#313795", "#4575b4", "#74add1", "#abd9e9", "#fdae61", "#f46d43", "#d72f27")) + 
    ylim(0, 1) + 
    labs(title="Downstream: Clustering Metric", y="Value") + 
    theme(text=element_text(size=12 ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.2, fill=NA))
  
  write.csv(df_result, file = paste0(path, "/df_result.csv"), row.names = FALSE)
  
  message("spatial cluster plot End")
  
  return(plot)
  
}

eva_scFeature <- function(real_sce, sim_sce, sampleID, feat_types, sc_type, sc_species, path, predict_prob){
  
  message("scfeature real Start")
  
  real_sce <- scater::logNormCounts(real_sce)
  
  scfeatures_real <- scFeatures::scFeatures(assay(real_sce , "logcounts"),
                                sample = sampleID,
                                # celltype = real_sce@colData$spatial.cluster,
                                spatialCoords = list(colData(real_sce)$row,colData(real_sce)$col),
                                feature_types = feat_types,
                                type = sc_type,
                                species = sc_species,
                                spotProbability = predict_prob) 
  
  message("scfeature real End")
  
  message("scfeature sim Start")
  
  sim_sce <- scater::logNormCounts(sim_sce)
  
  scfeatures_sim <- scFeatures::scFeatures(assay(sim_sce , "logcounts"),
                               sample = sampleID,
                               # celltype = sim_sce@colData$spatial.cluster,
                               spatialCoords = list(colData(real_sce)$row,colData(real_sce)$col),
                               feature_types = feat_types,
                               type = sc_type,
                               species = sc_species,
                               spotProbability = predict_prob) 

  message("scfeature sim End")
  
  saveRDS(scfeatures_real, paste0(path, "/scfeatures_real.rds"))
  saveRDS(scfeatures_sim, paste0(path, "/scfeatures_sim.rds"))
  
  scFeatures_plot <- scFeature_spatial_df_density(scfeatures_real,scfeatures_sim)
  
  return(scFeatures_plot)
  
  
}

sim_bayerspace <- function(real_sce, sim_sce, plat, path){
  
  
  sim_sce <- spatialPreprocess(
    sim_sce,
    platform = plat,
    n.PCs = 7,
    n.HVGs = 150,
    skip.PCA = FALSE,
    log.normalize = TRUE,
    assay.type = "logcounts"
  )
  
  single_sim_SRT <- BayesSpace::spatialCluster(sce=sim_sce, 
                                   q=length(unique(real_sce@colData$spatial.cluster)), 
                                   d=7, 
                                   platform=plat, 
                                   nrep=50000, 
                                   gamma=2)
  saveRDS(single_sim_SRT, paste0(path, "/final_sim.rds"))
  return(single_sim_SRT)
  
}

plot_bayerspace <- function(real_sce, sim_sce){

  palette <- c("#f8766d", "#c77cff", "#03bfc4", "#7cae00","#eca680")

  sim_sce <- reclassify_simsce(real_sce, sim_sce)
  clusterRealSingle <- clusterPlot(real_sce, palette=palette, color="black", size=0.1) + labs(title="Real")
  clusterSimSingle <- clusterPlot(sim_sce, palette=palette, color="black", size=0.1) +labs(title="Sim")
  
  combined_plot <- gridExtra::grid.arrange(clusterRealSingle, clusterSimSingle, ncol = 2)
  
  return(combined_plot)
  
  
}

