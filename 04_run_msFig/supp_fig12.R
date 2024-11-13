library(readr)
library(ggpubr)

sample_meanVar <- read_csv("/Users/cabiria/Downloads/symsim/sp_meanVar.csv")
TMMfactor <- read_csv("/Users/cabiria/Downloads/symsim/spTMM.csv")
libSize <- read_csv("/Users/cabiria/Downloads/symsim/spLibSize.csv")
fracZero <- read_csv("/Users/cabiria/Downloads/symsim/spfracZero.csv")
sampleCor_pearson <- read_csv("/Users/cabiria/Downloads/symsim/sampleCor_pearson_new.csv")

th <-   theme(text=element_text(size=12 ),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=0.2, fill=NA) )

# Add the zstat result as a text annotation in the plot
sampleVar <- ggplot() +
  geom_density(data = sample_meanVar, aes(x = scaled_simVar, fill = "simVar"), alpha = 0.5) +
  geom_density(data = sample_meanVar, aes(x = scaled_realVar, fill = "realVar"), alpha = 0.5) +
  labs(title = "spotScaleVar", x = "Value", y = "Density") +
  scale_fill_manual(name = "Type", values = c(realVar = "#b3202c", simVar = "#184275")) +
  scale_colour_manual(name = "Type", values = c(realVar = "#b3202c", simVar = "#184275")) +
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(sample_meanVar$scaled_realVar), 
                                                                               x2 = as.numeric(sample_meanVar$scaled_simVar))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") + 
  th  

sampleMean <- ggplot() +
  geom_density(data=sample_meanVar, aes(x=scaled_simMean, fill="simVar"), alpha=0.5) +
  geom_density(data=sample_meanVar, aes(x=scaled_realMean, fill="realVar"), alpha=0.5) +
  labs(title="spotScaleMean", x="Value", y="Density") +
  scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) +
  scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + 
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(sample_meanVar$scaled_simMean), 
                                                                               x2 = as.numeric(sample_meanVar$scaled_realMean))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") +th

libSizePlot <- ggplot() +
  geom_density(data=libSize, aes(x=libReal, fill="real"), alpha=0.5) +
  geom_density(data=libSize, aes(x=libSim, fill="sim"), alpha=0.5) +
  labs(title="libSizeSpot", x="Value", y="Density") +
  scale_fill_manual(name="data", values=c(real="#b3202c", sim="#184275")) +
  scale_colour_manual(name="data", values=c(real="#b3202c", sim="#184275")) + 
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(libSize$libReal), 
                                                                               x2 = as.numeric(libSize$libSim))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") + th

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
  ylim(0, ymax) + 
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(TMMfactor$realTMM), 
                                                                               x2 = as.numeric(TMMfactor$simTMM))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  th

effLibSimNew <- libSize$libSim * TMMfactor$simTMM
effLibRealNew <- libSize$libReal * TMMfactor$realTMM
effLibSize <- cbind(effLibSim = data.frame(effLibSimNew), effLibReal = data.frame(effLibRealNew))

effLib <- ggplot() +
  geom_density(data=effLibSize, aes(x=effLibRealNew, fill="realEffLib"), alpha=0.5) +
  geom_density(data=effLibSize, aes(x=effLibSimNew, fill="simEffLib"), alpha=0.5) +
  labs(title="EffLibSizeSpot", x="Value", y="Density") +
  scale_fill_manual(name="Type", values=c(realEffLib="#b3202c", simEffLib="#184275")) +
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(effLibSize$effLibRealNew), 
                                                                               x2 = as.numeric(effLibSize$effLibSimNew))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  th  # I replaced 'th' with theme_minimal().

fracZeroPlot <- ggplot() +
  geom_density(data=fracZero, aes(x=simCount, fill="simCount"), alpha=0.5) +
  geom_density(data=fracZero, aes(x=realCount, fill="realCount"), alpha=0.5) +
  labs(title="fracZeroSpot", x="Value", y="Density") +
  scale_fill_manual(name="Type", values=c(realCount="#b3202c", simCount="#184275")) +
  scale_colour_manual(name="Type", values=c(realCount="#b3202c", simCount="#184275")) + 
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(fracZero$simCount), 
                                                                               x2 = as.numeric(fracZero$realCount))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  th

merged_data <- cbind(libSize, fracZero)

libsize_fraczero_sp <- ggplot(merged_data) +
  geom_point(aes(x=libReal, y=realCount, color="real"), alpha=0.5) +
  geom_point(aes(x=libSim, y=simCount, color="sim"), alpha=0.5) +
  labs(title="libsize_fraczero",x="Library Size", y="Fraction Zero") +
  scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = merged_data$libReal * merged_data$realCount, x2 = merged_data$libSim * merged_data$simCount)$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") + 
  th

density_plot_pearson <- ggplot() +
  geom_density(data=sampleCor_pearson, aes(x=sim, fill="sim"), alpha=0.5) +
  geom_density(data=sampleCor_pearson, aes(x=real, fill="real"), alpha=0.5) +
  labs(title="fracZeroSpot", x="Value", y="Density") +
  scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
  scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
  annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(sampleCor_pearson$sim), 
                                                                               x2 = as.numeric(sampleCor_pearson$real))$zstat, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  th



fig_spatial_sample <- list()
fig_spatial_sample$libSizePlot <- libSizePlot
fig_spatial_sample$TMMPlot <- TMMPlot
fig_spatial_sample$effLib <- effLib
fig_spatial_sample$sampleVar <- sampleVar
fig_spatial_sample$sampleMean <- sampleMean
fig_spatial_sample$fracZeroPlot <- fracZeroPlot
fig_spatial_sample$libsize_fraczero <- libsize_fraczero_sp
fig_spatial_sample$CorSamPear <- density_plot_pearson


fig_spatial_sample_final <- ggarrange( plotlist =  fig_spatial_sample ,  common.legend = T)




ggsave(filename = "SRTsim_fig_sample_final.pdf", 
       plot = fig_spatial_sample_final, 
       width = 10, height = 8)


spot_final <- plot_spatial_sample(sp_meanVar, sp_distance, spTMM, spLibSize, spfracZero,sampleCor_pearson, spCor_spearman)



ggsave(filename = "symsim_sample_final.pdf", 
       plot = fig_spatial_sample_final, 
       width = 10, height = 8)


scfeatures_real <- readRDS("/Users/cabiria/Downloads/symsim/scfeatures_real.rds")
scfeatures_sim <- readRDS("/Users/cabiria/Downloads/symsim/scfeatures_sim.rds")
library(tidyr)

scFeature_spatial_df_density <- function(scfeatures_real, scfeatures_sim){
  
  common_cols <- intersect(colnames(scfeatures_real$L_stats),colnames(scfeatures_sim$L_stats))
  
  scfeatures_real$L_stats <- scfeatures_real$L_stats[, common_cols]
  scfeatures_sim$L_stats <- scfeatures_sim$L_stats[, common_cols]
  
  L_stat <- data.frame(real = unname(t(scfeatures_real$L_stats)[,1]),sim=unname(t(scfeatures_sim$L_stats)[,1]))
  rownames(L_stat) <- colnames(scfeatures_real$L_stats)
  
  celltype_interact <- data.frame(real = unname(t(scfeatures_real$celltype_interaction)[,1]),sim=unname(t(scfeatures_sim$celltype_interaction)[,1]))
  rownames(celltype_interact) <- colnames(scfeatures_real$celltype_interaction)
  
  
  nn_correlation <- data.frame(real = unname(t(scfeatures_real$nn_correlation)[,1]),sim=unname(t(scfeatures_sim$nn_correlation)[,1]))
  rownames(nn_correlation) <- colnames(scfeatures_real$nn_correlation)
  
  
  morans_I <- data.frame(real = unname(t(scfeatures_real$morans_I)[,1]),sim=unname(t(scfeatures_sim$morans_I)[,1]))
  rownames(morans_I) <- colnames(scfeatures_real$morans_I)
  
  morans_I_plot <- ggplot() +
    geom_density(data=morans_I, aes(x=sim, fill="sim"), alpha=0.5) +
    geom_density(data=morans_I, aes(x=real, fill="real"), alpha=0.5) +
    labs(title=" morans_I", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(morans_I$sim), 
                                                                                 x2 = as.numeric(morans_I$real))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  nn_correlation_graph <- ggplot() +
    geom_density(data=nn_correlation, aes(x=sim, fill="sim"), alpha=0.5) +
    geom_density(data=nn_correlation, aes(x=real, fill="real"), alpha=0.5) +
    labs(title="nn_cor", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(nn_correlation$sim), 
                                                                                 x2 = as.numeric(nn_correlation$real))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  celltype_interact_graph <- ggplot() +
    geom_density(data=celltype_interact, aes(x=sim, fill="sim"), alpha=0.5) +
    geom_density(data=celltype_interact, aes(x=real, fill="real"), alpha=0.5) +
    labs(title="ct_interact", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(celltype_interact$sim), 
                                                                                 x2 = as.numeric(celltype_interact$real))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  L_stat_graph <- ggplot() +
    geom_density(data=L_stat, aes(x=sim, fill="sim"), alpha=0.5) +
    geom_density(data=L_stat, aes(x=real, fill="real"), alpha=0.5) +
    labs(title="L_stat", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(L_stat$sim), 
                                                                                 x2 = as.numeric(L_stat$real))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  fig_scFeatures <- list()
  fig_scFeatures$libSizePlot <- morans_I_plot
  fig_scFeatures$TMMPlot <- nn_correlation_graph
  fig_scFeatures$effLib <- celltype_interact_graph
  fig_scFeatures$sampleVar <- L_stat_graph
  
  fig_scFeatures_final <- ggarrange( plotlist =  fig_scFeatures ,  common.legend = T)
  
  
  
  return(fig_scFeatures_final)
}

scFeatures_plot <- scFeature_spatial_df_density(scfeatures_real,scfeatures_sim)

ggsave(filename = "SRTsim_scFeatures_final.pdf", 
       plot = scFeatures_plot, 
       width = 10, height = 8)

fcFracZero <- read_csv("/Users/cabiria/Downloads/SRTsim/fcFracZero.csv")
feature_meanVar <- read_csv("/Users/cabiria/Downloads/SRTsim/feature_meanVar.csv")
combined_data_pearson <- read_csv("/Users/cabiria/Downloads/SRTsim/genePearson_new.csv")


plot_spatial_gene <- function(fcFracZero, feature_meanVar, combined_data_pearson){
  
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
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(fcFracZero$simCount), 
                                                                                 x2 = as.numeric(fcFracZero$realCount))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  featureVar_plot <- ggplot() +
    geom_density(data=feature_meanVar, aes(x=scaled_simVar, fill="simVar"), alpha=0.5) +
    geom_density(data=feature_meanVar, aes(x=scaled_realVar, fill="realVar"), alpha=0.5) +
    labs(title="featureScaleVar", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) +
    scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(feature_meanVar$scaled_simVar), 
                                                                                 x2 = as.numeric(feature_meanVar$scaled_realVar))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +th
  
  featureMean_plot <- ggplot() +
    geom_density(data=feature_meanVar, aes(x=scaled_simMean, fill="simVar"), alpha=0.5) +
    geom_density(data=feature_meanVar, aes(x=scaled_realMean, fill="realVar"), alpha=0.5) +
    labs(title="featureScaleMean", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275" )) +
    scale_colour_manual(name="Type", values=c(realVar="#b3202c", simVar="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(feature_meanVar$scaled_simMean), 
                                                                                 x2 = as.numeric(feature_meanVar$scaled_realMean))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  featureMeanVarScale_plot <- ggplot(feature_meanVar) +
    geom_point(aes(x=scaled_realMean, y=scaled_realVar, color="real"), alpha=0.5) +
    geom_point(aes(x=scaled_simMean, y=scaled_simVar, color="sim"), alpha=0.5) +
    labs(title="mean_variance(scale)",x="mean expression", y="var of gene expression") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = feature_meanVar$scaled_realMean * feature_meanVar$scaled_realVar, x2 = feature_meanVar$scaled_simMean * feature_meanVar$scaled_simVar)$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  featureMeanVar_plot <- ggplot(feature_meanVar) +
    geom_point(aes(x=realMean, y=realVar, color="real"), alpha=0.5) +
    geom_point(aes(x=simMean, y=simVar, color="sim"), alpha=0.5) +
    labs(title="mean_variance",x="mean expression", y="var of gene expression") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = feature_meanVar$realMean * feature_meanVar$realVar, x2 = feature_meanVar$simMean * feature_meanVar$simVar)$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") + th
  
  mean_fraczero <- cbind(fcFracZero, feature_meanVar)
  
  mean_fraczero_plot <-  ggplot(mean_fraczero) +
    geom_point(aes(x=realMean, y=realCount, color="real"), alpha=0.5) +
    geom_point(aes(x=simMean, y=simCount, color="sim"), alpha=0.5) +
    labs(title="mean_fraczero",x="mean expression", y="Fraction Zero") +
    scale_color_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = mean_fraczero$realMean * mean_fraczero$realCount, x2 = mean_fraczero$simMean * mean_fraczero$simCount)$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +th
  
  density_plot_pearson <- ggplot() +
    geom_density(data=combined_data_pearson, aes(x=sim, fill="sim"), alpha=0.5) +
    geom_density(data=combined_data_pearson, aes(x=real, fill="real"), alpha=0.5) +
    labs(title="fracZeroSpot", x="Value", y="Density") +
    scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
    scale_colour_manual(name="Type", values=c(real="#b3202c", sim="#184275")) + 
    annotate("text", x = Inf, y = Inf, label = paste("KS = ", round(ks::kde.test(x1 = as.numeric(combined_data_pearson$sim), 
                                                                                 x2 = as.numeric(combined_data_pearson$real))$zstat, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black") +
    th
  
  fig_spatial_gene <- list()
  fig_spatial_gene$fcFracZero <- fcFracZero_plot
  fig_spatial_gene$featureVar <- featureVar_plot
  fig_spatial_gene$featureMean <- featureMean_plot
  fig_spatial_gene$featureMeanVarScale <- featureMeanVarScale_plot
  fig_spatial_gene$featureMeanVar <- featureMeanVar_plot
  fig_spatial_gene$meanFraczero <- mean_fraczero_plot
  fig_spatial_gene$featurePear <- density_plot_pearson
  
  fig_spatial_gene_final <- ggarrange( plotlist =  fig_spatial_gene ,  common.legend = T)
  
  return(fig_spatial_gene_final)
}

final_spatial_gene <- plot_spatial_gene(fcFracZero, feature_meanVar, combined_data_pearson)

ggsave(filename = "symsim_feature_final.pdf", 
       plot = final_spatial_gene, 
       width = 10, height = 8)