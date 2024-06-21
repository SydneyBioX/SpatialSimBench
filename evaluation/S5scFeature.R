scFeature_spatial_df_density <- function(scfeatures_real, scfeatures_sim){
  
  common_cols <- intersect(colnames(scfeatures_real$L_stats),colnames(scfeatures_sim$L_stats))
  
  scfeatures_real$L_stats <- scfeatures_real$L_stats[, common_cols]
  scfeatures_sim$L_stats <- scfeatures_sim$L_stats[, common_cols]
  
  L_stat <- data.frame(real = unname(t(scfeatures_real$L_stats)[,1]),sim=unname(t(scfeatures_sim$L_stats)[,1]))
  rownames(L_stat) <- colnames(scfeatures_real$L_stats)
  L_stat_long <- L_stat %>%
    pivot_longer(cols = c(real, sim), names_to = "Type", values_to = "Value")
  L_stat_long["metric"] <- "L statistics"
  L_stat_long <- L_stat_long[, c("metric", "Type", "Value")]
  
  celltype_interact <- data.frame(real = unname(t(scfeatures_real$celltype_interaction)[,1]),sim=unname(t(scfeatures_sim$celltype_interaction)[,1]))
  rownames(celltype_interact) <- colnames(scfeatures_real$celltype_interaction)
  celltype_interact_long <- celltype_interact %>%
    pivot_longer(cols = c(real, sim), names_to = "Type", values_to = "Value")
  celltype_interact_long["metric"] <- "Cell type interaction"
  celltype_interact_long <- celltype_interact_long[, c("metric", "Type", "Value")]
  
  nn_correlation <- data.frame(real = unname(t(scfeatures_real$nn_correlation)[,1]),sim=unname(t(scfeatures_sim$nn_correlation)[,1]))
  rownames(nn_correlation) <- colnames(scfeatures_real$nn_correlation)
  nn_correlation_long <- nn_correlation %>%
    pivot_longer(cols = c(real, sim), names_to = "Type", values_to = "Value")
  nn_correlation_long["metric"] <- "Nearest neighbour correlation"
  nn_correlation_long <- nn_correlation_long[, c("metric", "Type", "Value")]
  
  morans_I <- data.frame(real = unname(t(scfeatures_real$morans_I)[,1]),sim=unname(t(scfeatures_sim$morans_I)[,1]))
  rownames(morans_I) <- colnames(scfeatures_real$morans_I)
  morans_I_long <- morans_I %>%
    pivot_longer(cols = c(real, sim), names_to = "Type", values_to = "Value")
  morans_I_long["metric"] <- "Moran's I"
  morans_I_long <- morans_I_long[, c("metric", "Type", "Value")]
  
  final_version <- rbind(L_stat_long,celltype_interact_long,nn_correlation_long,morans_I_long)
  
  # compare_plot <- ggplot(final_version, aes(x=Value, fill=Type)) +
  #  geom_density(alpha=0.5) + 
  #  facet_wrap(~ metric, scales = "free") +
  #  labs(y="Density", x="Value", fill="Type") +
  #  scale_fill_manual(name="Type", values=c(real="#b3202c", sim="#184275")) +
  #  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  compare_plot <- ggplot(final_version, aes(x = Value, fill = Type)) +
    geom_density(alpha = 0.5) + 
    facet_wrap(~ metric, scales = "free") +
    labs(y = "Density", x = "Value", fill = "Type") +
    scale_fill_manual(name = "Type", values = c(real = "#b3202c", sim = "#184275")) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets if needed
      plot.background = element_rect(fill = "white", colour = "white"))
  
  
  
  return(compare_plot)
}