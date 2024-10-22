library(Seurat)
library(limma)
library(scater)
library(pheatmap)
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

library(patchwork)
library(readr)
library(ggforce) 
library(tibble)


# panel a

scDesign2_sce <- readRDS("output/overall/scDesign2/MOBNEW/final_sim.rds")
scDesign2_sce <- reclassify_simsce(obj, scDesign2_sce)
scDesign3_gau_sce <- readRDS("output/overall/scDesign3_gau/MOBNEW/final_sim.rds")
scDesign3_gau_sce <- reclassify_simsce(obj, scDesign3_gau_sce)
scDesign3_nb_sce <- readRDS("output/overall/scDesign3_nb/MOBNEW/final_sim.rds")
scDesign3_nb_sce <- reclassify_simsce(obj, scDesign3_nb_sce)
scDesign3_poi_sce <- readRDS("output/overall/scDesign3_poi/MOBNEW/final_sim.rds")
scDesign3_poi_sce <- reclassify_simsce(obj, scDesign3_poi_sce)
SPAR_sce <- readRDS("output/overall/SPARsim/MOBNEW/final_sim.rds")
SPAR_sce <- reclassify_simsce(obj, SPAR_sce)
symsim_sce <- readRDS("output/overall/symsim/MOBNEW/final_sim.rds")
symsim_sce <- reclassify_simsce(obj, symsim_sce)
zinwave_sce <- readRDS("output/overall/zinbwave/MOBNEW/final_sim.rds")
zinwave_sce <- reclassify_simsce(obj, zinwave_sce)
splatter_sce <- readRDS("output/overall/splatter/MOBNEW/final_sim.rds")
SRTsim_sce <- readRDS("output/overall/SRTsim/MOBNEW/final_sim_new.rds")
SRTsim_sce <- reclassify_simsce(obj, SRTsim_sce)
SRTsim_sce@colData

remap_cluster <- function(cluster) {
  if (cluster == 1) return(3)
  else if (cluster == 2) return(1)
  else if (cluster == 3) return(2)
  else return(cluster)  # Keep other values unchanged
}

# Update the "spatial.cluster" metadata column
SRTsim_sce@colData$spatial.cluster <- sapply(SRTsim_sce@colData$spatial.cluster, remap_cluster)


clusterRealSingle_scDesign2 <- clusterPlot(scDesign2_sce, palette=palette) + labs(title="scDesign2")
clusterRealSingle_scDesign3_gau <- clusterPlot(scDesign3_gau_sce, palette=palette) +
  labs(title="scDesign3_gau")
clusterRealSingle_scDesign3_nb <- clusterPlot(scDesign3_nb_sce, palette=palette) +
  labs(title="scDesign3_nb")
clusterRealSingle_scDesign3_poi <- clusterPlot(scDesign3_poi_sce, palette=palette) +
  labs(title="scDesign3_poi")
clusterRealSingle_SPAR <- clusterPlot(SPAR_sce, palette=palette) +
  labs(title="SPARsim")
clusterRealSingle_symsim <- clusterPlot(symsim_sce, palette=palette) +
  labs(title="symsim")
clusterRealSingle_zinwave <- clusterPlot(zinwave_sce, palette=palette) +
  labs(title="zinwave")
clusterRealSingle_SRTsim_sce <- clusterPlot(SRTsim_sce, palette=palette) +
  labs(title="SRTsim")
clusterRealSingle_splatter <- clusterPlot(splatter_sce, palette=palette) +
  labs(title="splatter")


cluster_graph <- clusterRealSingle_real + clusterRealSingle_scDesign2 + clusterRealSingle_SPAR + clusterRealSingle_splatter + clusterRealSingle_symsim + clusterRealSingle_zinwave  

ggplot2::ggsave("final_result/fig2/cluster_graph.pdf", cluster_graph, width = 10, height = 10, dpi=300)

# panel b

models <- c("real","scDesign2", "SPARsim",  "splatter", "symsim","zinbwave")
plots_list <- list()

for(model_name in models) {
  # Construct the file path
  data_path <- sprintf("output/domain/%s/MOBNEW/final_sim_new.rds", model_name)
  
  # Read the RDS file
  sim_data <- readRDS(data_path)
  
  # Prepare the data frame
  celltype_pro <- as.data.frame(sim_data@metadata$celltype_prop)
  celltype_pro <- celltype_pro %>%
    tibble::rownames_to_column(var = "location") %>%
    mutate(
      ax = as.numeric(sub("x.*", "", location)),  
      bx = as.numeric(sub(".*x", "", location))
    )
  
  # Create the plot
  plot <- ggplot(celltype_pro, aes(x=ax, y=bx, group=location)) +
    geom_scatterpie(aes(x=ax, y=bx, group=location), data=celltype_pro, cols=c("GC", "PGC", "M/TC", "OSNs"), color=NA) +
    coord_equal() + 
    labs(x = "X", y = "Y", title = model_name) + theme_void() +  coord_equal() 
  
  # Add the plot to the list
  plots_list[[model_name]] <- plot
}


spot_prob <- ggpubr::ggarrange(plotlist =  plots_list, common.legend = T, ncol = 6, nrow = 1,legend="right")

ggplot2::ggsave("final_result/fig2/spot_prob.pdf", spot_prob, width = 15, height = 7, dpi=300)
