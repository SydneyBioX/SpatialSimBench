library(spaSim)
set.seed(2014)
generate_spasim_obj <- function(spatial_sce){
  n_cells <- nrow(colData(spatial_sce))
  color_palette <- c("#f8766d", "#c77cff", "#03bfc4", "#7cae00","#eca680")
  points_matrix <- as.matrix(colData(spatial_sce)[c("col", "row")])
  distance_matrix <- as.matrix(dist(points_matrix))
  diag(distance_matrix) <- Inf
  min_distance <- min(distance_matrix)
  n <- max(colData(spatial_sce)$spatial.cluster)
  idents <- as.character(seq(1, n))
  props <- unname(table(colData(spatial_sce)$spatial.cluster))/sum(unname(table(colData(spatial_sce)$spatial.cluster)))
  
  bg <- simulate_background_cells(n_cells = n_cells,
                                  width = 2000,
                                  height = 2000,
                                  # method = "Hardcore",
                                  min_d = min_distance)
  #oversampling_rate = 1.6,
  #Cell.Type = "Others")
  
  mix_bg <- simulate_mixing(bg_sample = bg,
                            idents = idents,
                            props = props, 
                            plot_image = FALSE,
                            plot_colours = color_palette[1:n])
  
  return(mix_bg)
  
}

datasets <- c("BREAST", "HOSTEOSARCOMA", "HPROSTATE", "MBRAIN", "MCATUMOR", "MCORTEX", "MGASTRULA", "MOBNEW", "MUSCLE", "PDAC")

# Loop over each dataset
for (dataset in datasets) {
  # Load the dataset
  data_path <- paste0("~/spatial_simulationV2/data/final/", dataset, ".rds")
  data_obj <- readRDS(data_path)
  
  # Generate simulation object
  simulated_location <- generate_spasim_obj(data_obj)
  
  # Create output directory if it doesn't exist
  output_dir <- paste0("output/others/spaSim/", dataset)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the simulated location object
  output_path <- paste0(output_dir, "/final_sim_new.rds")
  saveRDS(simulated_location, output_path)
}

# Print message when the process is complete
cat("All datasets have been processed and saved.\n")