source("S0packages.R")
source("S1setup.R")
source("S2feature.R")
source("S3sample.R")
source("S4clustering.R")
source("S5scFeature.R")
source("S6simulation.R")
source("S7evaluation.R")

# initally set up dic                       
base_path <- "data/final/"
# dataset_names <- c("BREAST.rds", "HOSTEOSARCOMA.rds", "HPROSTATE.rds",
#                    "MBRAIN.rds","MCATUMOR.rds", "MCORTEX.rds", "MGASTRULA.rds",  
#                    "MOBNEW.rds", "MUSCLE.rds", "PDAC.rds")

dataset_names <- c("MOBNEW.rds")

# dataset_names <- c("PDAC.rds","HOSTEOSARCOMA.rds")

ST_list <- c("MOBNEW", "MGASTRULA", "HOSTEOSARCOMA", "MCORTEX", "PDAC")
Human_list <- c("BREAST","PDAC","HOSTEOSARCOMA", "HPROSTATE")

data_dic_list <- paste0(base_path, dataset_names)
dir_path <- "output/domain/symsim"

# active python env
use_condaenv("calibenv", required = TRUE)
source_python("S8spatial.py")

# set up time
filename <- paste0(dir_path, "/time_taken.csv")
headers <- data.frame(dataset=character(), time=numeric())
write.csv(headers, filename, row.names = FALSE, quote = FALSE)


for (data_path in data_dic_list){
  
  data_name <- extract_name(data_path)
  message(paste0(data_name, " start"))
  start_time <- Sys.time()
  
  real_sce <- readRDS(data_path)
  #real_sce <- rebuilt(real_sce_all)
  sim_sce <- sim_symsim(real_sce, "domain")
  
  # Ensure the output directory exists
  store_result_dir <- paste0(dir_path, "/", data_name)
  dir.create(store_result_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (data_name %in% ST_list){
    platform <- "ST"
  } else{
    platform <- "Visium"
  }
  
  if (data_name %in% Human_list){
    sc_species <- "Homo sapiens"
  } else{
    sc_species <- "Mus musculus"
  }
  
  
  # evaluation
  ## feature & sample
  feature <- eva_feature(real_sce, sim_sce, store_result_dir)
  sample <- eva_sample(real_sce, sim_sce, store_result_dir)
  
  ## sample clustering
  sim_sce <- sim_bayerspace(real_sce, sim_sce, platform, store_result_dir)
  bayerspace <- plot_bayerspace(real_sce, sim_sce)
  clustering <- eva_spatialCluster(real_sce, sim_sce, store_result_dir)
  
  ## scFeature
  sampleID <- rep("sample1", ncol(real_sce))
  feat_types <- c("L_stats","celltype_interaction","nn_correlation","morans_I")
  
  sc_type <- "spatial_t"
  
  prob_matrix <- t(real_sce@metadata$celltype_prop)
  
  scFeature <- eva_scFeature(real_sce, sim_sce, sampleID, feat_types, sc_type, sc_species, store_result_dir, prob_matrix)
  
  ## spatial
  df <- data.frame(x = real_sce@colData$row, 
                   y = real_sce@colData$col, 
                   real = real_sce@colData$spatial.cluster,
                   sim = sim_sce@colData$spatial.cluster)
  final_spatial(df, paste0(store_result_dir, "/spatial.png"))
  
  # save
  # Create a list of ggplot objects
  plot_list <- list(feature = feature, 
                    sample = sample, 
                    bayerspace = bayerspace, 
                    clustering = clustering,
                    scFeature = scFeature)
  
  
  
  # Loop through the list to save each plot
  for (plot_name in names(plot_list)) {
    ggsave(filename = paste0(store_result_dir, "/", plot_name, ".png"),
           plot = plot_list[[plot_name]], 
           width = 10, height = 8)
  }
  
  
  # time record
  end_time <- Sys.time()
  time_taken_sec <- as.numeric(end_time - start_time, units = "secs")
  time_taken_df <- data.frame(dataset = data_name, time = time_taken_sec)
  write.table(time_taken_df, filename, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  message(paste0(data_name, " completed"))
  rm(feature, sample)
  gc()
  
}