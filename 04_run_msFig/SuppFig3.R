
# panel a

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/spfracZero.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")

    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$real
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }

    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$sim
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

sample_fracZero_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Fraction Zero", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_fracZero_group 


base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/sp_meanVar.csv"

# Initialize an empty data frame for merged data
merged_data <- data.frame(value = numeric(), type = character(), group = character())

# Loop through each path and model to read and merge data
for (path in base_paths) {
  for (model in models) {
    # Determine group based on the path and rename accordingly
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    
    # Reading real data for the first model to represent the 'real' type only once per group
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$realMean
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    
    # Reading simulated data
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$realMean
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}


dodge <- position_dodge(width=0.9)

sample_mean_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Mean", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_mean_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/sp_meanVar.csv"

# Initialize an empty data frame for merged data
merged_data <- data.frame(value = numeric(), type = character(), group = character())

# Loop through each path and model to read and merge data
for (path in base_paths) {
  for (model in models) {
    # Determine group based on the path and rename accordingly
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    
    # Reading real data for the first model to represent the 'real' type only once per group
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$realVar
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    
    # Reading simulated data
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$realVar
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

sample_var_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Variance", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_var_group 


base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/spLibSize.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$libReal
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$libSim
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

sample_libsize_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Library Size", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_libsize_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/sampleCor_pearson_new.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")

    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$real
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }

    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$sim
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}


dodge <- position_dodge(width=0.9)

sample_pearson_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Pearson Correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_pearson_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/spCor_spearman_new.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {

    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$real
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$sim
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}



dodge <- position_dodge(width=0.9)

sample_spearman_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Spearman Correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_spearman_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/spTMM.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())


for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$realTMM
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$simTMM
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}



dodge <- position_dodge(width=0.9)

sample_tmm_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample TMM", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

sample_tmm_group 

# panel b

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/fcFracZero.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$realCount
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$simCount
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

feature_count_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Fraction Zero", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

feature_count_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/feature_meanVar.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$realMean
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }
    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$simMean
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

feature_mean_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Mean", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

feature_mean_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/feature_meanVar.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$realVar
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }

    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$simVar
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

feature_var_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Variance", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

feature_var_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/genePearson_new.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$real
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }

    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$sim
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

feature_pearson_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Pearson Correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

feature_pearson_group 

base_paths <- c("output/domain/", "output/tissue/")
models <- c("scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
file_suffix <- "MOBNEW/geneSpearman_new.csv"

merged_data <- data.frame(value = numeric(), type = character(), group = character())

for (path in base_paths) {
  for (model in models) {
    group <- ifelse(grepl("domain", path), "Region-Based", "Region-Free")
    if (model == models[1]) {
      real_data <- read.csv(paste0(path, model, "/", file_suffix))$real
      merged_data <- rbind(merged_data, data.frame(value = real_data, type = "real", group = group))
    }

    sim_data <- read.csv(paste0(path, model, "/", file_suffix))$sim
    merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model, group = group))
  }
}

dodge <- position_dodge(width=0.9)

feature_spearman_group <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(aes(group = interaction(type, group)), position = dodge, trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Spearman Correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x")

feature_spearman_group 


# panel c

models <- c("real", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
base_paths <- list(domain = "output/domain/", tissue = "output/tissue/")

combined_data <- list() 

for (base_path_name in names(base_paths)) {
  for (model in models) {
    
    file_path <- sprintf("%s%s/MOBNEW/scfeatures_sim.rds", base_paths[[base_path_name]], model)
    scfeatures_sim <- readRDS(file_path)
    group <- ifelse(base_path_name == "domain", "Region-Based", "Region-Free")
    morans_I <- data.frame(model = model, value = unname(t(scfeatures_sim$morans_I)[,1]), group = group)
    combined_data[[paste(model, base_path_name, sep = "_")]] <- morans_I # Use unique names for list elements
  }
}

combined_data_df <- bind_rows(combined_data)

th <- theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 0.2, fill = NA),
            strip.text.x = element_text(size = 12))

moran_I_group <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Moran's I", x = "", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x", nrow = 1)

moran_I_group

models <- c("real", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
base_paths <- list(domain = "output/domain/", tissue = "output/tissue/")

combined_data <- list() # Initialize an empty list to store data frames

for (base_path_name in names(base_paths)) {
  for (model in models) {
    file_path <- sprintf("%s%s/MOBNEW/scfeatures_sim.rds", base_paths[[base_path_name]], model)
    scfeatures_sim <- readRDS(file_path)
    group <- ifelse(base_path_name == "domain", "Region-Based", "Region-Free")
    nn_correlation <- data.frame(model = model, value = unname(t(scfeatures_sim$nn_correlation)[,1]), group = group)
    combined_data[[paste(model, base_path_name, sep = "_")]] <- nn_correlation # Use unique names for list elements
  }
}

combined_data_df <- bind_rows(combined_data)

nn_correlation_group <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "nn_correlation", x = "", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x", nrow = 1)

nn_correlation_group

models <- c("real", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
base_paths <- list(domain = "output/domain/", tissue = "output/tissue/")

combined_data <- list() 

for (base_path_name in names(base_paths)) {
  for (model in models) {
    file_path <- sprintf("%s%s/MOBNEW/scfeatures_sim.rds", base_paths[[base_path_name]], model)
    scfeatures_sim <- readRDS(file_path)
    group <- ifelse(base_path_name == "domain", "Region-Based", "Region-Free")
    # Prepare the data frame with 'group' column
    L_stats <- data.frame(model = model, value = unname(t(scfeatures_sim$L_stats)[,1]), group = group)
    combined_data[[paste(model, base_path_name, sep = "_")]] <- L_stats # Use unique names for list elements
  }
}

combined_data_df <- bind_rows(combined_data)

L_stats_group <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "L_stats", x = "", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x", nrow = 1)

L_stats_group

models <- c("real", "scDesign3_gau", "scDesign3_nb", "scDesign3_poi", "SRTsim")
base_paths <- list(domain = "output/domain/", tissue = "output/tissue/")

combined_data <- list()

for (base_path_name in names(base_paths)) {
  for (model in models) {
    file_path <- sprintf("%s%s/MOBNEW/scfeatures_sim.rds", base_paths[[base_path_name]], model)
    scfeatures_sim <- readRDS(file_path)
    group <- ifelse(base_path_name == "domain", "Region-Based", "Region-Free")
    # celltype_interaction the data frame with 'group' column
    celltype_interaction <- data.frame(model = model, value = unname(t(scfeatures_sim$celltype_interaction)[,1]), group = group)
    combined_data[[paste(model, base_path_name, sep = "_")]] <- celltype_interaction # Use unique names for list elements
  }
}

combined_data_df <- bind_rows(combined_data)

celltype_interaction_group <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "celltype_interaction", x = "", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none") +
  facet_wrap(~group, scales = "free_x", nrow = 1)

celltype_interaction_group
