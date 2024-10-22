library(readr)
library(dplyr)
library(tidyr)
library(Hmisc)
library(magrittr)
library(janitor)
library(tidyr)
library(stringr)
library(funkyheatmap)
library(ggthemr)
library(readr)

# process svg

svg_result <- read_csv("~/spatial_simulationV2/spatialTask/SVG/svg_result.csv")
svg_result$...1 <- NULL
svg_result$...2 <- NULL

svg_cleaned <- svg_result %>%
  select(-rank)

svg_wide <- svg_cleaned %>%
  pivot_wider(names_from = model, values_from = SVG)

cindex_results <- list()


for (dataset_name in unique(svg_wide$dataset)) {
  
  df_NMI <- svg_wide %>%
    filter(dataset == dataset_name)
  
  real_values <- as.numeric(df_NMI$real)
  
  tau_results <- numeric()  
  model_names <- character()  
  
  for (model in colnames(df_NMI)[!(colnames(df_NMI) %in% c('dataset', 'real'))]) {
    sim_values <- as.numeric(df_NMI[[model]])
    tau_value <- rcorr.cens(real_values, sim_values)['C Index']
    
    tau_results <- c(tau_results, tau_value)
    model_names <- c(model_names, model)
  }
  
  cindex_results[[dataset_name]] <- setNames(tau_results, model_names)
}

svg_cindex_df <- data.frame(dataset = names(cindex_results))

for (i in seq_along(cindex_results)) {
  dataset_name <- names(cindex_results)[i]
  cindex_values <- cindex_results[[dataset_name]]
  
  for (model_name in names(cindex_values)) {
    svg_cindex_df[i, model_name] <- cindex_values[model_name]
  }
}

model_names <- unique(unlist(lapply(cindex_results, names))) 
dataset_names <- names(cindex_results)  

svg_cindex_df <- data.frame(model = model_names)

for (dataset in dataset_names) {
  svg_cindex_df[[dataset]] <- NA  
  for (model in model_names) {
    if (model %in% names(cindex_results[[dataset]])) {
      svg_cindex_df[svg_cindex_df$model == model, dataset] <- cindex_results[[dataset]][[model]]
    }
  }
}

svg_cindex_df <- svg_cindex_df[2:14, ]

# process spatial clustering

calculate_cindex <- function(ARI_result_path, dataname) {
  
  ARI_result <- read_csv(ARI_result_path)
  ARI_result$...1 <- NULL
  ARI_result <- ARI_result %>%
    filter(Dataset == "real" | Dataset != "real")
  
  real_values <- as.numeric(ARI_result %>%
                              filter(Dataset == "real") %>%
                              select(-Dataset))
  
  cindex_results <- list()
  
  for (dataset_name in ARI_result$Dataset[ARI_result$Dataset != "real"]) {
    sim_values <- as.numeric(ARI_result %>%
                               filter(Dataset == dataset_name) %>%
                               select(-Dataset))
    
    tau_value <- rcorr.cens(real_values, sim_values)['C Index']
    
    cindex_results[[dataset_name]] <- tau_value
  }
  
  cindex_df <- data.frame(Dataset = names(cindex_results), C_Index = unlist(cindex_results))
  colnames(cindex_df)[2] <- dataname
  return(cindex_df)
}

mobnew_cindex <- calculate_cindex("~/spatial_simulationV2/spatialTask/clustering/MOBNEW_ARI_result.csv", "MOBNEW")
mgastrula_cindex <- calculate_cindex("~/spatial_simulationV2/spatialTask/clustering/MGASTRULA_ARI_result.csv", "MGASTRULA")
pdac_cindex <- calculate_cindex("~/spatial_simulationV2/spatialTask/clustering/PDAC_ARI_result.csv", "PDAC")

mobnew_cindex_NMI <- calculate_cindex("~/spatial_simulationV2/spatialTask/clustering/MOBNEW_NMI_result.csv", "MOBNEW")
mgastrula_cindex_NMI <- calculate_cindex("~/spatial_simulationV2/spatialTask/clustering/MGASTRULA_NMI_result.csv", "MGASTRULA")
pdac_cindex_NMI <- calculate_cindex("~/spatial_simulationV2/spatialTask/clustering/PDAC_NMI_result.csv", "PDAC")

ARI_cindex <- mobnew_cindex %>%
  full_join(mgastrula_cindex, by = "Dataset") %>%
  full_join(pdac_cindex, by = "Dataset") %>%
  rename(model = Dataset)

NMI_cindex <- mobnew_cindex_NMI %>%
  full_join(mgastrula_cindex_NMI, by = "Dataset") %>%
  full_join(pdac_cindex_NMI, by = "Dataset") %>%
  rename(model = Dataset)


# heatmap

svg_cindex_df
ARI_cindex
NMI_cindex

calculate_overall_and_rank <- function(data) {
  ranked_data <- data %>%
    mutate(across(-model, ~ rank(-., ties.method = "min")))  
  ranked_data <- ranked_data %>%
    rowwise() %>%
    mutate(average_rank = mean(c_across(-model), na.rm = TRUE))
  
  data <- data %>%
    mutate(overall = ranked_data$average_rank)
  
  return(data)
}


svg_cindex <- calculate_overall_and_rank(svg_cindex_df)
ARI_cindex <- calculate_overall_and_rank(ARI_cindex)
NMI_cindex <- calculate_overall_and_rank(NMI_cindex)
ARI_cindex <- ARI_cindex %>% rename_with(~ paste0(., "_ARI"), -model)
NMI_cindex <- NMI_cindex %>% rename_with(~ paste0(., "_NMI"), -model)

final_cindex <- svg_cindex %>%
  full_join(ARI_cindex, by = "model") %>%
  full_join(NMI_cindex, by = "model") %>%
  mutate(model = factor(model, levels = c(
    "SRTsim", "SRTsim_rf", "SPARsim", "scDesign3_poi", "scDesign3_poi_rf", 
    "scDesign3_nb_rf", "scDesign3_nb", "zinbwave", "splatter", 
    "scDesign3_gau_rf", "scDesign3_gau", "scDesign2", "symsim"
  ))) %>%
  arrange(model) %>%
  rename(method = model)

final_cindex$adaptor <- c("Yes", "No", "Yes", "No", "Yes", "Yes", "No", "Yes", "Yes", "Yes", "No", "Yes", "Yes")
final_cindex$overall <- -final_cindex$overall
final_cindex$overall_ARI <- -final_cindex$overall_ARI
final_cindex$overall_NMI <- -final_cindex$overall_NMI
final_cindex$overall_final <- rowMeans(final_cindex[, c("overall_ARI", "overall_NMI", "overall")], na.rm = TRUE)

column_info_detailed <- tribble(
  ~id,               ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "method",           "method",             NA,                         "text",       NA,          list(hjust = 0, width = 3),
  "adaptor",           "method",             "Adaptor",                         "text",       NA,          list(hjust = 0, width = 1),
  
  "overall_final",  "summary",       "Overall",    "bar",  "palette4",  list(legend=TRUE),
  "overall",  "summary",       "SVG",    "bar",  "palette1",  list(legend=TRUE),
  "overall_ARI",      "summary",       "ARI",    "bar",  "palette2",  list(legend=TRUE),
  "overall_NMI",      "summary",       "NMI",    "bar",  "palette2",  list(legend=TRUE),
  
  "BREAST",     "svg",       "DATA 1",    "funkyrect",  "palette1",  list(legend=TRUE),
  "HOSTEOSARCOMA",      "svg",       "DATA 2",    "funkyrect",  "palette1",  list(legend=TRUE),
  "HPROSTATE",  "svg",       "DATA 3",    "funkyrect",  "palette1",  list(legend=TRUE),
  "MBRAIN",          "svg",       "DATA 4",    "funkyrect",  "palette1",  list(legend=TRUE),
  "MCATUMOR",    "svg",       "DATA 5",    "funkyrect",  "palette1",  list(legend=TRUE),
  "MCORTEX",   "svg",       "DATA 6",    "funkyrect",  "palette1",  list(legend=TRUE),
  "MGASTRULA",  "svg",       "DATA 7",    "funkyrect",  "palette1",  list(legend=TRUE),
  "MOBNEW",    "svg",       "DATA 8",    "funkyrect",  "palette1",  list(legend=TRUE),
  "MUSCLE",     "svg",       "DATA 9",    "funkyrect",  "palette1",  list(legend=TRUE),
  "PDAC",      "svg",       "DATA 10",    "funkyrect",  "palette1",  list(legend=TRUE),
  
  
  "MOBNEW_ARI",    "spatial_domain_ARI",       "DATA 4",    "funkyrect",  "palette2",  list(legend=TRUE),
  "MGASTRULA_ARI",     "spatial_domain_ARI",       "DATA 6",    "funkyrect",  "palette2",  list(legend=TRUE),
  "PDAC_ARI",      "spatial_domain_ARI",       "DATA 10",    "funkyrect",  "palette2",  list(legend=TRUE),
  
  
  "MOBNEW_NMI",    "spatial_domain_NMI",       "DATA 4",    "funkyrect",  "palette2",  list(legend=TRUE),
  "MGASTRULA_NMI",     "spatial_domain_NMI",       "DATA 6",    "funkyrect",  "palette2",  list(legend=TRUE),
  "PDAC_NMI",      "spatial_domain_NMI",       "DATA 10",    "funkyrect",  "palette2",  list(legend=TRUE),
  
  
)


column_groups <- tribble( 
  ~Experiment, ~Category,  ~group,         ~palette,
  "Method", "",  "method",      "palette4", 
  "Summary", "Rank-based overall", "summary",       "palette4",
  
  "SVG", "Proportion of SVG", "svg",       "palette1",
  
  "Cluster", "ARI", "spatial_domain_ARI",       "palette2",
  "Cluster", "NMI", "spatial_domain_NMI",       "palette2",
  
) 


generateColorPalette <- function(baseColor, numColors, rev=FALSE) {
  colors <- colorRampPalette(c(baseColor, "white"))(numColors)
  if (rev) {
    colors=rev(colors)
  }
  return(colors)
}

palettes <- tribble(
  ~palette,             ~colours,
  "palette1",           generateColorPalette("#E39C4B", 100),
  "palette2",           generateColorPalette("#EF6239", 100),
  "palette3",           generateColorPalette("#FFE5B4", 100),
  "palette4",           generateColorPalette("#000000", 100),
)

g1 <- funky_heatmap(
  data = final_cindex,
  column_info = column_info_detailed,
  column_groups = column_groups,
  palettes = palettes,
  # legends = legends_score,
  col_annot_offset = 3.2
)

ggplot2::ggsave("~/spatial_simulationV2/spatialTask/clustering/supp_fig7.pdf", g1, width = 12, height = 15, dpi=300)


