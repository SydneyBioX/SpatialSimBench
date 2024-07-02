
# panel a

## UMAP

obj <- readRDS("data/final/MOBNEW.rds")
obj <- runUMAP(obj)
umap <- as.data.frame(reducedDim(obj, "UMAP"))
umap$cluster <- as.factor(obj$spatial.cluster)
cluster_colors <- c('1' = '#f8766d', '2' = '#c77cff', '3' ='#03bfc4', '4' = "#7cae00")
UMAP <- ggplot(umap, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cluster_colors) + th
UMAP

ggplot2::ggsave("final_result/fig2/UMAP.pdf", UMAP, width = 8, height = 6, dpi=300)

## heatmap

design_matrix <- model.matrix(~ 0 + as.factor(obj$spatial.cluster))
colnames(design_matrix) <- c("cluster1", "cluster2", "cluster3", "cluster4")
contrast_matrix <- makeContrasts(cluster1VSall = cluster1 - (cluster2 + cluster3 + cluster4)/3,
                                 cluster2VSall = cluster2 - (cluster1 + cluster3 + cluster4)/3,
                                 cluster3VSall = cluster3 - (cluster1 + cluster2 + cluster4)/3,
                                 cluster4VSall = cluster4 - (cluster1 + cluster2 + cluster3)/3,
                                 levels = colnames(design_matrix))
y <- voom(counts(obj), design = design_matrix, plot = FALSE)
fit <- lmFit(y, design = design_matrix)
fit <- contrasts.fit(fit, contrasts = contrast_matrix)
fit <- eBayes(fit)

tt_cluster1 <- topTable(fit, coef = 1, n = Inf, adjust.method = "BH")
tt_cluster2 <- topTable(fit, coef = 2, n = Inf, adjust.method = "BH")
tt_cluster3 <- topTable(fit, coef = 3, n = Inf, adjust.method = "BH")
tt_cluster4 <- topTable(fit, coef = 4, n = Inf, adjust.method = "BH")

tt_cluster1_filtered <- tt_cluster1[tt_cluster1$adj.P.Val < 0.05, ]
tt_cluster1_ordered <- tt_cluster1_filtered[order(tt_cluster1_filtered$logFC, decreasing = TRUE), ]

tt_cluster2_filtered <- tt_cluster2[tt_cluster2$adj.P.Val < 0.05, ]
tt_cluster2_ordered <- tt_cluster2_filtered[order(tt_cluster2_filtered$logFC, decreasing = TRUE), ]

tt_cluster3_filtered <- tt_cluster3[tt_cluster3$adj.P.Val < 0.05, ]
tt_cluster3_ordered <- tt_cluster3_filtered[order(tt_cluster3_filtered$logFC, decreasing = TRUE), ]

tt_cluster4_filtered <- tt_cluster4[tt_cluster4$adj.P.Val < 0.05, ]
tt_cluster4_ordered <- tt_cluster4_filtered[order(tt_cluster4_filtered$logFC, decreasing = TRUE), ]

cluster1_de_genes <- rownames(tt_cluster1_ordered)[1:10]
cluster2_de_genes <- rownames(tt_cluster2_ordered)[1:10]
cluster3_de_genes <- rownames(tt_cluster3_ordered)[1:10]
cluster4_de_genes <- rownames(tt_cluster4_ordered)[1:10]
de_genes <- c(cluster1_de_genes, cluster2_de_genes, cluster3_de_genes, cluster4_de_genes)
heatmap_counts <- obj[de_genes, ]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

heatmap_counts <- heatmap_counts[, order(colData(heatmap_counts)$spatial.cluster)]
drgs_exprs_mtx <- logcounts(heatmap_counts)
drgs_exprs_mtx <- t(apply(drgs_exprs_mtx, 1, cal_z_score))
df <- data.frame(cluster = as.factor(heatmap_counts$spatial.cluster))
rownames(df) <- make.unique(rownames(heatmap_counts@colData))

cluster_colors <- c('1' = "#f8766d", '2' = "#c77cff", '3' ='#03bfc4', '4' = "#7cae00")
breaks <- seq(-2, 2,length.out = 101)
heatmap <- pheatmap(drgs_exprs_mtx,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    annotation_col = df,
                    show_colnames = FALSE,
                    fontsize_row = 9,
                    fontsize = 5,
                    breaks = breaks,
                    use_raster = TRUE,
                    name = " ",
                    annotation_colors = list(cluster = cluster_colors))

palette <- c("#f8766d", "#c77cff", "#03bfc4","#7cae00")

clusterRealSingle_real <- clusterPlot(obj, palette=palette) +
  labs(title="Real")

# panel b

base_path <- "output/overall/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/fcFracZero.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realCount

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simCount
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

th <-   theme(text=element_text(size=12 ),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=0.2, fill=NA) )

feature_fracZero <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Fraction Zero", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/overall/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/feature_meanVar.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realMean

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simMean
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

feature_mean <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Mean", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/overall/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/feature_meanVar.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realVar

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simVar
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

feature_var <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Variance", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/overall/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/genePearson_new.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$real

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$sim
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

feature_cor <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Pearson Correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/overall/"
file_suffix <- "MOBNEW/feature_meanVar.csv"

models <- c("scDesign2", "SPARsim", "symsim", "splatter", "zinbwave")

merged_data <- data.frame(mean = numeric(), variance = numeric(), type = character())

for (model in models) {
  file_path <- paste0(base_path, model, "/", file_suffix)
  
  data <- read.csv(file_path)
  mean_data <- data$scaled_simMean
  var_data <- data$scaled_simVar
  
  merged_data <- rbind(merged_data, data.frame(mean = mean_data, variance = var_data, Model = model))
}

real_mean <- read.csv(paste0(base_path, models[1], "/", file_suffix))$scaled_realMean
real_var <- read.csv(paste0(base_path, models[1], "/", file_suffix))$scaled_realVar
merged_data <- rbind(merged_data, data.frame(mean = real_mean, variance = real_var, Model = "real"))

custom_colors <- c("real" = "#8dd3c7",      # blue
                   "scDesign2" = "#fffeb3",  # orange
                   "SPARsim" = "#bebada",    # green
                   "symsim" = "#81b1d3",     # red
                   "splatter" = "#fb8073",   # purple
                   "zinbwave" = "#fdb462")   # brown

# Create the scatter plot with custom colors and theme
meanvariance <- ggplot(merged_data, aes(x = mean, y = variance, color = Model)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  labs(title = "Feature Mean vs. Variance", x = "ScaledMean", y = "ScaledVariance") +
  th +
  theme(legend.position = "bottom")

base_path <- "output/overall/"
file_suffixmean <- "MOBNEW/feature_meanVar.csv" 
file_suffizero <- "MOBNEW/fcFracZero.csv"

models <- c("scDesign2", "SPARsim", "symsim", "splatter", "zinbwave")

merged_data <- data.frame(mean = numeric(), variance = numeric(), type = character())

for (model in models) {
  file_pathmean <- paste0(base_path, model, "/", file_suffixmean)
  file_pathzero <- paste0(base_path, model, "/", file_suffizero)
  
  datamean <- read.csv(file_pathmean)
  datazero <- read.csv(file_pathzero)
  mean_data <- datamean$scaled_simMean
  zero_data <- datazero$simCount
  
  merged_data <- rbind(merged_data, data.frame(mean = mean_data, zero = zero_data, Model = model))
}

real_mean <- read.csv(paste0(base_path, models[1], "/", file_suffixmean))$scaled_simMean
real_zero <- read.csv(paste0(base_path, models[1], "/", file_suffizero))$simCount
merged_data <- rbind(merged_data, data.frame(mean = real_mean, zero = real_zero, Model = "real"))


custom_colors <- c("real" = "#8dd3c7",      # blue
                   "scDesign2" = "#fffeb3",  # orange
                   "SPARsim" = "#bebada",    # green
                   "symsim" = "#81b1d3",     # red
                   "splatter" = "#fb8073",   # purple
                   "zinbwave" = "#fdb462")   # brown

# Create the scatter plot with custom colors and theme
meanfraczero <- ggplot(merged_data, aes(x = mean, y = zero, color = Model)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  labs(title = "Feature Mean vs. Fraction Zero", x = "ScaledMean", y = "Fraction zero") +
  th +
  theme(legend.position = "bottom")

# panel c

base_path <- "output/domain/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/spfracZero.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$real

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$sim
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

sample_fracZero <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Fraction Zero", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/domain/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/spLibSize.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$libReal

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$libSim
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

sample_libsize <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Library Size", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/domain/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/sampleCor_pearson_new.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$real

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$sim
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

sample_cor <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Pearson Correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none")

base_path <- "output/domain/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/sp_meanVar.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realVar

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simVar
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

sample_var <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Variance", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/domain/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/sp_meanVar.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realMean

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simMean
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

sample_mean <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample Mean", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

base_path <- "output/domain/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/spTMM.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realTMM

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simTMM
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

sample_TMM <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Sample TMM", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th +
  theme(legend.position = "none")

base_path <- "output/overall/"
models <- c("scDesign2", "SPARsim", "symsim",  "splatter", "zinbwave")
file_suffix <- "MOBNEW/fcFracZero.csv"

read_data <- read.csv(paste0(base_path, models[1], "/", file_suffix))$realCount

merged_data <- data.frame(value = read_data, type = "real")

for (model in models) {
  sim_data <- read.csv(paste0(base_path, model, "/", file_suffix))$simCount
  merged_data <- rbind(merged_data, data.frame(value = sim_data, type = model))
}

feature_fracZero <- ggplot(merged_data, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Feature Fraction Zero", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

## panel d

models <- c("real", "scDesign2", "SPARsim", "symsim", "splatter", "zinbwave")
combined_data <- list() # Initialize an empty list to store data frames

for (model in models) {
  file_path <- sprintf("output/domain/%s/MOBNEW/scfeatures_sim.rds", model) # Change path as needed
  scfeatures_sim <- readRDS(file_path)
  morans_I <- data.frame(model = model, value = unname(t(scfeatures_sim$morans_I)[,1]))
  combined_data[[model]] <- morans_I # Add the data frame to the list
}

# Combine all data frames into one
combined_data_df <- bind_rows(combined_data)

moran <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Moran's I", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none")

combined_data <- list() # Initialize an empty list to store data frames

for (model in models) {
  file_path <- sprintf("output/domain/%s/MOBNEW/scfeatures_sim.rds", model) # Change path as needed
  scfeatures_sim <- readRDS(file_path)
  nn_correlation <- data.frame(model = model, value = unname(t(scfeatures_sim$nn_correlation)[,1]))
  combined_data[[model]] <- nn_correlation # Add the data frame to the list
}

# Combine all data frames into one
combined_data_df <- bind_rows(combined_data)


nn_correlation <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "nn_correlation", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none") 

for (model in models) {
  file_path <- sprintf("output/domain/%s/MOBNEW/scfeatures_sim.rds", model) # Change path as needed
  scfeatures_sim <- readRDS(file_path)
  L_stats <- data.frame(model = model, value = unname(t(scfeatures_sim$L_stats)[,1]))
  combined_data[[model]] <- L_stats # Add the data frame to the list
}

# Combine all data frames into one
combined_data_df <- bind_rows(combined_data)


L_stats <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "L_stats", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none") 

combined_data <- list() # Initialize an empty list to store data frames

for (model in models) {
  file_path <- sprintf("output/domain/%s/MOBNEW/scfeatures_sim.rds", model) # Change path as needed
  scfeatures_sim <- readRDS(file_path)
  celltype_interaction <- data.frame(model = model, value = unname(t(scfeatures_sim$celltype_interaction)[,1]))
  combined_data[[model]] <- celltype_interaction # Add the data frame to the list
}

# Combine all data frames into one
combined_data_df <- bind_rows(combined_data)


celltype_interaction <- ggplot(combined_data_df, aes(x = model, y = value, fill = model)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "celltype_interaction", x = "Model", y = "Value") +
  scale_fill_brewer(palette = "Set3") + th  +
  theme(legend.position = "none") 

