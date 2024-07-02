library(tidyverse)
library(ggplot2)

time_summary <- read_csv("scalability/result/time_summary.csv")
memory_summary <- read_csv("scalability/result/memory_summary.csv")

memory_summary <- memory_summary %>%
  separate(size, into = c("spot_size", "gene_size"), sep = "_") %>%
  mutate(spot_size = as.numeric(spot_size),
         gene_size = as.numeric(gene_size),
         memory_in_gigabytes = as.numeric(gsub("G", "", memory_unit))) # Convert memory_unit to numeric

th <- theme(text=element_text(size=12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=0.2, fill=NA))

time_summary <- time_summary %>%
  separate(size, into = c("spot_size", "gene_size"), sep = "_") %>%
  mutate(spot_size = as.numeric(spot_size),
         gene_size = as.numeric(gene_size))

time_spot_plot_list <- list()
for(gene in c(200, 500, 1000)) {
  
  subset_data <- time_summary %>%
    filter(gene_size == gene)
  
  if(nrow(subset_data) > 0) { 
    plot <- subset_data %>%
      ggplot(aes(x = spot_size, y = mean_time, group = model, color = model)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Gene Size:", gene),
           x = "Spot Size",
           y = "Time (s)",
           color = "Model") +
      th
    
    time_spot_plot_list[[paste("gene_size", gene, sep = "_")]] <- plot
  }
}

memory_spot_plot_list <- list()

for(gene in c(200, 500, 1000)) {
  subset_data <- memory_summary %>%
    filter(gene_size == gene)
  
  if(nrow(subset_data) > 0) { 
    plot <- subset_data %>%
      ggplot(aes(x = spot_size, y = memory_in_gigabytes, group = model, color = model)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Gene Size:", gene),
           x = "Spot Size",
           y = "Memory Usage (G)",
           color = "Model") +
      th
    
    memory_spot_plot_list[[paste("gene_size", gene, sep = "_")]] <- plot
  }
}

memory_gene_plot_list <- list()

for(spot in c(200, 500, 1000, 3000, 5000)) {
  subset_data <- memory_summary %>%
    filter(spot_size == spot)
  
  if(nrow(subset_data) > 0) { 
    plot <- subset_data %>%
      ggplot(aes(x = gene_size, y = memory_in_gigabytes, group = model, color = model)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Spot Size:", spot),
           x = "Gene Size",
           y = "Memory Usage (G)",
           color = "Model") +
      th +
      scale_x_continuous(breaks = c(200, 500, 1000)) 
    
    memory_gene_plot_list[[paste("spot_size", spot, sep = "_")]] <- plot
  }
}

time_gene_plot_list <- list()

for(spot in c(200, 500, 1000, 3000, 5000)) {
  subset_data <- time_summary %>%
    filter(spot_size == spot)
  
  if(nrow(subset_data) > 0) { 
    plot <- subset_data %>%
      ggplot(aes(x = gene_size, y = mean_time, group = model, color = model)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Spot Size:", spot),
           x = "Gene Size",
           y = "Time (s)",
           color = "Model") +
      th +
      scale_x_continuous(breaks = c(200, 500, 1000)) 
    
    time_gene_plot_list[[paste("spot_size", spot, sep = "_")]] <- plot
  }
}

time_spot_plot_list
memory_spot_plot_list
time_gene_plot_list
memory_gene_plot_list