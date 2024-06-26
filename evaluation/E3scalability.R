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

generate_plots <- function(data, x_var, y_var, filter_var, filter_values, breaks, y_label) {
  plot_list <- list()
  
  for(value in filter_values) {
    subset_data <- data %>%
      filter(!!sym(filter_var) == value)
    
    if(nrow(subset_data) > 0) {
      plot <- subset_data %>%
        ggplot(aes_string(x = x_var, y = y_var, group = "model", color = "model")) +
        geom_line() +
        geom_point() +
        labs(title = paste(filter_var, ":", value),
             x = x_var,
             y = y_label,
             color = "Model") +
        th +
        scale_x_continuous(breaks = breaks)
      
      plot_list[[paste(filter_var, value, sep = "_")]] <- plot
    }
  }
  
  return(plot_list)
}

# Generate time spot plots
time_spot_plot_list <- generate_plots(time_summary, "spot_size", "mean_time", "gene_size", c(200, 500, 1000), NULL, "Time (s)")

# Generate memory spot plots
memory_spot_plot_list <- generate_plots(memory_summary, "spot_size", "memory_in_gigabytes", "gene_size", c(200, 500, 1000), NULL, "Memory Usage (G)")

# Generate memory gene plots
memory_gene_plot_list <- generate_plots(memory_summary, "gene_size", "memory_in_gigabytes", "spot_size", c(200, 500, 1000, 3000, 5000), c(200, 500, 1000), "Memory Usage (G)")

# Generate time gene plots
time_gene_plot_list <- generate_plots(time_summary, "gene_size", "mean_time", "spot_size", c(200, 500, 1000, 3000, 5000), c(200, 500, 1000), "Time (s)")

empty_plot <- ggplot() + th

time_spot_plot_list <- c(time_spot_plot_list, list(empty_plot))
memory_spot_plot_list <- c(memory_spot_plot_list, list(empty_plot))

spot_plots <- c(time_spot_plot_list, memory_spot_plot_list)
spot_combined_plot <- do.call(ggpubr::ggarrange, c(spot_plots, list(common.legend = TRUE, ncol = 5, nrow = 2, legend = "right")))
gene_plots <- c(time_gene_plot_list, memory_gene_plot_list)
gene_combined_plot <- do.call(ggpubr::ggarrange, c(gene_plots, list(common.legend = TRUE, ncol = 5, nrow = 2, legend = "right")))

plots <- c(time_spot_plot_list, time_gene_plot_list, memory_spot_plot_list,  memory_gene_plot_list)
combined_plot <- do.call(ggpubr::ggarrange, c(plots, list(common.legend = TRUE, ncol = 4, nrow = 4, legend = "right")))

ggplot2::ggsave("scalability/spot_combined_plot.pdf", spot_combined_plot, width = 10, height = 7, dpi=300)
ggplot2::ggsave("scalability/gene_combined_plot.pdf", gene_combined_plot, width = 16, height = 7, dpi=300)
ggplot2::ggsave("scalability/combined_plot.pdf", combined_plot, width = 16, height = 15, dpi=300)