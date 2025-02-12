# Load required libraries
library(tidyverse)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)

# Read data
gene_expression_data <- read.csv("Expression_Data.csv")

# Data transformation
long_expression_data <- gene_expression_data %>%
  select(Gene_Symbol, 
         PC5_1, PC5_2, PC5_3, PC5_4,
         PC15_1, PC15_2, PC15_3, PC15_4) %>%
  gather(key = "Condition", value = "Expression", -Gene_Symbol) %>%
  mutate(
    Group = ifelse(grepl("PC5", Condition), "PC5", "PC15"),
    Replicate = as.numeric(gsub("PC[0-9]+_", "", Condition))
  )

# Save output to a PDF file
pdf("Gene_Expression_Analysis.pdf", width = 15, height = 12)

# 1. Create overview plots

# 1.1 Expression patterns of all genes across PC5 and PC15
plot_all_genes <- ggplot(long_expression_data, aes(x = Replicate, y = Expression, 
                                                   color = Gene_Symbol, linetype = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_linetype_manual(values = c("PC5" = "solid", "PC15" = "dashed")) +
  scale_color_viridis_d() +
  theme_bw() +
  labs(
    title = "Expression Profiles of All Genes (PC5 vs PC15)",
    x = "Replicate Number",
    y = "FPKM",
    color = "Gene Symbol",
    linetype = "Group"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# 1.2 Expression patterns of PC5 group
plot_pc5_genes <- long_expression_data %>%
  filter(Group == "PC5") %>%
  ggplot(aes(x = Replicate, y = Expression, color = Gene_Symbol)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  theme_bw() +
  labs(
    title = "Expression Patterns of Genes in PC5 Group",
    x = "Replicate Number",
    y = "FPKM",
    color = "Gene Symbol"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# 1.3 Expression patterns of PC15 group
plot_pc15_genes <- long_expression_data %>%
  filter(Group == "PC15") %>%
  ggplot(aes(x = Replicate, y = Expression, color = Gene_Symbol)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  theme_bw() +
  labs(
    title = "Expression Patterns of Genes",
    x = "Replicate Number",
    y = "FPKM",
    color = "Gene Symbol"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Print the overview plots
print(plot_all_genes)
print(plot_pc5_genes)
print(plot_pc15_genes)

# 2. Create faceted summary plot for all genes
faceted_gene_plot <- ggplot(long_expression_data, aes(x = Replicate, y = Expression, 
                                                      color = Group, group = interaction(Gene_Symbol, Group))) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  scale_color_manual(values = c("PC5" = "red", "PC15" = "blue")) +
  facet_wrap(~Gene_Symbol, scales = "free_y", ncol = 6) +
  theme_bw() +
  labs(
    title = "Gene-specific Expression Profiles Across Groups",
    x = "Replicate Number",
    y = "FPKM",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )

print(faceted_gene_plot)

# 3. Detailed plots for individual genes
genes_per_page <- 12
gene_list <- unique(long_expression_data$Gene_Symbol)
total_genes <- length(gene_list)
total_pages <- ceiling(total_genes / genes_per_page)

for (page in 1:total_pages) {
  # Select genes for the current page
  gene_subset <- gene_list[((page - 1) * genes_per_page + 1):min(page * genes_per_page, total_genes)]
  
  # Create plots for each gene
  plot_list <- lapply(gene_subset, function(gene) {
    gene_data <- long_expression_data %>% filter(Gene_Symbol == gene)
    gam_models <- list()
    predictions <- list()
    
    for (group in c("PC5", "PC15")) {
      group_data <- gene_data %>% filter(Group == group)
      gam_models[[group]] <- gam(Expression ~ s(Replicate, k = 4), data = group_data)
      pred_data <- data.frame(Replicate = seq(1, 4, length.out = 100))
      pred_data$Fitted <- predict(gam_models[[group]], newdata = pred_data)
      pred_data$Group <- group
      predictions[[group]] <- pred_data
    }
    
    combined_predictions <- bind_rows(predictions)
    
    ggplot() +
      geom_line(data = combined_predictions, aes(x = Replicate, y = Fitted, color = Group), size = 1) +
      geom_point(data = gene_data, aes(x = Replicate, y = Expression, color = Group), size = 2) +
      scale_color_manual(values = c("PC5" = "red", "PC15" = "blue")) +
      labs(
        title = paste0("Expression of ", gene),
        x = "Replicate Number",
        y = "FPKM"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "top"
      )
  })
  
  # Arrange and print plots for the current page
  grid.arrange(grobs = plot_list, ncol = 3, 
               top = textGrob(paste("Gene Expression Profiles - Page", page),
                              gp = gpar(fontsize = 16, font = 2)))
}

# 4. Statistical summary for all genes
summary_stats <- long_expression_data %>%
  group_by(Gene_Symbol, Group) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SD_Expression = sd(Expression),
    Coefficient_of_Variation = (SD_Expression / Mean_Expression) * 100,
    .groups = "drop"
  ) %>%
  left_join(gene_expression_data %>% select(Gene_Symbol, pval, log2FoldChange), by = "Gene_Symbol")

# Save the statistical summary to a CSV file
write.csv(summary_stats, "gene_expression_statistics.csv", row.names = FALSE)

# Close the PDF device
dev.off()
