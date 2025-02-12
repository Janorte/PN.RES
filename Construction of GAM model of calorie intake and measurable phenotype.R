# Clear workspace
rm(list = ls())

# Load required libraries
library(mgcv)            # For Generalized Additive Models (GAM)
library(RColorBrewer)    # For color palettes
library(geometry)        # For convex hull calculations
library(sp)              # For spatial data manipulation

# Load the dataset
data <- read.csv("Phenotype_data.csv")

# Define traits and their titles for analysis
traits <- c("BodyWeight")  # Trait of interest
titles <- c("Body Weight") # Corresponding title for the trait

# Define nutrient variables and order
nutrient_order <- c("Protein_Calories", "Glucose_Calories", "Fat_Calories")

# Define prediction ranges
x_limits <- c(5, 50)  # Protein Calorie range
y_limits <- c(40, 85) # Glucose Calorie range

# Fix fat calorie value at 10%
fat_calories_fixed <- c(10)

# Axis labels
x_label <- "Protein Calories (%)"
y_label <- "Glucose Calories (%)"

# Generate prediction grid
fit_resolution <- 100
x_new <- seq(min(x_limits, na.rm = TRUE), max(x_limits, na.rm = TRUE), length.out = fit_resolution)
y_new <- seq(min(y_limits, na.rm = TRUE), max(y_limits, na.rm = TRUE), length.out = fit_resolution)
z_new <- fat_calories_fixed

# Create a prediction data frame
prediction_grid <- expand.grid(
  Protein_Calories = x_new,
  Glucose_Calories = y_new,
  Fat_Calories = z_new
)

# Define output files
output_csv <- "BodyWeight_output.csv"
output_pdf <- "BodyWeight_output.pdf"

# Set color palette and levels for visualization
num_colors <- 100
color_palette <- colorRampPalette(c("green", "yellow", "orange"))(num_colors)
num_contour_levels <- 10

# Define function to check if a point is above the convex hull boundary
is_point_above_hull <- function(x, y, hull_points) {
  if (x < min(hull_points[, 1]) || x > max(hull_points[, 1])) return(TRUE)
  
  for (i in 1:(nrow(hull_points) - 1)) {
    if (x >= hull_points[i, 1] && x <= hull_points[i + 1, 1]) {
      slope <- (hull_points[i + 1, 2] - hull_points[i, 2]) / 
        (hull_points[i + 1, 1] - hull_points[i, 1])
      boundary_y <- hull_points[i, 2] + slope * (x - hull_points[i, 1])
      return(y > boundary_y)
    }
  }
  return(TRUE)
}

# Start plotting to PDF
pdf(output_pdf, width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), oma = c(0, 0, 2, 0))

for (k in seq_along(traits)) {
  # Extract trait-specific data
  data$current_outcome <- data[, traits[k]]
  
  # Fit a Generalized Additive Model (GAM)
  gam_model <- gam(current_outcome ~ s(Protein_Calories, Glucose_Calories, k = 5), 
                   data = data, method = "REML")
  
  # Save model summaries to CSV
  write.table(traits[k], file = output_csv, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
  write.table(summary(gam_model)$p.table, file = output_csv, sep = ",", append = TRUE, col.names = NA)
  write.table(summary(gam_model)$s.table, file = output_csv, sep = ",", append = TRUE, col.names = NA)
  write.table(summary(gam_model)$dev.expl * 100, 
              file = output_csv, sep = ",", append = TRUE, row.names = "Dev. Expl.", col.names = FALSE)
  
  # Generate predictions from the GAM model
  predictions <- predict(gam_model, newdata = prediction_grid, type = "response")
  
  # Determine the range of predictions
  prediction_min <- min(predictions, na.rm = TRUE)
  prediction_max <- max(predictions, na.rm = TRUE)
  
  # Extract unique points for convex hull
  points <- unique(data[, c("Protein_Calories", "Glucose_Calories")])
  points <- points[complete.cases(points), ]
  
  if (nrow(points) >= 3) {
    hull_indices <- convhulln(points)
    hull_points <- points[hull_indices, ]
    hull_points <- hull_points[order(hull_points[, 1]), ]
  }
  
  for (z_val in z_new) {
    subset_predictions <- predictions[prediction_grid$Fat_Calories == z_val]
    response_surface <- matrix(subset_predictions, nrow = fit_resolution)
    
    # Apply convex hull mask
    mask <- matrix(FALSE, nrow = fit_resolution, ncol = fit_resolution)
    for (row in 1:fit_resolution) {
      for (col in 1:fit_resolution) {
        x_val <- x_new[col]
        y_val <- y_new[row]
        mask[row, col] <- is_point_above_hull(x_val, y_val, hull_points)
      }
    }
    masked_surface <- response_surface
    masked_surface[mask] <- NA
    
    # Plot response surface
    image(x_new, y_new, masked_surface, col = color_palette, 
          xlab = x_label, ylab = y_label, cex.lab = 2.2)
    contour(x_new, y_new, masked_surface, add = TRUE, 
            levels = pretty(c(prediction_min, prediction_max), num_contour_levels))
    
    # Highlight convex hull
    if (nrow(points) >= 3) {
      polygon(hull_points, border = "blue", lwd = 2, lty = 2)
    }
    
    title(paste(titles[k], "- Response Surface"), line = 2)
  }
}

# Close PDF device
dev.off()