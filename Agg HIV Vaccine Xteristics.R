library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Data definition
data_1e <- data.frame(
  x = c(0, 10, 25, 50, 75, 100, 150, 250, 400, 550, 700, 850, 1000, 1500, 2000),
  Dosage_ann = c(84.765, 78.58697917, 70.03864583, 60.71364583, 51.91572917, 47.36125, 42.58385417, 34.91666667, 26.32385417, 20.74708333, 16.003125, 16.72322917, 11.0809375, 9.920208333, 9.678541667),
  Dosage_one = c(81.06914894, 76.49319149, 71.1143617, 62.2643617, 55.30404255, 47.6762766, 41.44, 33.80244681, 24.27542553, 20.67212766, 16.53744681, 16.4437234, 14.52180851, 12.56904255, 12.30744681),
  Dosage_two = c(82.48894737, 76.77073684, 74.16336842, 65.94378947, 58.05578947, 51.19073684, 44.51294737, 40.57242105, 31.71252632, 25.978, 23.73863158, 22.11315789, 20.73378947, 19.34031579, 19.02705263),
  Mode_muc = c(86.46510638, 78.66, 72.23010638, 62.68053191, 53.52744681, 47.32574468, 40.1256383, 35.27531915, 27.05776596, 22.63106383, 19.48074468, 17.72542553, 16.48180851, 14.12521277, 14.86446809),
  Mode_oral = c(80.17373737, 74.07414141, 67.40808081, 56.78959596, 49.77232323, 43.49434343, 38.46959596, 32.85525253, 23.82919192, 19.03080808, 16.35353535, 15.83323232, 12.37878788, 11.08737374, 9.893030303),
  Mode_sub = c(81.84206522, 79.35380435, 75.98858696, 69.91173913, 62.37771739, 55.8348913, 50.34608696, 41.47021739, 31.7298913, 25.99391304, 20.60652174, 21.93706522, 17.64934783, 16.80152174, 16.48847826)
)

# Define the Koff function with condition-specific k
Koff <- function(x, alpha, Qo, condition, data) {
  if(!is.null(data[[condition]])) {
    k <- log10(max(data$x[!is.na(data[[condition]])])) - log10(min(data$x[data$x > 0 & !is.na(data[[condition]])]))
  } else {
    stop("Invalid condition")
  }
  Qo * 10^(k * (exp(-alpha * Qo * x) - 1))
}

# Function to calculate Pmax
calculate_pmax <- function(alpha, Qo, condition, data) {
  if(!is.null(data[[condition]])) {
    k <- log10(max(data$x[!is.na(data[[condition]])])) - log10(min(data$x[data$x > 0 & !is.na(data[[condition]])]))
  } else {
    stop("Invalid condition")
  }
  if(alpha <= 0 || Qo <= 0 || k <= 0) {
    return(NA)
  }
  pmax <- (1/log(10^k))/(alpha * Qo)
  if(pmax < min(data$x, na.rm = TRUE) || pmax > max(data$x, na.rm = TRUE)) {
    return(NA)
  }
  return(pmax)
}

# Function to fit model and calculate R-squared
fit_and_calculate_rsq <- function(condition, data) {
  valid_data <- data[!is.na(data[[condition]]) & data$x > 0, ]
  
  formula_str <- sprintf("%s ~ Koff(x, alpha, Qo, '%s', data)", condition, condition)
  formula <- as.formula(formula_str)
  
  fit <- tryCatch({
    nls(formula = formula,
        data = valid_data,
        start = list(alpha = 0.0000001, Qo = 100),
        algorithm = "port",
        lower = c(alpha = 0, Qo = 0),
        upper = c(alpha = 0.1, Qo = 100),
        control = nls.control(maxiter = 50000))
  }, error = function(e) {
    message("Error fitting model for ", condition, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    residuals <- residuals(fit)
    tss <- sum((valid_data[[condition]] - mean(valid_data[[condition]], na.rm = TRUE))^2, na.rm = TRUE)
    rss <- sum(residuals^2)
    r_squared <- 1 - (rss / tss)
    
    params <- coef(fit)
    pmax <- calculate_pmax(params["alpha"], params["Qo"], condition, valid_data)
    
    cat("R-squared for", condition, "=", round(r_squared, 4), "\n")
    if(!is.na(pmax)) {
      cat("Pmax for", condition, "=", round(pmax, 2), "\n")
    } else {
      cat("Pmax for", condition, "could not be calculated\n")
    }
    
    return(list(fit = fit, r_squared = r_squared, pmax = pmax))
  }
  
  return(NULL)
}

# Modified create_plot function to handle specific condition groups
create_plot <- function(data, fits, conditions, plot_title, condition_labels) {
  valid_data <- data[data$x > 0, ]
  x_range <- 10^seq(log10(1), log10(max(valid_data$x)), length.out = 100)
  plot_data <- data.frame(x = x_range)
  
  # Define line types and shapes
  linetypes <- c("solid", "dotted", "dashed")
  shapes <- c(16, 17, 15)
  
  # Create the base plot
  p <- ggplot()
  
  # Create a data frame for legend
  legend_data <- data.frame(
    condition = factor(conditions, levels = conditions),
    linetype = factor(linetypes, levels = linetypes),
    shape = shapes
  )
  
  # Add points and lines for each condition
  for(i in seq_along(conditions)) {
    condition <- conditions[i]
    if (!is.null(fits[[condition]]$fit)) {
      plot_data[[paste0(condition, "_pred")]] <- predict(fits[[condition]]$fit, 
                                                         newdata = data.frame(x = x_range))
      
      p <- p + 
        geom_point(data = valid_data, 
                   aes_string(x = "x", y = condition, shape = shQuote(condition_labels[i])), 
                   color = "black") +
        geom_line(data = plot_data, 
                  aes_string(x = "x", y = paste0(condition, "_pred"), 
                             linetype = shQuote(condition_labels[i])),
                  color = "black")
    }
  }
  
  p <- p +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 2000),
                  labels = scales::comma,
                  limits = c(1, 2000)) +
    scale_shape_manual(name = "Condition", values = shapes) +
    scale_linetype_manual(name = "Condition", values = linetypes) +
    labs(title = plot_title,
         x = "Price ($)",
         y = "Average Likelihood of \nAccepting HIV Vaccine (%)") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(face = "bold")) +
    guides(shape = guide_legend(override.aes = list(size = 3)),
           linetype = guide_legend(override.aes = list(size = 1))) +
    annotation_logticks(sides = "b")
  
  return(p)
}

# Define all conditions
conditions <- c("Dosage_ann", "Dosage_one", "Dosage_two", 
                "Mode_muc", "Mode_oral", "Mode_sub")

# Fit models for all conditions
fits <- list()
for(condition in conditions) {
  fits[[condition]] <- fit_and_calculate_rsq(condition, data_1e)
}

# Print parameter values and statistics for each condition
for(condition in conditions) {
  if (!is.null(fits[[condition]]$fit)) {
    cat("\n=================================")
    cat("\nParameters for", condition, ":\n")
    print(summary(fits[[condition]]$fit)$parameters)
    cat("R-squared =", round(fits[[condition]]$r_squared, 4), "\n")
    cat("Pmax =", round(fits[[condition]]$pmax, 2), "\n")
  }
}

# Create two separate plots with custom labels
# Plot 1: Dosage conditions
dosage_conditions <- c("Dosage_ann", "Dosage_one", "Dosage_two")
dosage_labels <- c("Annual", "One-time", "Two-time")
p1 <- create_plot(data_1e, fits, dosage_conditions, 
                  "Aggregate Demand Curves - Dosage Conditions", 
                  dosage_labels)

# Plot 2: Mode conditions
mode_conditions <- c("Mode_muc", "Mode_oral", "Mode_sub")
mode_labels <- c("Mucosal", "Oral", "Sublingual")
p2 <- create_plot(data_1e, fits, mode_conditions, 
                  "Aggregate Demand Curves - Mode of Administration",
                  mode_labels)

# Display both plots
print(p1)
print(p2)