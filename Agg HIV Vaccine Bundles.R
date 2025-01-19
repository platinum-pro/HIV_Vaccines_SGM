library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Data definition
data_1e <- data.frame(
  x = c(0, 10, 25, 50, 75, 100, 150, 250, 400, 550, 700, 850, 1000, 1500, 2000),
  PrEP = c(88.02, 76.43, 71.39, 66.50, 58.03, 49.83, 37.21, 29.92, 25.92, 22.68, 19.70, 17.39, 14.34, 14.41, 15.31),
  HIV_Test = c(91.31, 85.97, 78.79, 67.80, 54.78, 48.49, 41.05, 34.15, 25.23, 19.12, 13.54, 12.05, 8.75, 7.07, 7.62),
  Condoms = c(80.12, 79.83, 71.63, 64.53, 55.93, 48.04, 42.10, 37.18, 29.11, 24.29, 23.22, 18.77, 22.32, 19.44, 17.81),
  COVID_Vacc = c(78.58, 77.99, 75.99, 68.23, 59.18, 50.90, 40.50, 30.36, 22.83, 18.21, 17.98, 16.36, 14.20, 13.84, 12.10),
  Flu_Vacc = c(79.26, 70.50, 68.02, 61.69, 53.13, 51.61, 38.61, 38.88, 30.71, 27.22, 19.77, 16.04, 15.99, 13.23, 13.07),
  HPV_Vacc = c(90.91, 87.03, 81.57, 71.66, 64.41, 58.38, 51.22, 39.32, 31.92, 24.12, 20.46, 17.98, 17.37, 14.10, 13.02),
  Dental = c(81.07, 72.58, 71.10, 63.29, 57.86, 55.34, 46.92, 38.81, 35.07, 30.15, 21.66, 22.54, 18.10, 15.76, 12.06),
  Eye = c(74.49, 74.24, 72.46, 69.65, 66.65, 61.02, 55.38, 41.27, 22.71, 12.70, 6.47, 5.35, 5.35, 4.44, 8.65),
  Counseling = c(81.84, 80.75, 78.76, 75.71, 66.71, 59.24, 53.49, 44.23, 35.70, 29.59, 25.15, 23.29, 20.36, 15.49, 14.16)
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

# Modified create_plot function without colors
create_plot <- function(data, fits, conditions, plot_title) {
  valid_data <- data[data$x > 0, ]
  x_range <- 10^seq(log10(1), log10(max(valid_data$x)), length.out = 100)
  plot_data <- data.frame(x = x_range)
  
  # Define line types and shapes
  linetypes <- c("solid", "dashed", "dotted", "longdash", "twodash", 
                 "dotdash", "solid", "dashed", "dotted")
  shapes <- c(16, 17, 15, 18, 19, 20, 21, 22, 23)
  
  # Create the base plot
  p <- ggplot()
  
  # Add points and lines for each condition
  for(i in seq_along(conditions)) {
    condition <- conditions[i]
    if (!is.null(fits[[condition]]$fit)) {
      plot_data[[paste0(condition, "_pred")]] <- predict(fits[[condition]]$fit, 
                                                         newdata = data.frame(x = x_range))
      
      p <- p + 
        geom_point(data = valid_data, 
                   aes_string(x = "x", y = condition, shape = shQuote(condition)), 
                   color = "black") +
        geom_line(data = plot_data, 
                  aes_string(x = "x", y = paste0(condition, "_pred"), 
                             linetype = shQuote(condition)),
                  color = "black")
    }
  }
  
  p <- p +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 2000),
                  labels = scales::comma,
                  limits = c(1, 2000)) +
    scale_shape_manual(name = "Commodity Bundled \nwith HIV vaccine", values = shapes) +
    scale_linetype_manual(name = "Commodity Bundled \nwith HIV vaccine", values = linetypes) +
    labs(title = plot_title,
         x = "Price ($)",
         y = "Average Likelihood of Accepting \nBundled Commodities (%)") +
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
conditions <- c("PrEP", "HIV_Test", "Condoms", "COVID_Vacc", "Flu_Vacc", 
                "HPV_Vacc", "Dental", "Eye", "Counseling")

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

# Create and display the plot
p <- create_plot(data_1e, fits, conditions, 
                 "Aggregate Demand Curves - Health Services")
print(p)