# Load necessary libraries
library(dplyr)

# Load the data
data <- read.csv("~/Downloads/HIV vaccine demand (SGM)- Study 2 data.csv")

# Helper function to check Stein's criteria for a participant
check_stein_violations <- function(df) {
  df <- df %>% arrange(x)
  
  # Criterion A: Trend (requires log10(x) for x > 0)
  df_nonzero <- df %>% filter(x > 0)
  log_price_range <- log10(max(df_nonzero$x)) - log10(min(df_nonzero$x))
  
  model <- lm(y ~ log10(x), data = df_nonzero)
  trend_slope <- coef(model)[["log10(x)"]]
  trend_violation <- trend_slope > -0.025  # Violation if too flat or increasing
  
  # Criterion B: Bounce
  y0 <- df$y[df$x == min(df$x)][1]
  bounce_violation <- any(df$y > y0 * 1.25 & df$x != min(df$x))
  
  # Criterion C: Reversal from zero
  y_seq <- df$y
  zero_runs <- rle(y_seq == 0)
  zero_2plus_idx <- which(zero_runs$values == TRUE & zero_runs$lengths >= 2)
  reversal_violation <- FALSE
  if (length(zero_2plus_idx) > 0) {
    zero_run_end <- cumsum(zero_runs$lengths)[zero_2plus_idx]
    for (end_idx in zero_run_end) {
      if (end_idx < length(y_seq)) {
        post_zero <- y_seq[(end_idx + 1):length(y_seq)]
        if (any(post_zero > 0, na.rm = TRUE)) {
          reversal_violation <- TRUE
          break
        }
      }
    }
  }
  
  return(data.frame(
    trend_violation = trend_violation,
    bounce_violation = bounce_violation,
    reversal_violation = reversal_violation
  ))
}

# Apply the function to each participant (ResponseId)
stein_flags <- data %>%
  group_by(ResponseId) %>%
  group_modify(~ check_stein_violations(.x)) %>%
  ungroup()

# Add ResponseId and compute any_violation flag
stein_flags <- data %>%
  select(ResponseId) %>%
  distinct() %>%
  left_join(stein_flags, by = "ResponseId") %>%
  mutate(any_violation = trend_violation | bounce_violation | reversal_violation)

# View the per-participant flags
print(stein_flags)

# Optional: Save to CSV
write.csv(stein_flags, "~/Downloads/stein_flags_study2.csv", row.names = FALSE)

# Summary of violations
stein_summary <- stein_flags %>%
  summarize(
    n_total = n(),
    n_trend = sum(trend_violation),
    n_bounce = sum(bounce_violation),
    n_reversal = sum(reversal_violation),
    n_any = sum(any_violation)
  )

print(stein_summary)
