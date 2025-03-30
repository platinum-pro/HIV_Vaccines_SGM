# Load necessary libraries
library(dplyr)

# Read the data
data <- read.csv("~/Downloads/HIV vaccine demand (SGM) - Study 1 data.csv")

# Add dosage and mode labels
data <- data %>%
  mutate(
    dosage = case_when(
      Dosage_ann ~ "Annual",
      Dosage_one ~ "One Dose",
      Dosage_two ~ "Two Doses",
      TRUE ~ NA_character_
    ),
    mode = case_when(
      Mode_muc ~ "Mucosal",
      Mode_oral ~ "Oral",
      Mode_sub ~ "Subcutaneous",
      TRUE ~ NA_character_
    )
  )

# Updated breakpoint calculation: max x where y drops from >0 to 0
breakpoints <- data %>%
  arrange(ResponseId, x) %>%
  group_by(ResponseId) %>%
  mutate(y_lag = lag(y)) %>%
  filter(y == 0 & y_lag > 0) %>%
  summarize(
    breakpoint = max(x, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(breakpoint = ifelse(is.infinite(breakpoint), NA, breakpoint))

# Merge with dosage and mode labels
breakpoints_labeled <- breakpoints %>%
  left_join(data %>% select(ResponseId, dosage, mode) %>% distinct(), by = "ResponseId")

# Final table of breakpoints by participant
breakpoint_table <- breakpoints_labeled %>%
  select(ResponseId, dosage, mode, breakpoint) %>%
  arrange(dosage, mode, breakpoint) %>%
  mutate(breakpoint = round(breakpoint, 4))

# View breakpoint table
print(breakpoint_table)

# --- Mean breakpoint by dosage ---
mean_bp_dosage <- breakpoint_table %>%
  group_by(dosage) %>%
  summarize(
    mean_breakpoint = round(mean(breakpoint, na.rm = TRUE), 4),
    sd_breakpoint = round(sd(breakpoint, na.rm = TRUE), 4),
    n = sum(!is.na(breakpoint)),
    .groups = "drop"
  )

# --- Mean breakpoint by mode ---
mean_bp_mode <- breakpoint_table %>%
  group_by(mode) %>%
  summarize(
    mean_breakpoint = round(mean(breakpoint, na.rm = TRUE), 4),
    sd_breakpoint = round(sd(breakpoint, na.rm = TRUE), 4),
    n = sum(!is.na(breakpoint)),
    .groups = "drop"
  )

# View summary tables
print(mean_bp_dosage)
print(mean_bp_mode)
