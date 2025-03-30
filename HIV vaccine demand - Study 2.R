# Load necessary libraries
library(dplyr)

# Load the dataset
data <- read.csv("~/Downloads/HIV vaccine demand (SGM)- Study 2 data.csv")

# Convert condition indicator columns to logical (TRUE/FALSE)
condition_vars <- c("PrEP", "HIV_Test", "Condoms", "COVID_Vacc", "Flu_Vacc", 
                    "HPV_Vacc", "Dental", "Eye", "Counseling")

data[condition_vars] <- lapply(data[condition_vars], as.logical)

# Identify the assigned bundle condition
data <- data %>%
  mutate(
    condition = case_when(
      PrEP ~ "PrEP",
      HIV_Test ~ "HIV Test",
      Condoms ~ "Condoms",
      COVID_Vacc ~ "COVID Vaccine",
      Flu_Vacc ~ "Flu Vaccine",
      HPV_Vacc ~ "HPV Vaccine",
      Dental ~ "Dental Exam",
      Eye ~ "Eye Exam",
      Counseling ~ "Counseling",
      TRUE ~ NA_character_
    )
  )

# Calculate breakpoint: max x where y drops from >0 to 0
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

# Merge breakpoint info with assigned condition
breakpoints_labeled <- breakpoints %>%
  left_join(data %>% select(ResponseId, condition) %>% distinct(), by = "ResponseId")

# Table of breakpoints by participant
breakpoint_table <- breakpoints_labeled %>%
  select(ResponseId, condition, breakpoint) %>%
  arrange(condition, breakpoint) %>%
  mutate(breakpoint = round(breakpoint, 4))

# View participant-level breakpoint table
print(breakpoint_table)

# --- Mean breakpoint by condition ---
mean_bp_condition <- breakpoint_table %>%
  group_by(condition) %>%
  summarize(
    mean_breakpoint = round(mean(breakpoint, na.rm = TRUE), 4),
    sd_breakpoint = round(sd(breakpoint, na.rm = TRUE), 4),
    n = sum(!is.na(breakpoint)),
    .groups = "drop"
  )

# View summary
print(mean_bp_condition)

# Optional: write to CSV
# write.csv(breakpoint_table, "study2_breakpoints_by_participant.csv", row.names = FALSE)
# write.csv(mean_bp_condition, "study2_mean_breakpoints_by_condition.csv", row.names = FALSE)
