library(readr)
library(dplyr)

iceh_data <- read_csv(PROJECT_DATA_ICEH, show_col_types = FALSE) %>%
  mutate(
    estimate = estimate / 100,
    standard_error = standard_error / 100
  )
write_csv(iceh_data, "M9_iceh_data.csv")
