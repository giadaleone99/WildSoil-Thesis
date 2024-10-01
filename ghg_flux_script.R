# script to process ghg flux data

library(stringr)
library(dplyr)

# Import data 
flux_data <- read.csv("flux_data/combined_data.csv")

flux_data1 <- flux_data %>%
  mutate(
    base_code = NA,
    treatment = NA,
    NEERE = NA,
    date = NA
  )

# Split the uniqueID into its components
split_values <- str_split_fixed(flux_data$UniqueID, "_", 4)

# Use a loop to fill in the new columns with the split components
for (i in 1:nrow(flux_data)) {
  flux_data1$base_code[i] <- split_values[i, 1]  # First component goes into base_code
  flux_data1$treatment[i] <- split_values[i, 2]  # Second component goes into treatment
  flux_data1$NEERE[i] <- split_values[i, 3]  
  flux_data1$date[i] <- split_values[i, 4] 
}




