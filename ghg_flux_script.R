# script to process ghg flux data

library(stringr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(lubridate)

# Import data 
flux_data_raw <- read.csv("flux_data/combined_data.csv")

flux_data <- flux_data_raw %>%
  mutate(
    base_code = NA,
    treatment = NA,
    NEERE = NA,
    date = NA
  )

# Split the uniqueID into its components
split_values <- str_split_fixed(flux_data_raw$UniqueID, "_", 4)

# Use a loop to fill in the new columns with the split components
for (i in 1:nrow(flux_data_raw)) {
  flux_data$base_code[i] <- split_values[i, 1]  # First component goes into base_code
  flux_data$treatment[i] <- split_values[i, 2]  # Second component goes into treatment
  flux_data$NEERE[i] <- split_values[i, 3]  
  flux_data$date[i] <- split_values[i, 4] 
}

flux_data <- flux_data %>% relocate(date, NEERE, treatment, base_code, .after = UniqueID)

flux_data <- flux_data %>% 
  mutate(
    Animal = case_when(
      grepl("^C", base_code) ~ "Cow",
      grepl("^H", base_code) ~ "Horse"
    )
  )

flux_data <- flux_data %>% 
  mutate(
    Campaign = case_when(
      grepl("G", base_code) ~ "Gradient",
      grepl("D", base_code) ~ "Daily",
      grepl("^PIT", base_code) ~ "PIT"
    )
  )

flux_data <- flux_data %>% 
  mutate(plotID = paste(base_code, treatment, sep = "_"),
         longdate = dmy(paste0(date, " 2024")))

flux_data <- flux_data %>% 
  mutate(plotNEERE = paste(plotID, NEERE, sep = "_"))

        

# Calculating photosynthesis
flux_data <- flux_data %>% 
  mutate(photosynthesis = NA)
  
for(i in 1:(nrow(flux_data) - 1)) {  # Loop until n-1 to avoid going out of bounds
  if (flux_data$NEERE[i] == "NEE" && flux_data$NEERE[i+1] == "RE") {
    if (flux_data$gastype[i] == "CO2") {
      # Subtract the best.flux values from NEE and RE
      flux_data$photosynthesis[i] <-  flux_data$best.flux[i+1] - flux_data$best.flux[i]
    }
  }
}

flux_data <- flux_data %>% select(...1, UniqueID, base_code, treatment, NEERE, date, best.flux, model, gastype, Animal, Campaign, plotID, longdate, plotNEERE, photosynthesis)

gradients <-  flux_data %>%  filter(Campaign == "Gradient", Animal == "Cow", gastype == "CH4")
dailies <- flux_data %>%  filter(Campaign == "Daily", Animal == "Horse", gastype == "CH4")


scatterplot1 <- ggplot(dailies, aes(x = UniqueID, y = best.flux))+
  geom_point(aes(color = treatment)) +
  geom_smooth(method = lm) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
 
scatterplot1

zoomdaily <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul"), Animal == "Horse")

scatterplot2 <- ggplot(zoomdaily, aes(x = UniqueID, y = best.flux))+
  geom_point(aes(color = treatment)) +
  geom_smooth(method = lm) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
scatterplot2

zoomdaily <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul"), Animal %in% c("Horse", "Cow"))

ggplot(zoomdaily, aes(x = date, y = best.flux, color = treatment)) +
  geom_point() +
  facet_wrap(~ Animal, scales = "free_y") +
  labs(x = "Date", y = "CO2 Flux")

### GENERATING PLOTS
# Create the plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# List of gases to plot
gases <- list(
  list(gastype = "CO2", y_label = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}), filename_suffix = "CO2"),
  list(gastype = "CH4", y_label = expression("nmol CH4 m"^{-2} * " s"^{-1}), filename_suffix = "CH4"),
  list(gastype = "N2O", y_label = expression("nmol N2O m"^{-2} * " s"^{-1}), filename_suffix = "N2O")
)

# Function to generate plots
generate_plots <- function(animal_type, campaign_type, date_filter, campaign_code) {
  for (gas in gases) {
    # Filter data based on gas type and other conditions
    gas_data <- flux_data %>%
      filter(Campaign == campaign_type, gastype == gas$gastype, date %in% date_filter, Animal == animal_type)
    
    # Check if the gas is CH4 or N2O, and filter for RE measurements only
    if (gas$gastype %in% c("CH4", "N2O")) {
      gas_data <- gas_data %>% filter(NEERE == "RE")  # Filter to include only RE measurements
    }
    
    # Print the number of rows in gas_data
    cat("Processing gas:", gas$gastype, " - Rows:", nrow(gas_data), "\n")
    
    # Check if there is data to plot
    if (nrow(gas_data) > 0) {
      # Create the plot dynamically
      gas_plot <- ggplot(gas_data, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
        geom_point(size = 2) +
        geom_line() +
        facet_wrap(~base_code, scales = "free") +
        labs(title = paste(campaign_type, animal_type, gas$gastype), x = "Date", y = gas$y_label,
             color = "Light/Dark",         # Changed title for the color legend
             shape = "Treatment") +
        scale_color_manual(values = c("NEE" = "gray", "RE" = "black")) +
        scale_shape_manual(values = c(16, 17)) +
        theme_minimal() +
        theme(
          legend.position.inside = c(0.75, 0.25),
          axis.line = element_line(color = "black"), 
          axis.text = element_text(color = "black"),   
          axis.ticks = element_line(color = "black"),
          panel.grid.major = element_blank(),               # No major grid lines
          panel.background = element_rect(fill = "white")   # White background
        )
      
      # Generate the filename dynamically
      filename <- file.path("plots", paste0(substr(animal_type, 1, 1), campaign_code, "_", gas$filename_suffix, ".jpeg"))
      
      # Save the plot
      ggsave(filename = filename, plot = gas_plot, width = 10, height = 8)
      
      cat("Saved plot to:", filename, "\n")  # Confirm save location
    } else {
      cat("No data available for gas:", gas$gastype, "\n")
    }
  }
}

# Daily plots for Horse and Cow
daily_date_filter <- c("15jul", "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug")
generate_plots("Horse", "Daily", daily_date_filter, "D")
generate_plots("Cow", "Daily", daily_date_filter, "D")

# Gradient plots for Horse and Cow
gradient_date_filter <- c("13jun", "14jun", "16jul", "17jul", "18jul", "29jul", "30jul")
generate_plots("Horse", "Gradient", gradient_date_filter, "G")
generate_plots("Cow", "Gradient", gradient_date_filter, "G")


# Photosynthesis plot # NEED TO ADJUST PHOTOSYNTHESIS IN DUNG PLOTS TO ACCOUNT FOR 
# VEGETATION AREA
photosynthesis_data <- flux_data %>% filter(gastype == "CO2", Campaign == "Daily", NEERE == "NEE")

photosynthesis_plot <- ggplot(photosynthesis_data, aes(x = Animal, y = photosynthesis, color = treatment)) +
  geom_boxplot() +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Cow vs Horse Photosynthesis") +
  theme_minimal()

photosynthesis_plot
ggsave(filename = "plots/HorseCow_photosynthesis.jpeg", plot = photosynthesis_plot, width = 10, height = 8)
