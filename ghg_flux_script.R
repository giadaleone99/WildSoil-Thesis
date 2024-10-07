# script to process ghg flux data

library(stringr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(lubridate)
library(knitr)
library(gridExtra)
library(grid)
library(gtable)
library(lme4)
library(nlme)
library(emmeans)


# Import data 
flux_data_raw <- read.csv("flux_data/combined_data.csv")
dung_area_data <- read.csv("data/dung_soil_data.csv")

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

frame_area <- 3058.15 #cm2

flux_data <- flux_data %>% select(...1, UniqueID, base_code, treatment, NEERE, date, best.flux, model, gastype, Animal, Campaign, plotID, longdate, plotNEERE, photosynthesis)

flux_data <- left_join(flux_data, dung_area_data, by = "UniqueID")

flux_data <- flux_data %>% 
  mutate(corr_veg_dung = frame_area - dung_area_cm2,
         plottype = case_when(
           Campaign %in% c("Gradient", "Daily") ~ "Vegetated",  
           Campaign == "PIT" ~ "PIT",                          
           TRUE ~ NA_character_                                 
         ))

flux_data <- flux_data %>% 
  mutate(corr_photosynthesis = (frame_area/corr_veg_dung) * photosynthesis,
         period = case_when(
           longdate >= as.Date("2024-06-13") & longdate <= as.Date("2024-06-14") ~ 1,  # Period 1
           longdate >= as.Date("2024-07-15") & longdate <= as.Date("2024-07-19") ~ 2,  # Period 2
           longdate >= as.Date("2024-07-29") & longdate <= as.Date("2024-08-02") ~ 3,  # Period 3
           TRUE ~ NA_integer_  # Handle dates outside these ranges
         ))

flux_data <- flux_data %>% 
  mutate(Plotname_cleaned = sub("_[^_]*$", "", UniqueID),
         date_formatted = as.Date(paste0("2024-", 
                                         toupper(substr(date, 3, 5)), 
                                         "-", 
                                         substr(date, 1, 2)), 
                                  format = "%Y-%b-%d"))

first_gradient_date <- as.Date("2024-06-13")
first_daily_date_1_2 <- as.Date("2024-07-15")
first_daily_date_3_4 <- as.Date("2024-07-29")
first_HD5_date <- as.Date("2024-07-30")


flux_data <- flux_data %>%
  mutate(Days_Since_First = case_when(
    Campaign == "Gradient" ~ as.numeric(date_formatted - first_gradient_date),
    Campaign == "Daily" & base_code == "HD5" ~ as.numeric(date_formatted - first_HD5_date),
    Campaign == "Daily" & date_formatted <= first_daily_date_1_2 ~ as.numeric(date_formatted - first_daily_date_1_2),
    Campaign == "Daily" & date_formatted > first_daily_date_1_2 & date_formatted <= first_daily_date_3_4 ~ as.numeric(date_formatted - first_daily_date_1_2),  # Adjusted this condition
    Campaign == "Daily" & date_formatted > first_daily_date_3_4 ~ as.numeric(date_formatted - first_daily_date_3_4),  # Adjusted this condition
    TRUE ~ NA_real_  # Catch-all for any unexpected cases
  ))


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


# Photosynthesis plot 
photosynthesis_data <- flux_data %>% filter(gastype == "CO2", Campaign == "Daily", NEERE == "NEE")

photosynthesis_plot <- ggplot(photosynthesis_data, aes(x = Animal, y = corr_photosynthesis, color = treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(aes(color = treatment), position = position_dodge(width = 0.8))+
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Cow vs Horse Photosynthesis",
       colour = "Treatment") +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )
  

photosynthesis_plot
ggsave(filename = "plots/HorseCow_photosynthesis.jpeg", plot = photosynthesis_plot, width = 6, height = 5)


# Soil temp and SWC
SWC_plot <- ggplot(flux_data, aes(x = factor(period), y = SWC_., colour = plottype)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +  # Adjust boxplots
  geom_point(aes(color = plottype), position = position_dodge(width = 0.75), size = 2) +  # Align points with boxplots
  labs(x = "Sampling period", y = "Soil Water Content (%)", title = "Soil Water Content by Period",
       colour = "Plot type") +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

SWC_plot
ggsave("plots/SWC_plot.jpeg", plot = SWC_plot, width = 6, height = 5, dpi = 300, units = "in")

Stemp_plot <- ggplot(flux_data, aes(x = factor(period), y = S_temp, colour = plottype)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(aes(color = plottype), position = position_dodge(width = 0.75)) +
  labs(y = expression("Soil temperature"~(degree*C)), x = "Sampling period",
       title = "Soil temperature",
       colour = "Plot type") +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank()  # Removes outer border, if any
  )

Stemp_plot

ggsave("plots/Stemp_plot.jpeg", plot = Stemp_plot, width = 6, height = 5, dpi = 300, units = "in")

### SUMMARY STATISTICS THSI CRAP ISNT WORKING!! FIX NEXT TIME
Soil_temp_stats <- flux_data %>%
  group_by(period) %>%  # Group by 'period' column
  summarize(
    Mean_temp = mean(S_temp, na.rm = TRUE),  # Mean
    Median_temp = median(S_temp, na.rm = TRUE),  # Median
    Min_temp = min(S_temp, na.rm = TRUE),  # Minimum
    Max_temp = max(S_temp, na.rm = TRUE),  # Maximum
    SD_temp = sd(S_temp, na.rm = TRUE)  # Standard deviation
  )
print(Soil_temp_stats)

Soil_temp_stats[sapply(Soil_temp_stats, is.numeric)] <- lapply(Soil_temp_stats[sapply(Soil_temp_stats, is.numeric)], round, digits = 1)

# Use regular strings for column names in the data frame
colnames(Soil_temp_stats) <- c("Period", "Mean Soil Temp", "Median Soil Temp", "Min Soil Temp", "Max Soil Temp", "SD Soil Temp")

# Create the table with default labels
Soil_temp_table <- tableGrob(Soil_temp_stats)

# Save the table as a JPEG file
jpeg("tables/soiltemp_table.jpeg", width = 800, height = 400)  # Adjust dimensions as needed
grid.draw(Soil_temp_table)
dev.off()

# swc
SWC_stats <- flux_data %>%
  group_by(period) %>%  # Group by 'period' column
  summarize(
    Mean_SWC = mean(SWC_., na.rm = TRUE),  # Mean
    Median_SWC = median(SWC_., na.rm = TRUE),  # Median
    Min_SWC = min(SWC_., na.rm = TRUE),  # Minimum
    Max_SWC = max(SWC_., na.rm = TRUE),  # Maximum
    SD_SWC = sd(SWC_., na.rm = TRUE)  # Standard deviation
  )
print(SWC_stats)

SWC_stats[sapply(SWC_stats, is.numeric)] <- lapply(SWC_stats[sapply(SWC_stats, is.numeric)], round, digits = 1)

colnames(SWC_stats) <- colnames(SWC_stats) %>%
  gsub("_", " ", .) %>%                   # Replace underscores with spaces
  tools::toTitleCase() %>%                # Capitalize the first letter of each word
  gsub("Swc", "SWC (%)", ., ignore.case = TRUE)  # Append % to SWC columns

# Create the table
SWC_table <- tableGrob(SWC_stats)

# Save the table as a JPEG file
jpeg("tables/swc_table.jpeg", width = 800, height = 400)  # Adjust dimensions as needed
grid.draw(SWC_table)
dev.off()

# Repeated measures ANOVA
library(rstatix)
library(reshape)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(plyr)
library(datarium)





# Run repeated measures ANOVA
flux_data %>%
  group_by(Days_Since_First) %>%
  get_summary_stats(best.flux, type = "mean_sd")

flux_data$Days_Since_First <- as.factor(flux_data$Days_Since_First)



ggplot(flux_data, aes(x = gastype, y = best.flux, colour = Days_Since_First)) +
  geom_boxplot(aes(gastype))

# Plot for CO2
ggplot(subset(flux_data, gastype == "CO2"), aes(x = as.factor(Days_Since_First), y = best.flux, fill = Animal)) +
  geom_boxplot() +
  labs(x = "Days Since First", y = "Best Flux", title = "Box Plot of Best Flux for CO2") +
  facet_wrap(~ Campaign)
  theme_minimal()

# Plot for CH4 and N2O
ggplot(subset(flux_data, gastype %in% c("CH4", "N2O")), aes(x = as.factor(Days_Since_First), y = best.flux, fill = Animal)) +
  geom_boxplot() +
  labs(x = "Days Since First", y = "Best Flux", title = "Box Plot of Best Flux for CH4 and N2O") +
  theme_minimal()


# Subset data for aNOVA


# Function to create subsets based on specific conditions
create_subsets <- function(flux_data) {
  # Create a list to hold subsets
  subsets <- list()
  
  # Loop through combinations and create specific subsets
  if (any(flux_data$Animal == "Horse") && any(flux_data$Campaign == "Daily") && any(flux_data$gastype == "CO2")) {
    subsets[["Horse_Daily_CO2"]] <- flux_data %>%
      filter(Animal == "Horse", Campaign == "Daily", gastype == "CO2")
  }
  
  if (any(flux_data$Animal == "Horse") && any(flux_data$Campaign == "Daily") && any(flux_data$gastype == "CH4")) {
    subsets[["Horse_Daily_CH4"]] <- flux_data %>%
      filter(Animal == "Horse", Campaign == "Daily", gastype == "CH4")
  }
  
  if (any(flux_data$Animal == "Horse") && any(flux_data$Campaign == "Daily") && any(flux_data$gastype == "N2O")) {
    subsets[["Horse_Daily_N2O"]] <- flux_data %>%
      filter(Animal == "Horse", Campaign == "Daily", gastype == "N2O")
  }
  
  if (any(flux_data$Animal == "Cow") && any(flux_data$Campaign == "Daily") && any(flux_data$gastype == "CO2")) {
    subsets[["Cow_Daily_CO2"]] <- flux_data %>%
      filter(Animal == "Cow", Campaign == "Daily", gastype == "CO2")
  }
  
  if (any(flux_data$Animal == "Cow") && any(flux_data$Campaign == "Daily") && any(flux_data$gastype == "CH4")) {
    subsets[["Cow_Daily_CH4"]] <- flux_data %>%
      filter(Animal == "Cow", Campaign == "Daily", gastype == "CH4")
  }
  
  if (any(flux_data$Animal == "Cow") && any(flux_data$Campaign == "Daily") && any(flux_data$gastype == "N2O")) {
    subsets[["Cow_Daily_N2O"]] <- flux_data %>%
      filter(Animal == "Cow", Campaign == "Daily", gastype == "N2O")
  }
  
  # Add more conditions as needed for other combinations...
  
  return(subsets)
}

# Create subsets using your flux_data
subsets_list <- create_subsets(flux_data)

# Example of accessing a specific subset
# View a specific subset
print(subsets_list[["Horse_Daily_CO2"]])

# plot graphs


# Function to plot subsets
plot_subsets <- function(subsets) {
  for (subset_name in names(subsets)) {
    subset_data <- subsets[[subset_name]]
    
    # Create a plot
    p <- ggplot(subset_data, aes(x = Days_Since_First, y = best.flux)) +
      geom_boxplot() +
      geom_point()+
      labs(title = subset_name, x = "Days Since First Measurement", y = "best flux") +
      theme_minimal()
    
    # Print the plot
    print(p)
  }
}

# Plot the subsets
plot_subsets(subsets_list)





# Test normality
flux_data %>%
  group_by(Days_Since_First, gastype) %>%
  shapiro_test(best.flux)






print(anova_results)

flux_data$gastype <- as.factor(flux_data$gastype) # CO2, CH4, N2O
flux_data$Animal <- as.factor(flux_data$Animal) # Horse, Cow
flux_data$treatment <- as.factor(flux_data$treatment) # Control or Treatment
flux_data$Days_Since_First <- as.numeric(flux_data$Days_Since_First)


model <- lmer(best.flux ~ gastype * Animal * treatment + (1 | UniqueID), data = flux_data)

summary(model)




