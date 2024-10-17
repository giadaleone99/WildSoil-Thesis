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

#Removing outlier gradient CH4
flux_data <- flux_data[-92, ]

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

first_gradient_date_1_3 <- as.Date("2024-06-13")
first_gradient_date_4_5 <- as.Date("2024-06-14")
first_daily_date_1_2 <- as.Date("2024-07-15")
first_daily_date_3_4 <- as.Date("2024-07-29")
first_HD5_date <- as.Date("2024-07-30")


flux_data <- flux_data %>%
  mutate(Days_Since_First = case_when(
    Campaign == "Gradient" ~ as.numeric(date_formatted - first_gradient_date_1_3),
    Campaign == "Gradient" ~ as.numeric(date_formatted - first_gradient_date_4_5),
    Campaign == "Daily" & base_code == "HD5" ~ as.numeric(date_formatted - first_HD5_date),
    Campaign == "Daily" & date_formatted <= first_daily_date_1_2 ~ as.numeric(date_formatted - first_daily_date_1_2),
    Campaign == "Daily" & date_formatted >= first_daily_date_1_2 & date_formatted < first_daily_date_3_4 ~ as.numeric(date_formatted - first_daily_date_1_2),  # Adjusted this condition
    Campaign == "Daily" & date_formatted >= first_daily_date_3_4 & !(base_code %in% c("HD1", "HD2", "CD1", "CD2")) ~ as.numeric(date_formatted - first_daily_date_3_4),  # Adjusted this condition
    Campaign == "Daily" & date_formatted >= first_daily_date_3_4 & base_code %in% c("HD1", "HD2", "CD1", "CD2") ~ as.numeric(date_formatted - first_daily_date_1_2),  # Adjusted for the old dailies
    TRUE ~ NA_real_  # Catch-all for any unexpected cases
  ))


flux_data_ANOVA <- flux_data %>% 
  filter(!(NEERE == "NEE" & gastype %in% c("CH4", "N2O"))) %>% 
  mutate(Unique_ANOVA = paste(plotNEERE, gastype, sep = "_")) %>% 
  filter(plottype != "PIT")


# DO NOT RUN -------------------------------------------------------------------------


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

# replace the best flux of NEE measurements of CO2 with photosynthesis flux
flux_data <- flux_data %>%
  mutate(best.flux = ifelse(!is.na(corr_photosynthesis), corr_photosynthesis, best.flux))

# Function to generate plots ### WORKING!!!! with DATES
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
      gas_plot <- ggplot(gas_data, aes(x = Days_Since_First, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
        geom_point(size = 2) +
        geom_line() +
        facet_wrap(~base_code, scales = "free") +
        labs(title = paste(campaign_type, animal_type, gas$gastype), x = "Date", y = gas$y_label,
             color = "RE/Photosynthesis",         # Changed title for the color legend
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

photosynthesis_plot <- ggplot(photosynthesis_data, aes(x = Animal, y = corr_photosynthesis, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Photosynthesis",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
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
library(lmerTest)
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
  facet_wrap(~ Campaign)+
  theme_minimal()

# Plot for CH4 and N2O
ggplot(subset(flux_data, gastype %in% c("CH4", "N2O")), aes(x = as.factor(Days_Since_First), y = best.flux, fill = Animal)) +
  geom_boxplot() +
  labs(x = "Days Since First", y = "Best Flux", title = "Box Plot of Best Flux for CH4 and N2O") +
  theme_minimal()


# Plot the subsets
plot_subsets(subsets_list)

gradient_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_CO2_PS_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "NEE") %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_CO2_RE_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "RE") %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_CH4_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "CH4") %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_N2O_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "N2O") %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_horse_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(Animal == "Horse") %>% 
  convert_as_factor(Days_Since_First) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_cow_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(Animal == "Cow") %>% 
  convert_as_factor(Days_Since_First) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  convert_as_factor(Days_Since_First) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_CO2_PS_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "NEE") %>%
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_CO2_RE_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "RE") %>%
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_CH4_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(gastype == "CH4") %>%
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_N2O_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(gastype == "N2O") %>%
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_horse_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(Animal == "Horse") %>% 
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

daily_cow_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(Animal == "Cow") %>% 
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_CH4_boxplot <- ggplot(gradient_CH4_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol CH4 m"^{-2} * " s"^{-1}),
       title = "Best flux CH4 Gradients",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

gradient_CH4_boxplot
#ggsave(filename = "veg_plots/gradient_CH4_boxplot.jpeg", plot = gradientvegweight, width = 6, height = 4)


# Gradient N2O boxplots
gradient_N2O_boxplot <- ggplot(gradient_N2O_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol N2O m"^{-2} * " s"^{-1}),
       title = "Best flux N2O Gradients",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

gradient_N2O_boxplot

# Gradient CO2 boxplots
gradient_CO2_boxplot_PS <- ggplot(gradient_CO2_PS_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Best flux CO2 Photosynthesis Gradients",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

gradient_CO2_boxplot_PS

gradient_CO2_boxplot_RE <- ggplot(gradient_CO2_RE_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Best flux CO2 RE Gradients",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

gradient_CO2_boxplot_RE

# Daily CO2 boxplots
daily_CO2_boxplot_PS <- ggplot(daily_CO2_PS_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Best flux CO2 Photosynthesis Daily",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

daily_CO2_boxplot_PS

daily_CO2_boxplot_RE <- ggplot(daily_CO2_RE_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Best flux CO2 RE Daily",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

daily_CO2_boxplot_RE

# Daily N2O boxplots
daily_N2O_boxplot <- ggplot(daily_N2O_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol N2O m"^{-2} * " s"^{-1}),
       title = "Best flux N2O Daily",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

daily_N2O_boxplot

# Daily CH4 boxplots
# Daily N2O boxplots
daily_CH4_boxplot <- ggplot(daily_CH4_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol CH4 m"^{-2} * " s"^{-1}),
       title = "Best flux CH4 Daily",
       colour = "Treatment",
       fill = "Treatment & Animal") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "bottom"
  )

daily_CH4_boxplot

#res <- anova_test(data = gradient_subset, dv = best.flux, wid = Unique_ANOVA, within = Days_Since_First)

# per gas type (cannot do per animal as well for the gradient data)
gradients_CO2_PS_ANOVA <- anova_test(data = gradient_CO2_PS_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Animal, treatment, period))
get_anova_table(gradients_CO2_PS_ANOVA)

gradients_CO2_RE_ANOVA <- anova_test(data = gradient_CO2_RE_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(gradients_CO2_RE_ANOVA)

gradients_CH4_ANOVA <- anova_test(data = gradient_CH4_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(gradients_CH4_ANOVA)

# N2O data for gradients too similar?
gradients_N2O_ANOVA <- anova_test(data = gradient_N2O_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(gradients_N2O_ANOVA)

# per animal probably not relevant cuz they are not per gas type
gradient_horse_ANOVA <- anova_test(data = gradient_horse_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(gradient_horse_ANOVA)

gradient_cow_ANOVA <- anova_test(data = gradient_cow_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(gradient_cow_ANOVA)

# dailies

# per gas type and animal
dailies_CO2_PS_ANOVA <- anova_test(data = daily_CO2_PS_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment, Animal))
get_anova_table(dailies_CO2_PS_ANOVA)

dailies_CO2_RE_ANOVA <- anova_test(data = daily_CO2_RE_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment, Animal))
get_anova_table(dailies_CO2_RE_ANOVA)

dailies_CH4_ANOVA <- anova_test(data = daily_CH4_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment, Animal))
get_anova_table(dailies_CH4_ANOVA)

dailies_N2O_ANOVA <- anova_test(data = daily_N2O_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment, Animal))
get_anova_table(dailies_N2O_ANOVA)

# per animal
daily_horse_ANOVA <- anova_test(data = daily_horse_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(daily_horse_ANOVA)

daily_cow_ANOVA <- anova_test(data = daily_cow_subset, dv = best.flux, wid = Unique_ANOVA, between = c(Days_Since_First, treatment))
get_anova_table(daily_cow_ANOVA)


# Test normality

# remove some data because there were only 2 measurements on this day since first
gradient_CO2_PS_subset <- gradient_CO2_PS_subset %>% 
  filter(Unique_ANOVA != "CG5_C_PS_CO2") %>% 
  filter(Unique_ANOVA != "CG5_F_PS_CO2")

gradient_CO2_RE_subset <- gradient_CO2_RE_subset %>% 
  filter(Unique_ANOVA != "CG5_C_RE_CO2") %>% 
  filter(Unique_ANOVA != "CG5_F_RE_CO2")

gradient_CH4_subset <- gradient_CH4_subset %>% 
  filter(Unique_ANOVA != "CG5_C_RE_CH4") %>% 
  filter(Unique_ANOVA != "CG5_F_RE_CH4")

gradient_N2O_subset <- gradient_N2O_subset %>% 
  filter(Unique_ANOVA != "CG5_C_RE_N2O") %>% 
  filter(Unique_ANOVA != "CG5_F_RE_N2O")

norm_daily_CO2_PS <- daily_CO2_PS_subset %>%
  shapiro_test(best.flux)

norm_daily_CO2_PS


norm_daily_CO2_RE <- daily_CO2_RE_subset %>%
  shapiro_test(best.flux)

norm_daily_CO2_RE

norm_daily_CH4 <- daily_CH4_subset %>%
  shapiro_test(best.flux)

norm_daily_CH4

norm_daily_N2O <- daily_N2O_subset %>%
  shapiro_test(best.flux)

norm_daily_N2O

# nothing is normal!!!
hist(daily_N2O_subset$best.flux, main="Histogram best flux")
hist(daily_CH4_subset$best.flux, main="Histogram best flux")
hist(daily_CO2_PS_subset$best.flux, main="Histogram best flux")
hist(daily_CO2_RE_subset$best.flux, main="Histogram best flux")

norm_gradient_CO2_PS <- gradient_CO2_PS_subset %>%
  group_by(Days_Since_First) %>% 
  shapiro_test(best.flux)

norm_gradient_CO2_PS

norm_gradient_CO2_RE <- gradient_CO2_RE_subset %>%
  shapiro_test(best.flux)

norm_gradient_CO2_RE

norm_gradient_CH4 <- gradient_CH4_subset %>%
  shapiro_test(best.flux)

norm_gradient_CH4

norm_gradient_N2O <- gradient_N2O_subset %>%
  shapiro_test(best.flux)

norm_gradient_N2O

# CO2 not normal 

data_wide <- gradient_CO2_RE_subset %>%
  pivot_wider(names_from = Unique_ANOVA, values_from = best.flux)

gradient_CO2_RE_subset$Unique_ANOVA <-factor(gradient_CO2_RE_subset$Unique_ANOVA)

gradient_CO2_RE_subset$Days_Since_First <- factor(gradient_CO2_RE_subset$Days_Since_First)


gradient_CO2_RE_subset$period <-  factor(gradient_CO2_RE_subset$period)

### WORKING
res.fried <- gradient_CO2_RE_subset %>% friedman_test(best.flux ~ period |Unique_ANOVA)
res.fried


result <- friedman.test(best.flux ~ treatment | plotID, data = gradient_CO2_RE_subset)

friedman_result <- friedman.test(as.matrix(data_wide[, -1]), groups = data_wide$UniqueID)
print(friedman_result)

res.fried



data.frame(normality)
residuals <- residuals(gradients_CH4_ANOVA)
qqnorm(residuals)

model <- aov(best.flux ~  treatment + Error(Unique_ANOVA/Days_Since_First), data = gradient_CO2_PS_subset)
res <- residuals(model)

print(model)
qqnorm(res)
qqline(res)
print(anova_results)

flux_data$gastype <- as.factor(flux_data$gastype) # CO2, CH4, N2O
flux_data$Animal <- as.factor(flux_data$Animal) # Horse, Cow
flux_data$treatment <- as.factor(flux_data$treatment) # Control or Treatment
flux_data$Days_Since_First <- as.numeric(flux_data$Days_Since_First)


model <- lmer(best.flux ~ gastype * Animal * treatment + (1 | UniqueID), data = flux_data)


model <- lmer(best.flux ~ treatment * period + 
                (treatment == "F") * Animal * period + 
                (1 | Unique_ANOVA), data = gradient_CO2_RE_subset)

lmm <- lmer(value)


summary(model)



residuals <- resid(model)
qqnorm(residuals)
qqline(residuals) 
shapiro.test(residuals)

plot(fitted(model), residuals)
abline(h = 0, col = "red") 

plot(model)

