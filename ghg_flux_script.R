# script to process ghg flux data

library(stringr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(lubridate)
library(knitr)
library(gridExtra)
library(patchwork)
library(grid)
library(gtable)
library(lme4)
library(nlme)
library(emmeans)
library(car)
library(glmmTMB)
library(rstatix)
library(Hmisc)


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

flux_data <- flux_data %>% dplyr::select(...1, UniqueID, base_code, treatment, NEERE, date, best.flux, model, gastype, Animal, Campaign, plotID, longdate, plotNEERE, photosynthesis)

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
    Campaign == "Gradient" & !(base_code %in% c("CG4", "CG5", "HG4", "HG5")) ~ as.numeric(date_formatted - first_gradient_date_1_3),
    Campaign == "Gradient" & base_code %in% c("CG4", "CG5", "HG4", "HG5") ~ as.numeric(date_formatted - first_gradient_date_4_5),
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

# replace the best flux of NEE measurements of CO2 with photosynthesis flux
flux_data_ANOVA <- flux_data_ANOVA %>%
  mutate(best.flux = ifelse(!is.na(corr_photosynthesis), corr_photosynthesis, best.flux))

flux_data <- flux_data %>%
  mutate(best.flux = ifelse(!is.na(corr_photosynthesis), corr_photosynthesis, best.flux))


# DO NOT RUN -------------------------------------------------------------------------


### GENERATING PLOTS
# Create the plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# List of gases to plot
gases <- list(
  list(gastype = "CO2", y_label = expression(mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}), filename_suffix = "CO2"),
  list(gastype = "CH4", y_label = expression("nmol CH"[4] * " m"^{-2} * " s"^{-1}), filename_suffix = "CH4"),
  list(gastype = "N2O", y_label = expression("nmol N"[2] * "O m"^{-2} * " s"^{-1}), filename_suffix = "N2O")
)

# Make the summarised plots
avgflux_data <- flux_data %>% 
  mutate(animain = substr(base_code, 1, nchar(base_code)-1)) %>% 
  filter(plottype == 'Vegetated') %>% 
  filter(!(gastype %in% c("CH4", "N2O") & NEERE == "NEE"))

## Cow GPP 
cowGPPgraph <- ggplot(avgflux_data[avgflux_data$gastype == "CO2" & avgflux_data$NEERE == "NEE" & avgflux_data$Animal == "Cow", ],
                      aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), , alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "none")+
  labs(title = 'A', x = 'Days since first measurement', y = expression(mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}),
       color = "Treatment", fill = "Treatment")

cowGPPgraph
ggsave(filename = "plots/cowGPPgraph.jpeg", plot = cowGPPgraph, width = 6, height = 4)

## Horse GPP
horseGPPgraph <- ggplot(avgflux_data[avgflux_data$gastype == "CO2" & avgflux_data$NEERE == "NEE" & avgflux_data$Animal == "Horse", ],
                        aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "right")+
  labs(title = 'B', x = 'Days since first measurement', y = element_blank(),
       color = "Treatment", fill = "Treatment")
horseGPPgraph
ggsave(filename = "plots/horseGPPgraph.jpeg", plot = horseGPPgraph, width = 6, height = 4)

cowhorseGPP <- cowGPPgraph + horseGPPgraph
cowhorseGPP
ggsave(filename = "plots/cowhorseGPPgraph.jpeg", plot = cowhorseGPP, width = 8, height = 4)

## Cow Reco graph
cowRecograph <- ggplot(avgflux_data[avgflux_data$gastype == "CO2" & avgflux_data$NEERE == "RE" & avgflux_data$Animal == "Cow", ],
                       aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "none")+
  labs(title = 'A', x = 'Days since first measurement', y = expression(mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}),
       color = "Treatment", fill = "Treatment")
cowRecograph
ggsave(filename = "plots/cowRecograph.jpeg", plot = cowRecograph, width = 6, height = 4)

## Horse Reco graph
horseRecograph <- ggplot(avgflux_data[avgflux_data$gastype == "CO2" & avgflux_data$NEERE == "RE" & avgflux_data$Animal == "Horse", ],
                        aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "right")+
  labs(title = 'B', x = 'Days since first measurement', y = element_blank(),
       color = "Treatment", fill = "Treatment")
horseRecograph
ggsave(filename = "plots/horseRecograph.jpeg", plot = horseRecograph, width = 6, height = 4)

CowhorseRecograph <- cowRecograph + horseRecograph
CowhorseRecograph
ggsave(filename = "plots/CowhorseRecograph.jpeg", plot = CowhorseRecograph, width = 8, height = 4)

## Cow CH4 graph
cowCH4graph <- ggplot(avgflux_data[avgflux_data$gastype == "CH4" & 
                                     avgflux_data$Animal == "Cow", ],
                      aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
    scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  labs(title = 'A', 
       x = 'Days since first measurement', 
       y = expression("nmol CH"[4] * " m"^{-2} * " s"^{-1}),
       color = "Treatment", fill = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "none")

cowCH4graph
ggsave(filename = "plots/cowCH4graph.jpeg", plot = cowCH4graph, width = 6, height = 4)

## Horse CH4 graph
horseCH4graph <- ggplot(avgflux_data[avgflux_data$gastype == "CH4" & avgflux_data$Animal == "Horse", ],
                        aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "right")+
  labs(title = 'B', x = 'Days since first measurement', y = element_blank(),
       color = "Treatment", fill = "Treatment")
horseCH4graph
ggsave(filename = "plots/horseCH4graph.jpeg", plot = horseCH4graph, width = 6, height = 4)

CowhorseCH4graph <- cowCH4graph + horseCH4graph
CowhorseCH4graph
ggsave(filename = "plots/CowhorseCH4graph.jpeg", plot = CowhorseCH4graph, width = 8, height = 4)


## Cow N2O 
cowN2Ograph <- ggplot(avgflux_data[avgflux_data$gastype == "N2O" & avgflux_data$Animal == "Cow", ],
                      aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "none")+
  labs(title = 'A', x = 'Days since first measurement', y = expression("nmol N"[2] * "O m"^{-2} * " s"^{-1}),
       color = "Treatment", fill = "Treatment")
cowN2Ograph
ggsave(filename = "plots/cowN2Ograph.jpeg", plot = cowN2Ograph, width = 6, height = 4)

horseN2Ograph <- ggplot(avgflux_data[avgflux_data$gastype == "N2O" & avgflux_data$Animal == "Horse", ],
                        aes(x = factor(Days_Since_First), y = best.flux, color = treatment)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5), alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "point",
               aes(group = treatment, color = treatment),
               position = position_dodge(width = 1), 
               size = 2.5) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position = position_dodge(width = 1), 
               width = 0.5) +
  scale_color_manual(values = c("C" = "grey", "F" = "black"), 
                     labels = c("Control", "Dung")) +
  scale_x_discrete(breaks = levels(factor(avgflux_data$Days_Since_First))) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "right")+
  labs(title = 'B', x = 'Days since first measurement', y = element_blank(),
       color = "Treatment", fill = "Treatment")
horseN2Ograph
ggsave(filename = "plots/horseN2Ograph.jpeg", plot = horseN2Ograph, width = 6, height = 4)

CowhorseN2Ograph <- cowN2Ograph + horseN2Ograph
CowhorseN2Ograph
ggsave(filename = "plots/CowhorseN2Ograph.jpeg", plot = CowhorseN2Ograph, width = 8, height = 4)


generate_avgplots <- function(animal_type, campaign_type, date_filter, campaign_code) {
  for (gas in gases) {
    # Filter data based on gas type and other conditions
    gas_data <- flux_data %>%
      filter(Campaign == campaign_type, gastype == gas$gastype, date %in% date_filter, Animal == animal_type)
    
    # Check if the gas is CH4 or N2O, and filter for RE measurements only
    if (gas$gastype %in% c("CH4", "N2O")) {
      gas_data <- gas_data %>% filter(NEERE == "RE")  # Filter to include only RE measurements
    }
    
    #summarise the data
    
    # Print the number of rows in gas_data
    cat("Processing gas:", gas$gastype, " - Rows:", nrow(gas_data), "\n")
    
    # Check if there is data to plot
    if (nrow(gas_data) > 0) {
      # Create the plot dynamically
      gas_plot <- ggplot(gas_data, aes(x = Days_Since_First, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
        geom_point(size = 2) +
        geom_line() +
        facet_wrap(~base_code, scales = "free") +
        labs(title = paste(campaign_type, animal_type, gas$gastype), x = "Days since first measurement", y = gas$y_label,
             color = "Measurement type",         # Changed title for the color legend
             shape = "Treatment") +
        scale_color_manual(values = c("NEE" = "gray", "RE" = "black"),
                           labels = c("Photosynthesis",
                                      "Respiration")) +
        scale_shape_manual(values = c("C" = 16, "F" = 17),  # 16 and 17 are shape codes for circles and triangles
                           labels = c("Control", "Dung")) +
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
        labs(title = paste(campaign_type, animal_type, gas$gastype), x = "Days since first measurement", y = gas$y_label,
             color = "Measurement type",         # Changed title for the color legend
             shape = "Treatment") +
        scale_color_manual(values = c("NEE" = "gray", "RE" = "black"),
                           labels = c("Photosynthesis",
                                      "Respiration")) +
        scale_shape_manual(values = c("C" = 16, "F" = 17),  # 16 and 17 are shape codes for circles and triangles
                           labels = c("Control", "Dung")) +
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
photosynthesis_data <- flux_data %>% filter(gastype == "CO2", Campaign %in% c("Daily", "Gradient"), NEERE == "NEE")

photosynthesis_plot <- ggplot(photosynthesis_data, aes(x = Animal, y = corr_photosynthesis, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  facet_wrap(~Campaign) +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
       title = "Photosynthesis",
       colour = "Treatment",
       fill = "Plot type") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh")) +
  theme_minimal() +
  scale_y_continuous() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    panel.grid.major.x = element_blank()
  )
  

photosynthesis_plot
ggsave(filename = "plots/HorseCow_photosynthesis.jpeg", plot = photosynthesis_plot, width = 6, height = 5)


# Soil temp and SWC
SWC_plot <- ggplot(flux_data, aes(x = factor(period), y = SWC_., colour = plottype)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +  # Adjust boxplots
  geom_point(aes(color = plottype), position = position_dodge(width = 0.75), size = 2) +  # Align points with boxplots
  labs(x = "Sampling period", y = "Soil Water Content (%)", title = "Soil Water Content by Period",
       colour = "Plot type") +
  scale_color_manual(values = c("Vegetated" = "black", "PIT" = "gray"))+
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
  scale_color_manual(values = c("Vegetated" = "black", "PIT" = "gray"))+
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank()  # Removes outer border, if any
  )

Stemp_plot

ggsave("plots/Stemp_plot.jpeg", plot = Stemp_plot, width = 6, height = 5, dpi = 300, units = "in")

### SUMMARY STATISTICS THSI CRAP ISNT WORKING!! FIX NEXT TIME ADD Standard ERROR!
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

### RUN FROM HERE -----------------------------------------------------------------------

# Repeated measures ANOVA
library(rstatix)
library(ggbreak)
library(reshape)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(plyr)
library(datarium)
library(lmerTest)
library(bestNormalize)
library(plotrix)
library(DHARMa)


# Plot with all the gasses and stuff
ggplot(flux_data, aes(x = gastype, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(aes(gastype)) +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24")) +
  theme_minimal()

# Plot for CO2
ggplot(subset(flux_data, gastype == "CO2"), aes(y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot() +
  labs(x = "Days Since First", y = "Best Flux", title = "Box Plot of Best Flux for CO2") +
  facet_wrap(~ Campaign)+
  theme_minimal() +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"))

# Plot for CH4 and N2O
ggplot(subset(flux_data, gastype %in% c("CH4", "N2O")), aes(x = as.factor(Days_Since_First), y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot() +
  facet_wrap(~ gastype) +
  labs(x = "Days Since First", y = "Best Flux", title = "Box Plot of Best Flux for CH4 and N2O") +
  theme_minimal() +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"))

### CREATING SUBSETS -------------------------------------------

## SUBSETS WITH BOTH CAMPAIGNS
CO2_PS_subset <- flux_data_ANOVA %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "NEE") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

CO2_RE_subset <- flux_data_ANOVA %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "RE") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         log_best.flux = log(best.flux))

CH4_subset <- flux_data_ANOVA %>% 
  filter(gastype == "CH4") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>%
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         normalized_best.flux = bestNormalize(best.flux)$x.t) %>% 
  mutate(ranked_best.flux = rank(best.flux))

hist(CH4_subset$ranked_best.flux, breaks = 20)

N2O_subset <- flux_data_ANOVA %>% 
  filter(gastype == "N2O") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         normalized_best.flux = bestNormalize(best.flux)$x.t) %>% 
  mutate(ranked_best.flux = rank(best.flux)) %>% 
  mutate(logged_best.flux = log(0.48 + best.flux))

## GRADIENT SUBSETS --------------------
gradient_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_CO2_PS_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "NEE") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA))

gradient_CO2_RE_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "CO2") %>%
  filter(NEERE == "RE") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
        log_best.flux = log(best.flux))

gradient_CH4_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "CH4") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>%
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         normalized_best.flux = bestNormalize(best.flux)$x.t)

gradient_N2O_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Gradient") %>% 
  filter(gastype == "N2O") %>%
  convert_as_factor(period) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         normalized_best.flux = bestNormalize(best.flux)$x.t)

## DAILY SUBSETS -----------------

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
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         normalized_best.flux = bestNormalize(best.flux)$x.t)

daily_N2O_subset <- flux_data_ANOVA %>% 
  filter(Campaign == "Daily") %>% 
  filter(gastype == "N2O") %>%
  filter(!(
    base_code %in% c("CD1", "CD2", "HD1", "HD2") | as.integer(Days_Since_First) > 10
  )) %>%
  convert_as_factor(Days_Since_First) %>% 
  convert_as_factor(gastype) %>% 
  convert_as_factor(treatment) %>% 
  convert_as_factor(Animal) %>% 
  mutate(Unique_ANOVA = as.factor(Unique_ANOVA),
         normalized_best.flux = log(0.4 + best.flux))

# summary stats

gradient_CO2_PS_summary <- gradient_CO2_PS_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
gradient_CO2_RE_summary <- gradient_CO2_RE_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
gradient_CH4_summary <- gradient_CH4_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
gradient_N2O_summary <- gradient_N2O_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
daily_CO2_PS_summary <- daily_CO2_PS_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
daily_CO2_RE_summary <- daily_CO2_RE_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
daily_CH4_summary <- daily_CH4_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
daily_N2O_summary <- daily_N2O_subset %>%
  group_by(Animal, treatment) %>% 
  dplyr::summarise(
    Nperc_mean = mean(best.flux, na.rm = TRUE),
    Nperc_se = std.error(best.flux),)
### PLOTTING SUBSETS --------------------------------------------------

## Gradients
gradient_CH4_boxplot <- ggplot(gradient_CH4_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5) +  # Suppress outliers in the boxplot
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol CH"[4] * " m"^{-2} * " s"^{-1}),
       title = "B",
       colour = "Treatment",
       fill = "Plot type") +
  scale_y_break(c(4, 25),c(30, 115), scales = c(0.2)) + 
  # Breaks with scaling
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank()
  )

gradient_CH4_boxplot
ggsave(filename = "plots/gradient_CH4_boxplot.jpeg", plot = gradient_CH4_boxplot, width = 6, height = 5)

# Gradient N2O boxplots
gradient_N2O_boxplot <- ggplot(gradient_N2O_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol N"[2] * "O m"^{-2} * " s"^{-1}),
       title = "B",
       colour = "Treatment",
       fill = "Plot type") +
  scale_y_break(c(0.8, 1.15), scales = (0.4)) +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank()
  )

gradient_N2O_boxplot
ggsave(filename = "plots/gradient_N2O_boxplot.jpeg", plot = gradient_N2O_boxplot, width = 6, height = 5)

# Gradient CO2 boxplots
gradient_CO2_boxplot_PS <- ggplot(gradient_CO2_PS_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO"[2] *" m"^{-2} * " s"^{-1}),
       title = "B",
       colour = "Treatment",
       fill = "Plot type") +
  scale_y_continuous(limits = c(-1, 36)) +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank()
  )

gradient_CO2_boxplot_PS
ggsave(filename = "plots/gradient_CO2_PS_boxplot.jpeg", plot = gradient_CO2_boxplot_PS, width = 6, height = 5)

gradient_CO2_boxplot_RE <- ggplot(gradient_CO2_RE_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO"[2] *" m"^{-2} * " s"^{-1}),
       title = "B",
       colour = "Treatment",
       fill = "Plot type") +
  ylim(0,30) +
  scale_y_continuous(limits = c(4, 30)) +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank()
  )

gradient_CO2_boxplot_RE
ggsave(filename = "plots/gradient_CO2_RE_boxplot.jpeg", plot = gradient_CO2_boxplot_RE, width = 6, height = 5)
## DAILY BOXPLOTS -------------------------------------------------

# Daily CO2 boxplots
daily_CO2_boxplot_PS <- ggplot(daily_CO2_PS_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression(mu * "mol CO"[2] *" m"^{-2} * " s"^{-1}),
       title = "A",
       colour = "Treatment",
       fill = "Plot type") +
  scale_y_continuous(limits = c(-1, 36)) +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "none"
  )

daily_CO2_boxplot_PS
ggsave(filename = "plots/daily_CO2_PS_boxplot.jpeg", plot = daily_CO2_boxplot_PS, width = 6, height = 5)

daily_CO2_boxplot_RE <- ggplot(daily_CO2_RE_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  scale_y_continuous(limits = c(4, 30)) +
  labs(y = expression(mu * "mol CO"[2] *" m"^{-2} * " s"^{-1}),
       title = "A",
       colour = "Treatment",
       fill = "Plot type") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "none"
  )

daily_CO2_boxplot_RE
ggsave(filename = "plots/daily_CO2_RE_boxplot.jpeg", plot = daily_CO2_boxplot_RE, width = 6, height = 5)

# Daily N2O boxplots
daily_N2O_boxplot <- ggplot(daily_N2O_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol N"[2] * "O m"^{-2} * " s"^{-1}),
       title = "A",
       colour = "Treatment",
       fill = "Plot type") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "none"
  )

daily_N2O_boxplot
ggsave(filename = "plots/daily_N2O_boxplot.jpeg", plot = daily_N2O_boxplot, width = 6, height = 5)

# Daily CH4 boxplots
daily_CH4_boxplot <- ggplot(daily_CH4_subset, aes(x = Animal, y = best.flux, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8)) +
  labs(y = expression("nmol CH"[4] * " m"^{-2} * " s"^{-1}),
       title = "A",
       colour = "Treatment",
       fill = "Plot type") +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "F.Cow" = "#656D4A", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Adds axis lines for both x and y
    panel.border = element_blank(),
    legend.position = "none"
  )

daily_CH4_boxplot

combined_plot_data <- daily_CH4_subset[daily_CH4_subset$treatment %in% c("C", "F") & 
                                         !(daily_CH4_subset$Animal == "Cow" & daily_CH4_subset$treatment == "F"), ]

# Combined plot for Cow control, Horse control, and Fresh
NoCow_CH4_plot <- ggplot(combined_plot_data, aes(x = interaction(treatment, Animal), y = best.flux)) +
  geom_boxplot(aes(fill = interaction(treatment, Animal)), position = position_dodge(), width = 0.7) +
  geom_jitter(position = position_dodge(width = 0.8)) +  # Jittered points
  labs(x = "Animal",  # Change the x-axis label here
       y = expression("nmol CH"[4] * " m"^{-2} * " s"^{-1}),
       title = "A",
       fill = "Plot type") +
  scale_y_break(c(1, 40), scales = (0.3)) +
  scale_fill_manual(values = c("C.Cow" = "#A4AC86", 
                               "C.Horse" = "#A68A64", 
                               "F.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Horse control", "Horse dung")) +
  theme_minimal() +
  #scale_x_discrete(labels = c("C.Cow" = "Cow", "C.Horse" = "Horse", "F.Horse" = "Horse")) +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        legend.position = "none")

# Display the plot
NoCow_CH4_plot
ggsave(filename = "plots/NoCow_CH4_boxplot.jpeg", plot = NoCow_CH4_plot, width = 6, height = 5)

# Separate plot for Cow fresh
# Subset for Cow Fresh data
cow_fresh_plot_data <- daily_CH4_subset[daily_CH4_subset$treatment == "F" & daily_CH4_subset$Animal == "Cow", ]

cow_fresh_CH4_plot <- ggplot(cow_fresh_plot_data, aes(x = treatment, y = best.flux)) +
  geom_boxplot(aes(fill = interaction(treatment, Animal)), position = position_dodge(), width = 0.25) +
  geom_jitter(position = position_dodge(width = 0.8)) + 
  labs(x = "Animal",
       y = expression("nmol CH"[4] * " m"^{-2} * " s"^{-1}),
       title = "",
       fill = "Plot type") +
  scale_fill_manual(values = c("F.Cow" = "#656D4A"),
                    labels = c("Cow dung")) +
  scale_x_discrete(labels = c("F" = "Cow")) +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

cow_fresh_CH4_plot
ggsave(filename = "plots/Cow_fresh_CH4_boxplot.jpeg", plot = cow_fresh_CH4_plot, width = 4, height = 4.5)

# combine the plots per campaign
CH4_boxplots = NoCow_CH4_plot + cow_fresh_CH4_plot + gradient_CH4_boxplot
N2O_boxplots = daily_N2O_boxplot + gradient_N2O_boxplot
CO2_RE_boxplots = daily_CO2_boxplot_RE + gradient_CO2_boxplot_RE
CO2_PS_boxplots = daily_CO2_boxplot_PS + gradient_CO2_boxplot_PS

# CH4_boxplots
# N2O_boxplots
# CO2_RE_boxplots
# CO2_PS_boxplots

ggsave(filename = "plots/CH4_boxplots.jpeg", plot = CH4_boxplots, width = 10, height = 5)
ggsave(filename = "plots/N2O_boxplots.jpeg", plot = N2O_boxplots, width = 10, height = 5)
ggsave(filename = "plots/CO2_RE_boxplots.jpeg", plot = CO2_RE_boxplots, width = 10, height = 5)
ggsave(filename = "plots/CO2_PS_boxplots.jpeg", plot = CO2_PS_boxplots, width = 10, height = 5)
# MODELLING ----------------------------------------------------
# modelling Lasses way and combining the campaigns
# HERE
CH4_subset$Days_Since_First <- as.factor(CH4_subset$Days_Since_First)
CH4_subset$Campaign <- as.factor(CH4_subset$Campaign)

CO2_PS_model <- glmmTMB(best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), data = CO2_PS_subset)
CO2_RE_model <- glmmTMB(best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), family = tweedie, data = CO2_RE_subset)
CH4_model1 <- glmmTMB(ranked_best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), 
                      data = CH4_subset)
CH4_model2 <- glmmTMB(normalized_best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), 
                      data = CH4_subset)
N2O_model1 <- glmmTMB(ranked_best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), 
                     data = N2O_subset)
N2O_model2 <- glmmTMB(normalized_best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), 
                      data = N2O_subset)
N2O_model3 <- glmmTMB(logged_best.flux ~ Animal * treatment * Campaign + (1|Days_Since_First), 
                      data = N2O_subset)

run_model <- function(dataset, model) {
  #print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  #testDispersion(simulationOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$treatment)
  plotResiduals(simuOutput, form = dataset$Campaign)
  test <- emmeans(model, ~ Campaign|treatment|Animal)
  contrast(test, method = "pairwise") %>% as.data.frame()
}
run_model(N2O_subset, N2O_model3)


# creating the models for all the subsets
gradient_CO2_PS_model <- lmer(best.flux ~ Animal + treatment + (1|period), data = gradient_CO2_PS_subset) #best model to use. using * instead of + was not significant when comparing with anova(), so we use the more simple model
gradient_CO2_RE_model <- lmer(log_best.flux ~ Animal + treatment + (1|period), data = gradient_CO2_RE_subset) #same thing, use +
gradient_CO2_RE_model <- art(best.flux ~ Animal * treatment + (1|period), data = gradient_CO2_RE_subset)
anova(gradient_CO2_RE_model)
gradient_CH4_model <- lmer(normalized_best.flux ~ Animal * treatment + (1|period), data = gradient_CH4_subset) #this time the anova() comparision was significant, so use *
gradient_CH4_model <- art(best.flux ~ Animal * treatment + (1|period), data = gradient_CH4_subset)
anova(gradient_CH4_model)
gradient_N2O_model <- lmer(normalized_best.flux ~ Animal + treatment + (1|period), data = gradient_N2O_subset) #not significant so use +
gradient_N2O_model <- art(best.flux ~ Animal * treatment + (1|period), data = gradient_N2O_subset) 
anova(gradient_N2O_model)

daily_CO2_PS_model <- lmer(best.flux ~ Animal + treatment + (1|Days_Since_First), data = daily_CO2_PS_subset) #not significant so use +
daily_CO2_RE_model <- lmer(best.flux ~ Animal + treatment + (1|Days_Since_First), data = daily_CO2_RE_subset) #not significant so use +
daily_CH4_model <- lmer(normalized_best.flux ~ Animal * treatment + (1|Days_Since_First), data = daily_CH4_subset) #significant so use *
daily_CH4_model <- art(best.flux ~ Animal * treatment + (1|Days_Since_First), data = daily_CH4_subset)
anova(daily_CH4_model)
daily_N2O_model <- lmer(normalized_best.flux ~ Animal + treatment + (1|Days_Since_First), data = daily_N2O_subset) #not significant so use +
daily_N2O_model <- art(best.flux ~ Animal * treatment + (1|Days_Since_First), data = daily_N2O_subset)
anova(daily_N2O_model)

res <- residuals(gradient_CO2_RE_model)

plot(res)
qqnorm(res)
shapiro.test(res)
hist(res)
summary(gradient_CO2_RE_model)

plot(simulateResiduals(gradient_CO2_PS_model1))
AIC(gradient_CO2_RE_model)
plot(simulateResiduals(gradient_CO2_PS_model2))
AIC(gradient_CO2_PS_model2)
anova(gradient_CO2_PS_model, model)

# alternative model from chatGPT, bit sus
model <- lmer(best.flux ~ treatment + 
                (treatment == "F") * Animal  + 
                (1 | period), data = gradient_CO2_PS_subset)

# take out all Controls for all the subsets to test if animal becomes significant
gradient_CO2_PS_subset_F <- gradient_CO2_PS_subset %>% 
  filter(treatment == "F")
gradient_CO2_RE_subset_F <- gradient_CO2_RE_subset %>% 
  filter(treatment == "F")
gradient_CH4_subset_F <- gradient_CH4_subset %>% 
  filter(treatment == "F")
gradient_N2O_subset_F <- gradient_N2O_subset %>% 
  filter(treatment == "F")
daily_CO2_PS_subset_F <- daily_CO2_PS_subset %>% 
  filter(treatment == "F")
daily_CO2_RE_subset_F <- daily_CO2_RE_subset %>% 
  filter(treatment == "F")
daily_CH4_subset_F <- daily_CH4_subset %>% 
  filter(treatment == "F")
daily_N2O_subset_F <- daily_N2O_subset %>% 
  filter(treatment == "F")

gradient_CO2_PS_model <- lmer(best.flux ~ Animal + (1|period), data = gradient_CO2_PS_subset_F)
gradient_CO2_RE_model <- lmer(best.flux ~ Animal + (1|period), data = gradient_CO2_RE_subset_F)
gradient_CH4_model <- lmer(normalized_best.flux ~ Animal + (1|period), data = gradient_CH4_subset_F)
gradient_N2O_model <- lmer(normalized_best.flux ~ Animal + (1|period), data = gradient_N2O_subset_F)

daily_CO2_PS_model <- lmer(best.flux ~ Animal + (1|Days_Since_First), data = daily_CO2_PS_subset_F)
daily_CO2_RE_model <- lmer(best.flux ~ Animal + (1|Days_Since_First), data = daily_CO2_RE_subset_F)
daily_CH4_model <- lmer(normalized_best.flux ~ Animal + (1|Days_Since_First), data = daily_CH4_subset_F)
daily_N2O_model <- lmer(normalized_best.flux ~ Animal + (1|Days_Since_First), data = daily_N2O_subset_F)

### Linear Mixed Models ---------------------------------------------------------------
## Comment from stats dude: if it is a gradient subset, use the period as the random effect, 
## if it is a daily, then use the days_since_first as the random effect.

## Daily models --------------------------------------
# Daily CO2 PS subset

model <- lmer(normalized_best.flux ~ treatment + Animal + (1 | Days_Since_First), data = daily_CO2_PS_subset)
summary(model)
anova(model)
residuals <- residuals(model)
plot(residuals)
qqnorm(residuals)
qqline(residuals)
shapiro.test(residuals)

model <- lmer(normalized_best.flux ~ treatment + 
                (treatment == "F") * Animal  + 
                (1 | Days_Since_First), data = daily_CO2_PS_subset)


# Daily CO2 RE subset
model <- lmer(normalized_best.flux ~ treatment + Animal + (1 | Days_Since_First), data = daily_CO2_RE_subset)
summary(model)
anova(model)
residuals <- residuals(model)
plot(residuals)
qqnorm(residuals)
qqline(residuals)
shapiro.test(residuals)

model <- lmer(normalized_best.flux ~ treatment + 
                (treatment == "F") * Animal  + 
                (1 | Days_Since_First), data = daily_CO2_RE_subset)

summary(model)

# Daily CH4 subset
model <- lmer(normalized_best.flux ~ treatment + Animal + (1 | Days_Since_First), data = daily_CH4_subset)
summary(model)
anova(model)
residuals <- residuals(model)
plot(residuals)
qqnorm(residuals)
qqline(residuals)
shapiro.test(residuals)

model <- lmer(normalized_best.flux ~ treatment + 
                (treatment == "F") * Animal  + 
                (1 | Days_Since_First), data = daily_CH4_subset)

summary(model)


# Daily N2O subset
model <- lmer(normalized_best.flux ~ treatment + Animal + (1 | Days_Since_First), data = daily_N2O_subset)
summary(model)
anova(model)
residuals <- residuals(model)
plot(residuals)
qqnorm(residuals)
qqline(residuals)
shapiro.test(residuals)

model <- lmer(normalized_best.flux ~ treatment + 
                (treatment == "F") * Animal  + 
                (1 | Days_Since_First), data = daily_N2O_subset)

summary(model)

#### GRADIENT models ----------------------------------------------------------------------

## Gradient PS

## Gradient RE

## Gradient CH4

# Gradient N2O


### Consider running the models without control treatment rows and just look at Animal (cow, horse)

# remove trash and save for combined analysis
gas_data <- flux_data_ANOVA %>% 
  subset(select = -c(model, photosynthesis, plotNEERE, date, plottype, corr_photosynthesis, date_formatted, Plotname_cleaned, Unique_ANOVA))

saveRDS(gas_data, file = "flux_data/clean_flux_data.rds")
