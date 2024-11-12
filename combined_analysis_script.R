library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(interactions)

# Import data
allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
veg_growth_data <- readRDS("plant_data/veg_growth_data.rds") %>% 
  rename("plotID"= "plot_id") %>% 
  dplyr::select(-2:-9)
soil_data_raw <- read.csv("data/soil_data_raw.csv")
veg_combined_data <- readRDS("plant_data/veg_combined_data.rds")%>% 
  rename("plotID"= "plot_id",
         "plant_CN" = "CN_ratio") %>% 
  dplyr::select(-5)
  

# Combine data
CO2_PS_data <- allflux_data %>% 
  filter(gastype == "CO2" & NEERE == "NEE") %>% 
  dplyr::rename(CO2_PS_flux = best.flux)
CO2_RE_data <- allflux_data %>% 
  filter(gastype == "CO2" & NEERE == "RE") %>%
  dplyr::rename(CO2_RE_flux = best.flux)
CH4_data <- allflux_data %>% 
  filter(gastype == "CH4") %>% 
  dplyr::rename(CH4_flux = best.flux)
N2O_data <- allflux_data %>% 
  filter(gastype == "N2O") %>% 
  dplyr::rename(N2O_flux = best.flux)
# Join gas data per date
flux_data <- CO2_PS_data %>% 
  left_join(CO2_RE_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-18:-22, -24:-32)
flux_data <- flux_data %>% 
  left_join(CH4_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-19:-23, -25:-33)
flux_data <- flux_data %>% 
  left_join(N2O_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-1:-5, -7:-9, -12:-17, -20, -26)
flux_data <- flux_data[, c(6, 7, 2, 8, 9, 11, 12, 13, 3, 17, 18, 1, 4, 5, 10, 14, 15, 16)]

# create a df summing the gas data on the diff days for each plot
plot_flux_data <- flux_data %>% 
  group_by(plotID) %>% 
  dplyr::summarise(CO2_PS_sum = sum(CO2_PS_flux), 
                   CO2_RE_sum = sum(CO2_RE_flux), 
                   CH4_sum = sum(CH4_flux), 
                   N2O_sum = sum(N2O_flux), .groups = 'drop')

# join with veg growth
plot_flux_data <- plot_flux_data %>% 
  left_join(veg_growth_data, by = "plotID")


soil_data <- soil_data_raw %>% 
  filter(sample_type %in% c("Fresh", "Control"))

plot_flux_data <- plot_flux_data %>% 
  left_join(veg_combined_data, by = "plotID") %>% 
  left_join(soil_data, by = "plotID")

growth_model_data <- plot_flux_data %>% 
  filter(!is.na(height_value))

mod <- glmmTMB(P ~ PO4.P, data = growth_model_data)
summary(mod)
print(Anova(mod))

model <- glmmTMB(height_value ~ treatment * Animal * CN_ratio * PO4.P, data = growth_model_data)
summary(model)

  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  testDispersion(simuOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = growth_model_data$Animal)
  plotResiduals(simuOutput, form = growth_model_data$treatment)
  plotResiduals(simuOutput, form = growth_model_data$CN_ratio)
  plotResiduals(simuOutput, form = growth_model_data$PO4.P)
  plotResiduals(simuOutput, form = growth_model_data$P)
  test <- emmeans(model, ~ Animal|treatment*PO4.P*CN_ratio)
  contrast(test, method = "pairwise") %>% as.data.frame()
  interactions::interact_plot(model, pred = PO4.P, modx = Animal)

# Run gas flux models with soil moisture and temperature as random factors
CO2_PS_model <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data)
CO2_RE_model <- glmmTMB(CO2_RE_flux ~ Animal * treatment * Campaign * SWC_. + (1|Days_Since_First), data = flux_data)

run_model <- function(dataset, model) {
  #print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  #testDispersion(simulationOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$treatment)
  plotResiduals(simuOutput, form = dataset$Campaign)
  plotResiduals(simuOutput, form = dataset$S_temp)
  test <- emmeans(model, ~ Campaign|treatment|Animal|S_temp)
  contrast(test, method = "pairwise") %>% as.data.frame()
}

run_model(flux_data, CO2_RE_model)
