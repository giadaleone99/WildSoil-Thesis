library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)

# Import data
allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")

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

# Run gas flux models with soil moisture and temperature as random factors
CO2_PS_model <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign + SWC_. + S_temp + (1|Days_Since_First), data = flux_data)
CO2_RE_model <- glmmTMB(CO2_RE_flux ~ Animal * treatment * Campaign + SWC_. + S_temp + (1|Days_Since_First), data = flux_data)

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

run_model(flux_data, CO2_RE_model)
