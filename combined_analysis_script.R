library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)

# Import data
allflux_data <- readRDS("flux_data/clean_flux_data.rds")

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
            