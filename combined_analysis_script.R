library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)

# Import data
flux_data <- readRDS("flux_data/clean_flux_data.rds")

# Combine data
CO2_PS_data <- flux_data %>% 
  filter(gastype == "CO2" && NEERE == "NEE")
CO2_RE_data <- flux_data %>% 
  filter(gastype == "CO2" && NEERE == "RE")
CH4_data <- flux_data %>% 
  filter(gastype == "CH4")
N2O_data <- flux_data %>% 
  filter(gastype == "N2O")