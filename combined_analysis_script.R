library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(interactions)
library(plotrix)
library(ggpmisc)
library(bestNormalize)
library(effects)
library(sjPlot)


# Import data
allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
par_data <- read.csv("data/PAR_file.csv", sep = ";")
start_data <- read.csv("data/Fieldwork_data_final.csv") %>% 
  select(UniqueID, Start_timehhmmss)
veg_growth_data <- readRDS("plant_data/veg_growth_data.rds") %>% 
  dplyr::rename("plotID"= "plot_id") %>% 
  dplyr::select(-2:-9)
dung_lab <- read.csv("data/Plant_Dung_CN_pH_elements.csv", sep = ";", check.names = FALSE) %>% 
  filter(Kode_3 %in% "Dung") %>% 
  rename("base_code" = "Kode_2",
         "sample_type" = "Kode_3",
         "plotID" = "Kode_1") %>% 
  dplyr::select(-4,-5,-9)
veg_combined_data <- readRDS("plant_data/veg_combined_data.rds")%>% 
  dplyr::rename("plotID"= "plot_id",
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

# PAR stuff
CO2_PS_data <- CO2_PS_data %>% 
  left_join(start_data, by = "UniqueID") 
CO2_PS_data <- CO2_PS_data %>% 
  mutate(minute = str_sub(Start_timehhmmss, 1, 5),
         par = 0)

par_data <- par_data %>% 
  mutate(date = str_sub(min.time, 1, 10),
         min = str_sub(min.time, 12, 16))
par_data <- par_data %>% 
  mutate(longdate = as.Date(par_data$date, format = "%d-%m-%Y"))

for (i in 1:ncol(par_data)) {
  for (j in 1:ncol(CO2_PS_data)) {
    if (par_data$longdate[i] == CO2_PS_data$longdate[j] && par_data$min[i] == CO2_PS_data$minute[j]) {
      CO2_PS_data$par[j] <- mean(par_data$PAR_umol_m2_s[i], par_data$PAR_umol_m2_s[i+1], par_data$PAR_umol_m2_s[i+2], par_data$PAR_umol_m2_s[i+3], par_data$PAR_umol_m2_s[i+4])
    }
  }
} 

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
  left_join(veg_combined_data, by = "plotID")
plot_flux_data <- plot_flux_data %>% 
  full_join(soil_data, by = "plotID")

soil_data_test <- soil_data %>% 
  left_join(plot_flux_data, by = "plotID", suffix = c("_soil", "_veg")) %>% 
  mutate(Animal = case_when(
    grepl("^C", plotID) ~ "Cow",
    grepl("^H", plotID) ~ "Horse"))

growth_model_data <- plot_flux_data %>% 
  filter(!is.na(height_value))

# add bulk density to gas data
bulkdensity <- soil_data %>% 
  select(plotID, bulk_density)

flux_data <- flux_data %>% 
  left_join(bulkdensity, by = "plotID")

# create data transformed columns for CH4 and N2O fluxes
flux_data <- flux_data %>% 
  mutate(ranked_CH4_flux = rank(CH4_flux),
         normalized_CH4_flux = bestNormalize(CH4_flux)$x.t) %>% 
  mutate(ranked_N2O_flux = rank(N2O_flux),
         normalized_N2O_flux = bestNormalize(N2O_flux)$x.t)
hist(flux_data$normalized_N2O_flux)
normalizeCH4 <- orderNorm(flux_data$N2O_flux)$x.t
hist(normalizeCH4)
#create df with dung soil data
dung_soil_data <- soil_data_raw %>% 
  filter(sample_type %in% "Dung soil") %>% 
  mutate(dung_name =  case_when(
    grepl("^CG", base_code) ~ "CG",
    grepl("^HG", base_code) ~ "HG",
    base_code %in% c("CD1", "CD2") ~ "COD",
    base_code %in% c("CD3", "CD4") ~ "CND",
    base_code %in% c("HD1", "HD2") ~ "HOD",
    base_code %in% c("HD3", "HD4", "HD5") ~ "HND"))

summary_dung_soil <- dung_soil_data %>% 
  group_by(dung_name) %>% 
  summarise(
    pH_mean = mean(pH, na.rm = TRUE),
    pH_se = std.error(pH),
    PO4.P_mean = mean(PO4.P, na.rm = TRUE),
    PO4.P_se = std.error(PO4.P),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio)) %>% 
  ungroup() %>% 
  dplyr::rename("base_code" = "dung_name")
  

summary_dung_soil <- summary_dung_soil %>% 
  full_join(dung_lab, by = "base_code")


dungsoil_dung_data <- bind_rows(dung_soil_data, dung_lab) %>% 
  dplyr::select(-4,-8,-14:-18, -20:-38) %>% 
  mutate(Animal = case_when(
    grepl("^C", base_code) ~ "Cow",
    grepl("^H", base_code) ~ "Horse")) 

# regression plot of bulk density and %C
soil_data_nooutlier <- soil_data_raw %>% 
  filter(bulk_density < 4)

C_bulkdensity <- ggplot(soil_data_nooutlier, aes(bulk_density, TC)) +
  geom_point() +
  stat_poly_line(colour = "black") +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_minimal() +
  labs(x = "Bulk density (g/cm3)", y = "Total carbon")
C_bulkdensity
ggsave(C_bulkdensity, file = "plots/C_bulkdensity.jpeg",  width = 6, height = 4)

# comparing plant and soil CN
CN_horse_dung <- lm(plant_CN ~ CN_ratio_soil, data = horsedata, subset = treatment == "Fresh")
CN_horse_control <- lm(plant_CN ~ CN_ratio_soil, data = horsedata, subset = treatment == "Control")
CN_cow_dung <- lm(plant_CN ~ CN_ratio_soil, data = cowdata, subset = treatment == "Fresh")
CN_cow_control <- lm(plant_CN ~ CN_ratio_soil, data = cowdata, subset = treatment == "Control")
plot(CN_cow_dung)
ggplot(horsedata, aes(x = CN_ratio_soil, y = plant_CN, colour = treatment)) +
  geom_point() + 
  stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  geom_smooth(method = "lm") +
  labs(
    x = "CN soil",        
    y = "CN veg",        
    color = "treatment") +
  theme_minimal() +               
  scale_color_manual(values = c("Control" = "#656D4A",
                                "Fresh" = "#7F4F24"),
                     labels = c("Control", "Fresh")) +
  theme(
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CO2 RE flux by dung area and plot type")

soilCNvegCN <- ggplot(soil_data_test, aes(CN_ratio_soil, plant_CN)) +
  geom_point(aes(color = interaction(treatment, Animal)), size = 2) +
  # stat_poly_line(color = "black") +
  theme_minimal() +
  stat_poly_eq(use_label(c( "R2", "p"))) +
  geom_smooth(method = "lm", aes(group = factor(treatment), color = interaction(treatment, Animal))) +
  # stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  scale_y_continuous(expand = c(0, 5)) +
  labs(x = "Soil CN ratio", y = "Vegetation CN ratio",
       colour = "Plot type")+
  facet_wrap(~Animal) +
  scale_colour_manual(values = c("Control.Cow" = "#a4ac86",
                                 "Fresh.Cow" = "#656d4a",
                                 "Control.Horse" = "#a68a64",
                                 "Fresh.Horse" = "#7f4f24"),
                      labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))
soilCNvegCN
ggsave(soilCNvegCN, file = "plots/soilCNvegCN.jpeg",  width = 6, height = 4)

# Creating model with soil CN and plant CN
mod_soilveg_CN <- glmmTMB(plant_CN ~ CN_ratio_soil  * treatment * Animal , data = soil_data_test)

horsedata <- soil_data_test %>% 
  filter(Animal == "Horse")


cowdata <- soil_data_test %>% 
  filter(Animal == "Cow")

modHorse <- glmmTMB(plant_CN ~ CN_ratio_soil  * treatment, data = horsedata)
modCow <- glmmTMB(plant_CN ~ CN_ratio_soil  * treatment, data = cowdata)

#soil CN and total soil N model
mod_soilCN_soilN <- glmmTMB(CN_ratio_soil ~ TN_soil + treatment, data = soil_data_test)

# Model validation
simulationOutput <- simulateResiduals(fittedModel = mod_soilCN_soilN, n = 1000)
testDispersion(simulationOutput)
print(Anova(mod_soilCN_soilN))
plot(allEffects(mod_soilCN_soilN))
plot(simulationOutput)
plotResiduals(simulationOutput, form = mod_soilCN_soilN$CN_ratio_soil)
plotResiduals(simulationOutput, form = mod_soilCN_soilN$TN_soil)
plotResiduals(simulationOutput, form = mod_soilCN_soilN$Animal)
plotResiduals(simulationOutput, form = mod_soilCN_soilN$treatment)

# Post-hoc-test
test <- emmeans(mod_soilCN_soilN, ~ Animal|TN_soil|treatment)
contrast(test, method = "pairwise") %>% as.data.frame()

# Making a plot
soilvegCN_plot <- plot_model(mod_soilveg_CN, type = "pred", 
                                   terms = c("CN_ratio_soil", "treatment", "Animal")) +
  theme_minimal() +
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray"), labels = c("Control", "Dung")) +
scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  xlab("Soil CN ratio") +
  ylab("Vegetation CN ratio") +
  labs(colour = "Treatment", title = "Predicted values of plant CN ratio") +
  geom_point(aes(x = CN_ratio_soil, y = plant_CN, color = treatment),
             data = soil_data_test,
             inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(10, 18, by = 1)) +
  scale_y_continuous(breaks = seq(10, 50, by = 5)) +
  theme(panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

HorsevegCN_plot <- plot_model(modHorse, type = "pred", 
                             terms = c("CN_ratio_soil", "treatment")) +
  theme_minimal() +
  scale_color_manual(values = c("Fresh" = "#7f4f24", "Control" = "#a68a64"), labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  xlab("Soil C:N ratio") +
  ylab(element_blank()) +
  labs(colour = "Treatment", title = element_blank()) +
  geom_point(aes(x = CN_ratio_soil, y = plant_CN, color = treatment),
             data = horsedata,
             inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(10, 15, by = 1)) +
  scale_y_continuous(limits = c(10, 40)) +
  theme(panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

HorsevegCN_plot

CowvegCN_plot <- plot_model(modCow, type = "pred", 
                              terms = c("CN_ratio_soil", "treatment")) +
  theme_minimal() +
  scale_color_manual(values = c("Fresh" = "#656d4a", "Control" = "#a4ac86"), labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  xlab("Soil C:N ratio") +
  ylab("Vegetation C:N ratio") +
  labs(colour = "Treatment", title = element_blank()) +
  geom_point(aes(x = CN_ratio_soil, y = plant_CN, color = treatment),
             data = cowdata,
             inherit.aes = FALSE) +
  scale_x_continuous() +
  scale_y_continuous(limits = c(10, 40)) +
  theme(panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none")

CowvegCN_plot

combinedCNplot <- CowvegCN_plot + HorsevegCN_plot
combinedCNplot

soilvegCN_plot
ggsave(combinedCNplot, file = "veg_plots/soilvegCN_cowhorse_plot.jpeg",  width = 6, height = 4)

#Comparing pH in the dung and in the soil beneath the dung
ggplot(dungsoil_dung_data, aes(x = plotID, y = pH, fill = interaction(sample_type, Animal))) +
  geom_col() +  
  labs(
    x = "Plot ID",        
    y = "pH",        
    color = "Plot type",
    fill = "Sample type") +
  theme_minimal() +               
  scale_fill_manual(values = c("Dung.Cow" = "#656D4A",
                               "Dung soil.Cow"= "#333D29", 
                               "Dung.Horse" = "#7F4F24",
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung", "Soil below cow dung", "Horse dung", "Soil below horse dung"))+
  theme(
    panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
    panel.border = element_blank(), axis.line = element_line(),
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = c(0, 0))+
  ggtitle("Dung and Dung Soil pH")+
  facet_wrap(~ Animal, scales = "free_x")

#Comparing CN ratio in dungs and in soils beneath the dung
ggplot(dungsoil_dung_data, aes(x = plotID, y = CN_ratio, fill = interaction(sample_type, Animal))) +
  geom_col() +  
  labs(
    x = "Plot ID",        
    y = "CN ratio",        
    color = "Plot type",
    fill = "Sample type") +
  theme_minimal() +               
  scale_fill_manual(values = c("Dung.Cow" = "#656D4A",
                               "Dung soil.Cow"= "#333D29", 
                               "Dung.Horse" = "#7F4F24",
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung", "Soil below cow dung", "Horse dung", "Soil below horse dung"))+
  theme(
    panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
    panel.border = element_blank(), axis.line = element_line(),
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = c(0, 0))+
  ggtitle("Dung and Dung Soil CN ratio")+
  facet_wrap(~ Animal, scales = "free_x")


# Modeling
model <- glmmTMB(P ~ PO4.P, data = growth_model_data)
summary(mod)
print(Anova(mod))

model <- glmmTMB(height_value ~ Animal * CN_ratio * PO4.P * treatment, data = growth_model_data)
summary(model)

  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  testDispersion(simuOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = growth_model_data$Animal)
  plotResiduals(simuOutput, form = growth_model_data$treatment)
  plotResiduals(simuOutput, form = growth_model_data$CN_ratio)
  plotResiduals(simuOutput, form = growth_model_data$PO4.P)
  test1 <- emmeans(model, ~ Animal|treatment|PO4.P|CN_ratio)
  test2 <- emmeans(model, ~ treatment|PO4.P|CN_ratio|Animal)
  contrast(test2, method = "pairwise") %>% as.data.frame()
  interactions::interact_plot(model, pred = PO4.P, modx = treatment)
  modeleffects <- allEffects(model)
  plot(modeleffects)
  
# Run gas flux models with soil moisture and temperature as random factors
  run_model <- function(dataset, model) {
    #print(summary(model))
    print(Anova(model))
    simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
    testDispersion(simuOutput)
    plot(simuOutput)
    plotResiduals(simuOutput, form = dataset$Animal)
    plotResiduals(simuOutput, form = dataset$treatment)
    plotResiduals(simuOutput, form = dataset$Campaign)
    plotResiduals(simuOutput, form = dataset$S_temp)
    plotResiduals(simuOutput, form = dataset$SWC_.)
    #plotResiduals(simuOutput, form = dataset$bulk_density)
    test1 <- emmeans(model, ~ treatment|Animal|Campaign|S_temp|SWC_.)
    test2 <- emmeans(model, ~ Animal|treatment|Campaign|S_temp|SWC_.)
    test3 <- emmeans(model, ~ Campaign|treatment|Animal|S_temp|SWC_.)
    contrast(test1, method = "pairwise") %>% as.data.frame()
  }
  
CO2_PS_model <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data)
CO2_RE_model <- glmmTMB(CO2_RE_flux ~ Animal * treatment + Campaign + SWC_. + bulk_density + (1|Days_Since_First), data = flux_data)

#CH4_model1 <- glmmTMB(ranked_CH4_flux ~ Animal * treatment * Campaign * SWC_. + (1|Days_Since_First), data = flux_data) #close
#CH4_model2 <- glmmTMB(ranked_CH4_flux ~ Animal * treatment + Campaign + SWC_. + (1|Days_Since_First), data = flux_data) #everything works except for SWC
#CH4_model3 <- glmmTMB(ranked_CH4_flux ~ Animal * treatment * Campaign * bulk_density + (1|Days_Since_First), data = flux_data) #close
CH4_model <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * Campaign * bulk_density + (1|Days_Since_First), data = flux_data) #best
#CH4_model5 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * Campaign + bulk_density + S_temp + (1|Days_Since_First), data = flux_data) #good

#N2O_model1 <- glmmTMB(ranked_N2O_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data) #good
#N2O_model2 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data) #good
N2O_model <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Campaign * S_temp * SWC_. + (1|Days_Since_First), data = flux_data) #best

# models with dung age - IMPROVED FOR PAPER
CO2_PS_d1 <- glmmTMB(CO2_PS_flux ~ Animal * treatment * S_temp + (1|Days_Since_First), data = flux_data)

run_model(flux_data, N2O_model)

CO2_PS_effects <- allEffects(CO2_PS_model)
CO2_RE_effects <- allEffects(CO2_RE_model)
CH4_effects <- allEffects(CH4_model)
N2O_effects <- allEffects(N2O_model)
plot(N2O_effects)
interactions::interact_plot(N2O_model, pred = Animal, modx = Campaign)
plot(ggeffect(N2O_model))
N2O_model |>
  ggeffects::ggeffect() |>
  lapply(plot) |>
  patchwork::wrap_plots()
# Relate fluxes to dung dimensions
CO2_reg <- lm(CO2_RE_flux ~ dung_area_cm2,  data = flux_data, subset = treatment == "F")
CH4_reg <- lm(CH4_flux ~ dung_area_cm2,  data = flux_data, subset = treatment == "F")
N2O_reg <- lm(N2O_flux ~ dung_area_cm2,  data = flux_data, subset = treatment == "F")
plot(N2O_reg)
ggplot(flux_data %>% filter(treatment == "F"),
        aes(x = dung_area_cm2, y = CH4_flux, colour = interaction(treatment, Animal))) +
  geom_point() + 
  stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  geom_smooth(method = "lm") +
  labs(
    x = "Dung Area (cm²)",        
    y = "CH4 flux",        
    color = "Plot type"           
  ) +
  theme_minimal() +               
  scale_color_manual(values = c("F.Cow" = "#656D4A",
                               "F.Horse" = "#7F4F24"),
                     labels = c("Cow dung", "Horse dung")) +
  theme(
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CH4 flux by dung area and plot type")


ggplot(flux_data %>% filter(treatment == "F"),
       aes(x = dung_area_cm2, y = N2O_flux, colour = interaction(treatment, Animal))) +
  geom_point() + 
  stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  geom_smooth(method = "lm") +
  labs(
    x = "Dung Area (cm²)",        
    y = "N2O flux",        
    color = "Plot type") +
  theme_minimal() +               
  scale_color_manual(values = c("F.Cow" = "#656D4A",
                                "F.Horse" = "#7F4F24"),
                     labels = c("Cow dung", "Horse dung")) +
  theme(
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("N2O flux by dung area and plot type")

ggplot(flux_data %>% filter(treatment == "F"),
       aes(x = dung_area_cm2, y = CO2_RE_flux, colour = interaction(treatment, Animal))) +
  geom_point() + 
  stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  geom_smooth(method = "lm") +
  labs(
    x = "Dung Area (cm²)",        
    y = "CO2 RE flux",        
    color = "Plot type") +
  theme_minimal() +               
  scale_color_manual(values = c("F.Cow" = "#656D4A",
                                "F.Horse" = "#7F4F24"),
                     labels = c("Cow dung", "Horse dung")) +
  theme(
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CO2 RE flux by dung area and plot type")

# Modelling the fluxes in relation to dung size
dung_flux_data <- flux_data %>% filter(treatment == "F")

N2O_dungsize_model <- glmmTMB(N2O_flux ~ dung_area_cm2 * Animal, dung_flux_data)
CH4_dungsize_model <- glmmTMB(CH4_flux ~ dung_area_cm2 * Animal, dung_flux_data)
CO2_RE_dungsize_model <- glmmTMB(CO2_RE_flux ~ dung_area_cm2 * Animal, dung_flux_data)

print(Anova(CO2_RE_dungsize_model))
simuOutput <- simulateResiduals(fittedModel = CO2_RE_dungsize_model, n = 1000)
testDispersion(simuOutput)
plot(simuOutput)
plotResiduals(simuOutput, form = dung_flux_data$Animal)
plotResiduals(simuOutput, form = dung_flux_data$dung_area_cm2)
test <- emmeans(CO2_RE_dungsize_model, ~ Animal|dung_area_cm2)
contrast(test, method = "pairwise") %>% as.data.frame()

# Making N2O  plot
N2O_dungsize_plot <- plot_model(N2O_dungsize_model, type = "pred", 
                             terms = c("dung_area_cm2", "Animal")) +
  theme_minimal() +
  scale_color_manual(values = c("Cow" = "#656d4a", "Horse" = "#7f4f24"), labels = c("Cow", "Horse")) +
  scale_fill_manual(values = c("Horse" = "#696969", "Cow" = "#696969")) +
  xlab(expression("Dung area (cm"^2*")")) +
  ylab(expression("nmol N"[2] * "O m"^{-2} * " s"^{-1})) +
  labs(title = "Predicted values of N2O flux") +
  scale_y_continuous(breaks = seq(-5,10, by = 5)) + 
  scale_x_continuous(breaks = seq(0, 1500, by = 250)) +
  theme( 
  panel.grid.major.x = element_blank(),  
  axis.line = element_line(color = "black")) +
  geom_point(aes(x = dung_area_cm2, y = N2O_flux, color = Animal),
                data = subset(dung_flux_data, N2O_flux <= 10),
                inherit.aes = FALSE)


N2O_dungsize_plot
ggsave(N2O_dungsize_plot, file = "plots/N2O_dungsize_plot.jpeg",  width = 6, height = 4)

# Making CH4  plot
CH4_dungsize_plot <- plot_model(CH4_dungsize_model, type = "pred", 
                                terms = c("dung_area_cm2", "Animal")) +
  theme_minimal() +
  scale_color_manual(values = c("Cow" = "#656d4a", "Horse" = "#7f4f24"), labels = c("Cow", "Horse")) +
  scale_fill_manual(values = c("Horse" = "#696969", "Cow" = "#696969")) +
  xlab(expression("Dung area (cm"^2*")")) +
  ylab(expression("nmol CH"[4] * " m"^{-2} * " s"^{-1})) +
  labs(title = "Predicted values of CH4 flux") +
  scale_y_continuous(breaks = seq(-5,160, by = 25)) + 
  scale_x_continuous(breaks = seq(0, 1500, by = 250)) +
  theme(  
        panel.grid.major.x = element_blank(),  
        axis.line = element_line(color = "black")) +
  geom_point(aes(x = dung_area_cm2, y = CH4_flux, color = Animal),
             data = dung_flux_data,
             inherit.aes = FALSE)


CH4_dungsize_plot
ggsave(CH4_dungsize_plot, file = "plots/CH4_dungsize_plot.jpeg",  width = 6, height = 4)

# Making CO2 RE  plot
CO2_RE_dungsize_plot <- plot_model(CO2_RE_dungsize_model, type = "pred", 
                                terms = c("dung_area_cm2", "Animal")) +
  theme_minimal() +
  scale_color_manual(values = c("Cow" = "#656d4a", "Horse" = "#7f4f24"), labels = c("Cow", "Horse")) +
  scale_fill_manual(values = c("Horse" = "#696969", "Cow" = "#696969")) +
  xlab(expression("Dung area (cm"^2*")")) +
  ylab(expression(mu * "mol CO"[2] * " m"^{-2} * " s"^{-1})) +
  labs(title = "Predicted values of CO2 RE flux") +
  scale_y_continuous(breaks = seq(0,40, by = 5)) + 
  scale_x_continuous(breaks = seq(0, 1500, by = 250)) +
  theme( 
        panel.grid.major.x = element_blank(),  
        axis.line = element_line(color = "black")) +
  geom_point(aes(x = dung_area_cm2, y = CO2_RE_flux, color = Animal),
             data = dung_flux_data,
             inherit.aes = FALSE)


CO2_RE_dungsize_plot
ggsave(CO2_RE_dungsize_plot, file = "plots/CO2_RE_dungsize_plot.jpeg",  width = 6, height = 4)

