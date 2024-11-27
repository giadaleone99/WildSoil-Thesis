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


# Import data
allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
veg_growth_data <- readRDS("plant_data/veg_growth_data.rds") %>% 
  rename("plotID"= "plot_id") %>% 
  dplyr::select(-2:-9)
dung_lab <- read.csv("data/Plant_Dung_CN_pH_elements.csv", sep = ";", check.names = FALSE) %>% 
  filter(Kode_3 %in% "Dung") %>% 
  rename("base_code" = "Kode_2",
         "sample_type" = "Kode_3",
         "plotID" = "Kode_1") %>% 
  dplyr::select(-4,-5,-9)
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
  left_join(veg_combined_data, by = "plotID")
plot_flux_data <- plot_flux_data %>% 
  full_join(soil_data, by = "plotID")

soil_data_test <- soil_data %>% 
  left_join(plot_flux_data, by = "plotID", suffix = c("_soil", "_veg"))

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
  rename("base_code" = "dung_name")
  

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

ggsave(C_bulkdensity, file = "plots/C_bulkdensity.jpeg",  width = 6, height = 4)

# comparing plant and soil CN
soilCNvegCN <- ggplot(soil_data_test, aes(CN_ratio_soil, plant_CN)) +
  geom_point(aes(color = interaction(treatment, Animal)), size = 2) +
  stat_poly_line(color = "black") +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 5)) +
  labs(x = "Soil CN ratio", y = "Vegetation CN ratio",
       colour = "Plot type")+
  scale_colour_manual(values = c("Control.Cow" = "#a4ac86",
                                 "Fresh.Cow" = "#656d4a",
                                 "Control.Horse" = "#a68a64",
                                 "Fresh.Horse" = "#7f4f24"),
                      labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))
soilCNvegCN
ggsave(soilCNvegCN, file = "plots/soilCNvegCN.jpeg",  width = 6, height = 4)
  
# Creating model with soil CN and plant CN
mod_soilveg_CN <- glmmTMB(plant_CN ~ CN_ratio_soil  * treatment * Animal , data = soil_data_test)

# Model validation
simulationOutput <- simulateResiduals(fittedModel = mod_soilveg_CN, n = 1000)
testDispersion(simulationOutput)

plot(simulationOutput)
plotResiduals(simulationOutput, form = mod_soilveg_CN$CN_ratio_soil)
plotResiduals(simulationOutput, form = mod_soilveg_CN$plant_CN)
plotResiduals(simulationOutput, form = mod_soilveg_CN$treatment)
plotResiduals(simulationOutput, form = mod_soilveg_CN$Animal)

# Post-hoc-test
test <- emmeans(mod_soilveg_CN, ~ treatment)
contrast(test, method = "pairwise") %>% as.data.frame()

# Making a plot
soilvegCN_plot <- plot_model(mod_soilveg_CN, type = "pred", 
                                   terms = c("CN_ratio_soil", "treatment", "Animal")) +
  theme_minimal() +
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray"), labels = c("Control", "Dung")) +
scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  xlab("Soil CN ratio") +
  ylab("Vegetation CN ratio") +
  labs(colour = "Treatment") +
  scale_x_continuous(breaks = seq(10, 16, by = 1)) +
  scale_y_continuous(breaks = seq(10, 50, by = 5)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_point(aes(x = CN_ratio_soil, y = plant_CN, color = treatment),
             data = soil_data_test,
             inherit.aes = FALSE)


soilvegCN_plot


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
  run_model <- function(dataset, model) {
    #print(summary(model))
    print(Anova(model))
    simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
    testDispersion(simuOutput)
    plot(simuOutput)
    plotResiduals(simuOutput, form = dataset$Animal)
    plotResiduals(simuOutput, form = dataset$treatment)
    plotResiduals(simuOutput, form = dataset$Campaign)
    plotResiduals(simuOutput, form = dataset$bulk_density)
    test1 <- emmeans(model, ~ treatment|Animal|Campaign|bulk_density)
    test2 <- emmeans(model, ~ Animal|Campaign|treatment|bulk_density)
    test3 <- emmeans(model, ~ Campaign|treatment|Animal|bulk_density)
    contrast(test1, method = "pairwise") %>% as.data.frame()
  }
  
CO2_PS_model <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data)
CO2_RE_model <- glmmTMB(CO2_RE_flux ~ Animal * treatment + Campaign + SWC_. + bulk_density + (1|Days_Since_First), data = flux_data)

#CH4_model1 <- glmmTMB(ranked_CH4_flux ~ Animal * treatment * Campaign * SWC_. + (1|Days_Since_First), data = flux_data) #close
#CH4_model2 <- glmmTMB(ranked_CH4_flux ~ Animal * treatment + Campaign + SWC_. + (1|Days_Since_First), data = flux_data) #everything works except for SWC
#CH4_model3 <- glmmTMB(ranked_CH4_flux ~ Animal * treatment * Campaign * bulk_density + (1|Days_Since_First), data = flux_data) #close
CH4_model4 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * Campaign * bulk_density + (1|Days_Since_First), data = flux_data) #best
CH4_model5 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * Campaign + bulk_density + S_temp + (1|Days_Since_First), data = flux_data) #good

#N2O_model1 <- glmmTMB(ranked_N2O_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data) #good
N2O_model2 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data) #good
N2O_model3 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Campaign * S_temp * SWC_. + (1|Days_Since_First), data = flux_data) #best

run_model(flux_data, CH4_model4)


plot_model(CH4_model4, type = "pred", 
           terms = c("Animal", "treatment"), 
           title = "Gradient", show.p = TRUE)

dung_fluxes <- flux_data %>% filter(treatment == "F")
RE_dungarea_model <- glmmTMB(CO2_RE_flux ~ Animal * dung_area_cm2 + (1|Days_Since_First), data = dung_fluxes)
CH4_dungarea_model <- glmmTMB(CH4_flux ~ Animal * dung_area_cm2 + Campaign + (1|Days_Since_First), data = dung_fluxes)
N2O_dungarea_model <- glmmTMB(N2O_flux ~ Animal * dung_area_cm2 + Campaign + (1|Days_Since_First), data = dung_fluxes)


run_model2 <- function(dataset, model) {
  #print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  testDispersion(simuOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  #plotResiduals(simuOutput, form = dataset$Campaign)
  plotResiduals(simuOutput, form = dataset$dung_area_cm2)
  test <- emmeans(model, ~ Animal|dung_area_cm2)
  contrast(test, method = "pairwise") %>% as.data.frame()
  interactions::interact_plot(model, pred = dung_area_cm2, modx = Animal)
  
}

run_model2(dung_fluxes, RE_dungarea_model)

# Relate fluxes to dung dimensions
ggplot(flux_data %>% filter(treatment == "F"),
        aes(x = dung_area_cm2, y = CH4_flux, colour = interaction(treatment, Animal))) +
  geom_point() + 
  geom_smooth(method = "lm")+
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
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("CH4 flux by dung area and plot type")


ggplot(flux_data %>% filter(treatment == "F"),
       aes(x = dung_area_cm2, y = N2O_flux, colour = interaction(treatment, Animal))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(
    x = "Dung Area (cm²)",        
    y = "N2O flux",        
    color = "Plot type"           
  ) +
  theme_minimal() +               
  scale_color_manual(values = c("F.Cow" = "#656D4A",
                                "F.Horse" = "#7F4F24"),
                     labels = c("Cow dung", "Horse dung"))+
  theme(
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("N2O flux by dung area and plot type")

ggplot(flux_data %>% filter(treatment == "F"),
       aes(x = dung_area_cm2, y = CO2_RE_flux, colour = interaction(treatment, Animal))) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(
    x = "Dung Area (cm²)",        
    y = "CO2 RE flux",        
    color = "Plot type"           
  ) +
  theme_minimal() +               
  scale_color_manual(values = c("F.Cow" = "#656D4A",
                                "F.Horse" = "#7F4F24"),
                     labels = c("Cow dung", "Horse dung"))+
  theme(
    legend.position = "right",    
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("CO2 RE flux by dung area and plot type")

# Modelling the fluxes in relation to dung size
dung_flux_data <- flux_data %>% filter(treatment == "F")


N2O_dungsize_model <- glmmTMB(N2O_flux ~ dung_area_cm2 * Animal, dung_flux_data)

# Making a plot
N2O_dungsize_plot <- plot_model(N2O_dungsize_model, type = "pred", 
                             terms = c("dung_area_cm2", "Animal")) +
  theme_minimal() +
  scale_color_manual(values = c("Cow" = "#656d4a", "Horse" = "#7f4f24"), labels = c("Cow", "Horse")) +
  scale_fill_manual(values = c("Horse" = "#696969", "Cow" = "#696969")) +
  xlab(expression("Dung area (cm"^2*")")) +
  ylab("N2O flux (units)")
  # finish editing 
  # scale_y_continuous(limits = c(10,16)) + 
  # scale_x_continuous(limits = c(10, 40)) +
  # theme(panel.grid.minor = element_blank(),  
  #       panel.grid.major.x = element_blank(),  
  #       axis.line = element_line(color = "black")) +
  # geom_point(aes(x = plant_CN, y = CN_ratio_soil, color = treatment),
  #            data = soil_data_test,
  #            inherit.aes = FALSE)


N2O_dungsize_plot


