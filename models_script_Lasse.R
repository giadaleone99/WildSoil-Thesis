library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(bestNormalize)
library(effects)
library(MuMIn)
library(sjPlot)

allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
par_data <- read.csv("data/PAR_file.csv", sep = ",")
start_data <- read.csv("data/Fieldwork_data_final.csv") %>% 
  select(UniqueID, Start_timehhmmss)

# add the GHG data to dfs
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

# PAR stuff
CO2_PS_data <- CO2_PS_data %>% 
  left_join(start_data, by = "UniqueID") 
CO2_PS_data <- CO2_PS_data %>% 
  mutate(minute = str_sub(Start_timehhmmss, 1, 5),
         par = NA)
for (i in 1:nrow(CO2_PS_data)) {
  if (substring(CO2_PS_data$minute[i], 1, 1) == "0") {
    CO2_PS_data$minute[i] = substring(CO2_PS_data$minute[i], 2)
  }
}

# Join gas data per date
flux_data <- CO2_PS_data %>% 
  left_join(CO2_RE_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-21:-25, -27:-35)
flux_data <- flux_data %>% 
  left_join(CH4_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-22:-26, -28:-36)
flux_data <- flux_data %>% 
  left_join(N2O_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-1:-5, -7:-9, -12:-17, -23, -29)
flux_data <- flux_data[, c(9, 10, 2, 11, 12, 14, 15, 16, 3, 20, 21, 1, 7, 8, 13, 17, 18, 19, 4, 5, 6)]

# more PAR stuff
par_data <- par_data %>% 
  mutate(date = str_sub(min.time, 1, 10),
         min = str_sub(min.time, 12, 16))
par_data <- par_data %>% 
  mutate(longdate = as.Date(par_data$date, format = "%d-%m-%Y"))
for (j in 1:nrow(flux_data)) {
  for (i in 1:nrow(par_data)) {
    if (par_data$longdate[i] == flux_data$longdate[j] && par_data$min[i] == flux_data$minute[j]) {
      flux_data$par[j] <- mean(c(par_data$PAR_umol_m2_s[i], par_data$PAR_umol_m2_s[i+1], par_data$PAR_umol_m2_s[i+2], par_data$PAR_umol_m2_s[i+3], par_data$PAR_umol_m2_s[i+4]))
    }
  }
} 

# add bulk density to df
soil_data <- soil_data_raw %>% 
  filter(sample_type %in% c("Fresh", "Control"))
bulkdensity <- soil_data %>% 
  select(plotID, bulk_density)
flux_data <- flux_data %>% 
  left_join(bulkdensity, by = "plotID")

# create data transformed columns for CH4 and N2O fluxes
flux_data <- flux_data %>% 
  mutate(ranked_CH4_flux = rank(CH4_flux),
         normalized_CH4_flux = bestNormalize(CH4_flux)$x.t,
         addition_CH4_flux = CH4_flux + 2.904) %>% 
  mutate(ranked_N2O_flux = rank(N2O_flux),
         normalized_N2O_flux = bestNormalize(N2O_flux)$x.t,
         addition_N2O_flux = N2O_flux + 1.471)

#------- Models ---------
run_model <- function(dataset, model) {
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 10000)
  testDispersion(simuOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$treatment)
  plotResiduals(simuOutput, form = dataset$Days_Since_First)
  #plotResiduals(simuOutput, form = dataset$S_temp)
  #plotResiduals(simuOutput, form = dataset$SWC_.)
  #plotResiduals(simuOutput, form = dataset$bulk_density)
  test1 <- emmeans(model, ~ treatment|Animal)
  #test2 <- emmeans(model, ~ Animal|treatment|S_temp)#|SWC_.)
  #test3 <- emmeans(model, ~ S_temp|treatment|Animal)#|SWC_.)
  contrast(test1, method = "pairwise") %>% as.data.frame()
}

# old, working models we used in our thesis
CO2_PS_model <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data)
CO2_PS_parmodel <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|par) + (1|Days_Since_First), data = flux_data)
CO2_RE_model <- glmmTMB(CO2_RE_flux ~ Animal * treatment + Campaign + SWC_. + bulk_density + (1|Days_Since_First), data = flux_data)
CH4_model <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * Campaign * bulk_density + (1|Days_Since_First), data = flux_data) 
N2O_model <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Campaign * S_temp * SWC_. + (1|Days_Since_First), data = flux_data) 

# models with dung age - IMPROVED FOR PAPER
CO2_PS_d1 <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Days_Since_First + (1|par) + (1|base_code), data = flux_data) # with PAR data, works fine
CO2_RE_d1 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * Days_Since_First * S_temp + SWC_. + (1|base_code), data = flux_data) # works fine
N2O_d1 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Days_Since_First * S_temp * SWC_. + (1|base_code), data = flux_data) # works fine
N2O_d2 <- glmmTMB(addition_N2O_flux ~ Animal * treatment * Days_Since_First + (1|base_code), family = Gamma(link="log"), data = flux_data) # Frederiks way, DHARMa hates it
CH4_d1 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment + Days_Since_First + (1|base_code), data = flux_data) # not approved by DHARMa
CH4_d2 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment + Days_Since_First + bulk_density * S_temp + SWC_. + (1|base_code), data = flux_data) # with extra variables, not approved by DHARMa
CH4_d3 <- glmmTMB(addition_CH4_flux ~ Animal * treatment * Days_Since_First + (1|base_code), family = Gamma(link="log"), data = flux_data) # Frederiks way, DHARMa hates it
CH4_d4 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment + poly(Days_Since_First, 3) + poly(SWC_., 2) * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # works (but has no interaction effects of days since first)
CH4_d5 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * poly(Days_Since_First, 3) + poly(SWC_., 2) * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # no good
CH4_d6 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * poly(Days_Since_First, 2) + poly(SWC_., 2) * poly(Days_Since_First, 2) + (1|base_code), data = flux_data) # significant outliers, fine otherwise


run_model(flux_data, CH4_d6)

model.sel(CH4_d4, CH4_d5, CH4_d6)

testOutliers(CH4_d6)

#model plot GPP FOR PAPER
modplot <- get_model_data(CO2_PS_d1, type = "pred")
ggplot(bind_rows(modplot, .id="data_frame"),
       aes(Animal, treatment, Days_Since_First, colour=modplot)) +
  geom_line()

pred_GPPdata <- get_model_data(CO2_PS_d1,
                            type = "pred",
                            terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

pred_GPPdata <- as.data.frame(pred_GPPdata)

predictedGPPcow <- ggplot(pred_GPPdata[pred_GPPdata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = CO2_PS_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted "* mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Cow control", "Cow dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Cow control", "Cow dung")) +
  theme(legend.position = "bottom")
predictedGPPcow

predictedGPPhorse <- ggplot(pred_GPPdata[pred_GPPdata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = CO2_PS_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted "* mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedGPPhorse

predictedGPP <- predictedGPPcow + predictedGPPhorse
predictedGPP
ggsave(filename = "plots/predictedGPPpoints_paper.jpeg", plot = predictedGPP, width = 8, height = 4)

#model plot Reco FOR PAPER
pred_Recodata <- get_model_data(CO2_RE_d1,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

pred_Recodata <- as.data.frame(pred_Recodata)

predictedRecocow <- ggplot(pred_Recodata[pred_Recodata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = CO2_RE_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted "* mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Cow control", "Cow dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Cow control", "Cow dung")) +
  theme(legend.position = "bottom")
predictedRecocow

predictedRecohorse <- ggplot(pred_Recodata[pred_Recodata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = CO2_RE_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted "* mu * "mol CO"[2] * " m"^{-2} * " s"^{-1}), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedRecohorse

predictedReco <- predictedRecocow + predictedRecohorse
predictedReco
ggsave(filename = "plots/predictedRecopoints_paper.jpeg", plot = predictedReco, width = 8, height = 4)

#model plot N2O FOR PAPER
pred_N2Odata <- get_model_data(N2O_d1,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

pred_N2Odata <- as.data.frame(pred_N2Odata)

predictedN2Ocow <- ggplot(pred_N2Odata[pred_N2Odata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = normalized_N2O_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted nmol N"[2] * "O m"^{-2} * " s"^{-1}), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Cow control", "Cow dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Cow control", "Cow dung")) +
  theme(legend.position = "bottom")
predictedN2Ocow

predictedN2Ohorse <- ggplot(pred_N2Odata[pred_N2Odata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = normalized_N2O_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted nmol N"[2] * "O m"^{-2} * " s"^{-1}), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedN2Ohorse

predictedN2O <- predictedN2Ocow + predictedN2Ohorse
predictedN2O
ggsave(filename = "plots/predictedN2Opoints_paper.jpeg", plot = predictedN2O, width = 8, height = 4)

#model plot CH4 FOR PAPER
modplot <- get_model_data(CH4_d6, type = "pred")
ggplot(bind_rows(modplot, .id="data_frame"),
       aes(Animal, treatment, Days_Since_First, SWC_., colour=modplot)) +
  geom_line()

pred_CH4data <- get_model_data(CH4_d6,
                            type = "pred",
                            terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

pred_CH4data <- as.data.frame(pred_CH4data)

predictedCH4cow <- ggplot(pred_CH4data[pred_CH4data$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = normalized_CH4_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted CH"[4] * " flux (nmol" * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                               "F.Cow" = "#656d4a",
                               "C.Horse" = "#a68a64",
                               "F.Horse" = "#7f4f24"),
                     labels = c("Cow control", "Cow dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                                "F.Cow" = "#4b5037",
                                "C.Horse" = "#c3a276",
                                "F.Horse" = "#66401a"),
                     labels = c("Cow control", "Cow dung")) +
  theme(legend.position = "bottom")
predictedCH4cow

predictedCH4horse <- ggplot(pred_CH4data[pred_CH4data$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = normalized_CH4_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Days since first measurement", y = expression("Predicted CH"[4] * " flux (nmol" * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_minimal() +
  scale_color_manual(values = c("C.Cow" = "#a4ac86",
                                "F.Cow" = "#656d4a",
                                "C.Horse" = "#a68a64",
                                "F.Horse" = "#7f4f24"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#c4cda0",
                               "F.Cow" = "#4b5037",
                               "C.Horse" = "#c3a276",
                               "F.Horse" = "#66401a"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedCH4horse

predictedCH4 <- predictedCH4cow + predictedCH4horse
predictedCH4
ggsave(filename = "plots/predictedCH4points_paper.jpeg", plot = predictedCH4, width = 8, height = 4)
       