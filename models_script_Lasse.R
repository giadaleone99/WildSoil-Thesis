library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(bestNormalize)
library(effects)
library(MuMIn)
library(sjPlot)
library(splines)
library(effectsize)

allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
par_data <- read.csv("data/PAR_file.csv", sep = ",")
start_data <- read.csv("data/Fieldwork_data_final.csv") %>% 
  dplyr::select(UniqueID, Start_timehhmmss)

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
  dplyr::select(plotID, bulk_density)
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
  #test1 <- emmeans(model, ~ treatment|Animal|Days_Since_First)
  test2 <- emmeans(model, ~ Animal|treatment|Days_Since_First)#|SWC_.)
  #test3 <- emmeans(model, ~ Days_Since_First|treatment|Animal)#|SWC_.)
  contrast(test2, method = "pairwise") %>% as.data.frame()
  eff_size(test2, sigma = sigma(model), edf = df.residual(model))
}

# old, working models we used in our thesis
CO2_PS_model <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|Days_Since_First), data = flux_data)
CO2_PS_parmodel <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Campaign * S_temp + (1|par) + (1|Days_Since_First), data = flux_data)
CO2_RE_model <- glmmTMB(CO2_RE_flux ~ Animal * treatment + Campaign + SWC_. + bulk_density + (1|Days_Since_First), data = flux_data)
CH4_model <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * Campaign * bulk_density + (1|Days_Since_First), data = flux_data) 
N2O_model <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Campaign * S_temp * SWC_. + (1|Days_Since_First), data = flux_data) 

# models with dung age
CO2_PS_d1 <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Days_Since_First + (1|par) + (1|base_code), data = flux_data) # with PAR data, works fine
CO2_RE_d1 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * Days_Since_First * S_temp + SWC_. + (1|base_code), data = flux_data) # no good
CO2_RE_d2 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # works fine
CO2_RE_d3 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * poly(Days_Since_First, 3) + poly(S_temp, 2) * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # works but has some weird warnings
CO2_RE_d4 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * poly(Days_Since_First, 3) + poly(SWC_., 2) * poly(Days_Since_First, 2) + (1|base_code), data = flux_data) # works but has some weird warnings
N2O_d1 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Days_Since_First * S_temp * SWC_. + (1|base_code), data = flux_data) # works fine
N2O_d2 <- glmmTMB(addition_N2O_flux ~ Animal * treatment * Days_Since_First + (1|base_code), family = Gamma(link="log"), data = flux_data) # Frederiks way, DHARMa hates it
CH4_d1 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment + Days_Since_First + (1|base_code), data = flux_data) # not approved by DHARMa
CH4_d2 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment + Days_Since_First + bulk_density * S_temp + SWC_. + (1|base_code), data = flux_data) # with extra variables, not approved by DHARMa
CH4_d3 <- glmmTMB(addition_CH4_flux ~ Animal * treatment * Days_Since_First + (1|base_code), family = Gamma(link="log"), data = flux_data) # Frederiks way, DHARMa hates it
CH4_d4 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment + poly(Days_Since_First, 3) + poly(SWC_., 2) * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # works (but has no interaction effects of days since first)
CH4_d5 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * poly(Days_Since_First, 3) + poly(SWC_., 2) * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # no good
CH4_d6 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * poly(Days_Since_First, 2) + poly(SWC_., 2) * poly(Days_Since_First, 2) + (1|base_code), data = flux_data) # significant outliers, fine otherwise

# final models we are using in the paper
CO2_PS_d1 <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Days_Since_First + par + (1|base_code), data = flux_data) # with PAR data, works fine
CO2_PS_d0 <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Days_Since_First + ns(par, df = 1) + (1|base_code), data = flux_data) # with PAR data as covariate instead of random effect
CO2_RE_d2 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * poly(Days_Since_First, 3) + (1|base_code), data = flux_data) # works fine
N2O_d1 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Days_Since_First * S_temp * SWC_. + (1|base_code), data = flux_data) # works fine
CH4_d6 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * poly(Days_Since_First, 2) + poly(SWC_., 2) * poly(Days_Since_First, 2) + (1|base_code), data = flux_data) # significant outliers, fine otherwise

run_model(flux_data, N2O_d1)

standardize_parameters(CO2_PS_d0, method = "refit")
standardize_parameters(CO2_RE_d2, method = "refit")
standardize_parameters(N2O_d1, method = "refit")
standardize_parameters(CH4_d6, method = "refit")

model.sel(CO2_PS_d1, CO2_PS_d0)

CO2_PS_d1$sdr$pdHess

#confint(CH4_d6, level = 0.95)

testOutliers(CH4_d6)
res <- simulateResiduals(CO2_PS_d0)
plot(res)
testOutliers(res)
testUniformity(res)
testDispersion(res)
testZeroInflation(res)

CO2_PS_effects <- allEffects(CO2_PS_d0)
CO2_RE_effects <- allEffects(CO2_RE_d1)
N2O_effects <- allEffects(N2O_d1)
CH4_effects <- allEffects(CH4_d6)
plot(CO2_PS_effects)
plot(Effect(focal.predictors = c('treatment', 'SWC_.', 'Days_Since_First'), mod = N2O_d1))
N2O_d1 |>
  ggeffects::ggeffect() |>
  lapply(plot) |>
  patchwork::wrap_plots()

#-------------------------------model plots--------------------------------
#model plot GPP FOR PAPER

#modplot <- get_model_data(CO2_PS_d1, type = "pred")
#ggplot(bind_rows(modplot, .id="data_frame"),
#       aes(Animal, treatment, Days_Since_First, colour=modplot)) +
#  geom_line()

pred_GPPdata <- get_model_data(CO2_PS_d1,
                            type = "pred",
                            terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

pred_GPPdata <- as.data.frame(pred_GPPdata)

predictedGPPcow <- ggplot(pred_GPPdata[pred_GPPdata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = CO2_PS_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), linewidth = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("Corrected GPP ("* mu * "mol CO"[2] * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_y_continuous(limits = c(-38, 8)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Cattle control", "Cattle dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Cattle control", "Cattle dung")) +
  theme(legend.position = "bottom")
predictedGPPcow

predictedGPPhorse <- ggplot(pred_GPPdata[pred_GPPdata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = CO2_PS_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("Corrected GPP ("* mu * "mol CO"[2] * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  ylab(element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_y_continuous(limits = c(-38, 8)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedGPPhorse

predictedGPP <- predictedGPPcow + predictedGPPhorse
predictedGPP
ggsave(filename = "plots/predictedGPPpoints_paper.jpeg", plot = predictedGPP, width = 8, height = 4)

#model plot Reco FOR PAPER
pred_Recodata <- get_model_data(CO2_RE_d2,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

pred_Recodata <- as.data.frame(pred_Recodata)

predictedRecocow <- ggplot(pred_Recodata[pred_Recodata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = CO2_RE_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("R"[eco] * " (" * mu * "mol CO"[2] * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_y_continuous(limits = c(3, 34)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Cattle control", "Cattle dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Cattle control", "Cattle dung")) +
  theme(legend.position = "bottom")
predictedRecocow

predictedRecohorse <- ggplot(pred_Recodata[pred_Recodata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = CO2_RE_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("R"[eco] * " (" * mu * "mol CO"[2] * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  ylab(element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_y_continuous(limits = c(3, 34)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedRecohorse

predictedReco <- predictedRecocow + predictedRecohorse
predictedReco
ggsave(filename = "plots/predictedRecopoints_paper.jpeg", plot = predictedReco, width = 8, height = 4)

#model plot CH4 FOR PAPER
pred_CH4data <- get_model_data(CH4_d6,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))
#pred_CH4data <- as.data.frame(pred_CH4data)
normCH4 <- bestNormalize(flux_data$CH4_flux)
list1_CH4 <- dplyr::pull(pred_CH4data, predicted)
list2_CH4 <- dplyr::pull(pred_CH4data, std.error)
list3_CH4 <- dplyr::pull(pred_CH4data, conf.low)
list4_CH4 <- dplyr::pull(pred_CH4data, conf.high)
newpred1_CH4data <- predict(normCH4, newdata = list1_CH4, inverse = TRUE)
newpred2_CH4data <- predict(normCH4, newdata = list2_CH4, inverse = TRUE)
newpred3_CH4data <- predict(normCH4, newdata = list3_CH4, inverse = TRUE)
newpred4_CH4data <- predict(normCH4, newdata = list4_CH4, inverse = TRUE)
newpred_CH4data <- pred_CH4data[, -c(2, 3, 4, 5)] %>% 
  cbind(predicted=newpred1_CH4data) %>% 
  cbind(std.error=newpred2_CH4data) %>% 
  cbind(conf.low=newpred3_CH4data) %>% 
  cbind(conf.high=newpred4_CH4data)
newpred_CH4data <- newpred_CH4data[, c(1, 5, 6, 7, 8, 2, 3, 4)]

predictedCH4cow <- ggplot(newpred_CH4data[newpred_CH4data$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = CH4_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("CH"[4] * " (nmol" * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  #scale_y_continuous(limits = c(-1, 50)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Cattle control", "Cattle dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Cattle control", "Cattle dung")) +
  theme(legend.position = "bottom")
predictedCH4cow

predictedCH4horse <- ggplot(newpred_CH4data[newpred_CH4data$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = CH4_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("CH"[4] * " (nmol" * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  ylab(element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank()) +
  #scale_y_continuous(limits = c(-3, 3)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedCH4horse

predictedCH4 <- predictedCH4cow + predictedCH4horse
predictedCH4
ggsave(filename = "plots/newpredictedCH4points_paper.jpeg", plot = predictedCH4, width = 8, height = 4)

zoomCH4cow <- ggplot(newpred_CH4data[newpred_CH4data$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = CH4_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("CH"[4] * " (nmol" * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 14)) +
  coord_cartesian(ylim = c(-2, 2), xlim = c(0, 47)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Cattle control", "Cattle dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Cattle control", "Cattle dung")) +
  theme(legend.position = "bottom")
zoomCH4cow
ggsave(filename = "plots/zoomcowCH4points_paper.jpeg", plot = zoomCH4cow, width = 6, height = 3)

zoomCH4horse <- ggplot(newpred_CH4data[newpred_CH4data$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = CH4_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("CH"[4] * " (nmol" * " m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  ylab(element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 14)) +
  coord_cartesian(ylim = c(-2, 2), xlim = c(0, 47)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
zoomCH4horse
ggsave(filename = "plots/zoomhorseCH4points_paper.jpeg", plot = zoomCH4horse, width = 6, height = 3)

#model plot N2O FOR PAPER
pred_N2Odata <- get_model_data(N2O_d1,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))
#pred_N2Odata <- as.data.frame(pred_N2Odata)
normN2O <- bestNormalize(flux_data$N2O_flux)
list1_N2O <- dplyr::pull(pred_N2Odata, predicted)
list2_N2O <- dplyr::pull(pred_N2Odata, std.error)
list3_N2O <- dplyr::pull(pred_N2Odata, conf.low)
list4_N2O <- dplyr::pull(pred_N2Odata, conf.high)
newpred1_N2Odata <- predict(normN2O, newdata = list1_N2O, inverse = TRUE)
newpred2_N2Odata <- predict(normN2O, newdata = list2_N2O, inverse = TRUE)
newpred3_N2Odata <- predict(normN2O, newdata = list3_N2O, inverse = TRUE)
newpred4_N2Odata <- predict(normN2O, newdata = list4_N2O, inverse = TRUE)
newpred_N2Odata <- pred_N2Odata[, -c(2, 3, 4, 5)] %>% 
  cbind(predicted=newpred1_N2Odata) %>% 
  cbind(std.error=newpred2_N2Odata) %>% 
  cbind(conf.low=newpred3_N2Odata) %>% 
  cbind(conf.high=newpred4_N2Odata)
newpred_N2Odata <- newpred_N2Odata[, c(1, 5, 6, 7, 8, 2, 3, 4)]

predictedN2Ocow <- ggplot(newpred_N2Odata[newpred_N2Odata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = N2O_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("N"[2] * "O (nmol m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_y_continuous(limits = c(-2, 34)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Cattle control", "Cattle dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Cattle control", "Cattle dung")) +
  theme(legend.position = "bottom")
predictedN2Ocow

predictedN2Ohorse <- ggplot(newpred_N2Odata[newpred_N2Odata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = N2O_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("N"[2] * "O (nmol m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  ylab(element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_y_continuous(limits = c(-2, 34)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
predictedN2Ohorse

predictedN2O <- predictedN2Ocow + predictedN2Ohorse
predictedN2O
ggsave(filename = "plots/newpredictedN2Opoints_paper.jpeg", plot = predictedN2O, width = 8, height = 4)

zoomN2Ocow <- ggplot(newpred_N2Odata[newpred_N2Odata$facet == "Cow", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Cow", ], aes(x = Days_Since_First, y = N2O_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("N"[2] * "O (nmol m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 14)) +
  coord_cartesian(ylim = c(-1, 2.5), xlim = c(0, 47)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Cattle control", "Cattle dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Cattle control", "Cattle dung")) +
  theme(legend.position = "bottom")
zoomN2Ocow
ggsave(filename = "plots/zoomcowN2Opoints_paper.jpeg", plot = zoomN2Ocow, width = 6, height = 3)

zoomN2Ohorse <- ggplot(newpred_N2Odata[newpred_N2Odata$facet == "Horse", ], aes(x = x, y = predicted)) +
  geom_point(data = flux_data[flux_data$Animal == "Horse", ], aes(x = Days_Since_First, y = N2O_flux, color = interaction(treatment, Animal))) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill = interaction(group, facet)), alpha=0.2) +
  geom_line(aes(x = x, y = predicted, color = interaction(group, facet)), size = 1.4) +
  facet_grid(. ~facet) +
  labs(x = "Dung age (days)", y = expression("N"[2] * "O (nmol m"^{-2} * " s"^{-1} * ")"), colour = "Plot type", fill = "Plot type") +
  ylab(element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 14)) +
  coord_cartesian(ylim = c(-1, 2.5), xlim = c(0, 47)) +
  scale_color_manual(values = c("C.Cow" = "#6699CC",
                                "F.Cow" = "#004488",
                                "C.Horse" = "#EECC66",
                                "F.Horse" = "#997700"),
                     labels = c("Horse control", "Horse dung")) +
  scale_fill_manual(values = c("C.Cow" = "#6699CC",
                               "F.Cow" = "#004488",
                               "C.Horse" = "#EECC66",
                               "F.Horse" = "#997700"),
                    labels = c("Horse control", "Horse dung")) +
  theme(legend.position = "bottom")
zoomN2Ohorse
ggsave(filename = "plots/zoomhorseN2Opoints_paper.jpeg", plot = zoomN2Ohorse, width = 6, height = 3)
