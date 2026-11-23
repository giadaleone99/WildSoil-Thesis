### all data analyses
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(bestNormalize)
library(effects)
library(MuMIn)
library(splines)
library(effectsize)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)
library(dplyr)
library(tidyr)

# soil analysis ---------------------------------------------------------
com_soil_data <- read.csv("clean_data/soil_data.csv")

# summary stats
summary_stats <- comb_soil_data %>% 
  group_by(Animal, sample_type) %>%
  summarise(
    N_mean = mean(TN, na.rm = TRUE),
    N_se = std.error(TN),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    pH_mean = mean(pH, na.rm = TRUE),
    pH_se = std.error(pH),
    PO4.P_mean = mean(PO4.P, na.rm = TRUE),
    PO4.P_se = std.error(PO4.P)
  )

# boxplots
comb_soil_data <- comb_soil_data %>%
  mutate(
    sample_type_animal = interaction(sample_type, Animal),
    TN = TN * 10) #this changes the value from % to g/kg
comb_soil_data$dung_age <- ordered(comb_soil_data$dung_age, levels = c("3-4", "18", "49-50"))

comb_means <- comb_soil_data %>% 
  group_by(sample_type, Animal, dung_age, sample_type_animal) %>% 
  dplyr::summarise(CN_se = std.error(CN_ratio),
                   CN_ratio=mean(CN_ratio),
                   TN_se = std.error(TN),
                   TN=mean(TN),
                   TC_se = std.error(TC),
                   TC=mean(TC),
                   pH_se = std.error(pH),
                   pH=mean(pH),
                   PO4.P_se = std.error(PO4.P),
                   PO4.P=mean(PO4.P))

# CN ratio
cnplot <- ggplot(comb_soil_data, aes(x = Animal, y = CN_ratio, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = comb_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=CN_ratio-CN_se, ymax=CN_ratio+CN_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = comb_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = sample_type, color = Animal)
  ) +
  labs(shape = "Sampling location") +
  ylab("Soil C:N ratio\n") +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(), axis.title = element_text(size = 12),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) +
  scale_shape_manual(values = c(23,24,22), labels = c("Soil below dung", "Soil beside dung", "Control")) +
  scale_fill_manual(values = c("#004488", "#004488", "#004488", "#997700", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488", "#997700"),
                     labels = c("Cattle", "Horse"))+
  expand_limits(y = c(11, 16)) +
  scale_y_continuous(breaks = seq(11, 16, by = 1))
cnplot
ggsave(filename = "soil_plots/newcn.jpeg", plot = cnplot, width = 8, height = 5)

# Total N
tnplot <- ggplot(comb_soil_data, aes(x = Animal, y = TN, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = comb_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=TN-TN_se, ymax=TN+TN_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = comb_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = sample_type, color = Animal)
  ) +
  labs(shape = "Sampling location") +
  ylab(expression("Total N (g N kg"^{-1} * " soil)")) +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(), axis.title = element_text(size = 12),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) +
  scale_shape_manual(values = c(23,24,22), labels = c("Soil below dung", "Soil beside dung", "Control")) +
  scale_fill_manual(values = c("#004488", "#004488", "#004488", "#997700", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488", "#997700"),
                     labels = c("Cattle", "Horse")) +
  expand_limits(y = c(1.5, 4)) 
tnplot
ggsave(filename = "soil_plots/newtn.jpeg", plot = tnplot, width = 8, height = 5)

# Total C
tcplot <- ggplot(comb_soil_data, aes(x = Animal, y = TC, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = comb_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=TC-TC_se, ymax=TC+TC_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = comb_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = sample_type, color = Animal)
  ) +
  labs(shape = "Sampling location") +
  ylab("Total C% in soil\n") +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(), axis.title = element_text(size = 12),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) +
  scale_shape_manual(values = c(23,24,22), labels = c("Soil below dung", "Soil beside dung", "Control")) +
  scale_fill_manual(values = c("#004488", "#004488", "#004488", "#997700", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488", "#997700"),
                     labels = c("Cattle", "Horse")) +
  expand_limits(y = c(2, 6.5)) 
tcplot
ggsave(filename = "soil_plots/newtc.jpeg", plot = tcplot, width = 8, height = 5)

# pH
phplot <- ggplot(comb_soil_data, aes(x = Animal, y = pH, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = comb_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=pH-pH_se, ymax=pH+pH_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = comb_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = sample_type, color = Animal)
  ) +
  labs(shape = "Sampling location") +
  ylab("Soil pH") +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(), axis.title = element_text(size = 12),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) +
  scale_shape_manual(values = c(23,24,22), labels = c("Soil below dung", "Soil beside dung", "Control")) +
  scale_fill_manual(values = c("#004488", "#004488", "#004488", "#997700", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488", "#997700"),
                     labels = c("Cattle", "Horse"))+
  expand_limits(y = c(4, 7)) +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5))
phplot
ggsave(filename = "soil_plots/newph.jpeg", plot = phplot, width = 8, height = 5)

# PO4.P
po4plot <- ggplot(comb_soil_data, aes(x = Animal, y = PO4.P, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = comb_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=PO4.P-PO4.P_se, ymax=PO4.P+PO4.P_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = comb_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = sample_type, color = Animal)
  ) +
  labs(shape = "Sampling location") +
  ylab(expression("Plant available P (mg P kg"^{-1} * " soil)")) +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(), axis.title = element_text(size = 12),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 10), strip.text = element_text(size = 12)) +
  scale_shape_manual(values = c(23,24,22), labels = c("Soil below dung", "Soil beside dung", "Control")) +
  scale_fill_manual(values = c("#004488", "#004488", "#004488", "#997700", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488", "#997700"),
                     labels = c("Cattle", "Horse"))+
  expand_limits(y = c(50, 370)) +
  scale_y_continuous(breaks = seq(50, 370, by = 50))
po4plot
ggsave(filename = "soil_plots/newpo4p.jpeg", plot = po4plot, width = 8, height = 5)

# modelling 
run_model <- function(dataset, model) {
  print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  testDispersion(simuOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$sample_type)
  plotResiduals(simuOutput, form = dataset$dung_age)
  test <- emmeans(model, ~ sample_type|Animal|dung_age)
  test2 <- emmeans(model, ~ dung_age|sample_type|Animal)
  test3 <- emmeans(model, ~ Animal|dung_age|sample_type)
  contrast(test, method = "pairwise") %>% as.data.frame()
}

Soil_CNratio_d1 <- glmmTMB(CN_ratio ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

Soil_TN_d1 <- glmmTMB(TN ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

Soil_pH_d1 <- glmmTMB(pH ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

Soil_PO4.P_d1 <- glmmTMB(PO4.P ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

run_model(comb_soil_data, Soil_PO4.P_d1) # change accordingly

# model effects
Soil_CNratio_effects <- allEffects(Soil_CNratio_d1)
Soil_TN_effects <- allEffects(Soil_TN_d1)
Soil_pH_effects <- allEffects(Soil_pH_d1)
Soil_PO4.P_effects <- allEffects(Soil_PO4.P_d1)
print(Soil_pH_effects)
plot(Soil_pH_effects)

# vegetation analysis ---------------------------------------------------

veg_species <- read.csv("clean_data/veg_species.csv")
veg_combined <- read.csv("clean_data/veg_data.csv")
veg_growth <- read.csv("clean_data/veg_growth.csv")
presence_absence_matrix <- write.csv(presence_absence_matrix, file="clean_data/RDA_presence.csv")
species_df <- write.csv(species_df, file="clean_data/RDA_species.csv")
soil_df <- write.csv(soil_df, file="clean_data/RDA_soil.csv")


# boxplots
percentspecies <- ggplot(veg_species, aes(x = plot_id, y = dry_weight, fill = factor(veg_class)))  +
  geom_bar(position = "fill", stat = "identity") +
  xlab("\nPlot ID") + ylab("Proportion of total forb weight") +
  theme_minimal() +
  scale_fill_manual(values = c("#D95F30FF", "#A89E5EFF", "#8785B2FF", "#FE00CE", "#0DF9FF", "#F8D564FF", "#FF9616", "#479B55", "#EEA6FB" ,
                               "#DC587D", "#D626FF", "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", "#22FFA7", "#FED4C4", "#E3EE9E", "#86CE00",
                               "#BC7196", "#7E7DCD", "#FC6955", "#E48F72", "#846D86FF"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), strip.text = element_text(size=13),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), legend.position = "bottom") +
  labs(fill = "Forb species") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Proportions of forb species per age group")
percentspecies
ggsave(filename = "veg_plots/percentspeciesnoage.jpeg", plot = percentspecies, width = 10, height = 8)

veg_combined$dung_age <- ordered(veg_combined$dung_age, levels = c("3-4", "18", "49-50"))

vegcn_means <- veg_combined %>% 
  group_by(treatment, Animal, dung_age, Animal_treatment) %>% 
  summarise(CN_se = std.error(CN_ratio), CN_ratio=mean(CN_ratio))
vegcn_means$dung_age <- ordered(vegcn_means$dung_age, levels = c("3-4", "18", "49-50"))

veggrowth_means <- veg_growth %>% 
  group_by(treatment, Animal, dung_age, Animal_treatment) %>% 
  summarise(height_se = std.error(height_value), height_value=mean(height_value))
veggrowth_means$dung_age <- ordered(veggrowth_means$dung_age, levels = c("18", "49-50"))

vegtn_means <- veg_combined %>% 
  group_by(treatment, Animal, dung_age, Animal_treatment) %>% 
  summarise(TN_se = std.error(TN), TN=mean(TN))
vegcn_means$dung_age <- ordered(vegtn_means$dung_age, levels = c("3-4", "18", "49-50"))

vegcnplot <- ggplot(veg_combined, aes(x = Animal, y = CN_ratio, fill = Animal_treatment)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = vegcn_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=CN_ratio-CN_se, ymax=CN_ratio+CN_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = vegcn_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = treatment, color = Animal)
  ) +
  ylab("Vegetation C:N ratio") +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) +
  labs(shape = "Treatment") +
  scale_shape_manual(values = c(23,24,22), labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("#004488", "#004488", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488","#997700"),
                     labels = c("Cattle", "Horse")) +
  expand_limits(y = c(10, 40)) 
vegcnplot
ggsave(filename = "veg_plots/vegcn_new.jpeg", plot = vegcnplot, width = 6, height = 4)

veggrowthplot <- ggplot(veg_growth, aes(x = Animal, y = height_value, fill = Animal_treatment)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = veggrowth_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=height_value-height_se, ymax=height_value+height_se)) +
  scale_color_manual(values = c("#6699CC", "#EECC66"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = veggrowth_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = treatment, color = Animal)
  ) +
  ylab("Vegetation growth (cm)") +
  theme_minimal() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) +
  labs(shape = "Treatment") +
  scale_shape_manual(values = c(23,24,22), labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("#004488", "#004488", "#997700", "#997700")) +
  scale_color_manual("Herbivore species", values = c("#004488","#997700"),
                     labels = c("Cattle", "Horse")) 
veggrowthplot
ggsave(filename = "veg_plots/veggrowth_new.jpeg", plot = veggrowthplot, width = 6, height = 4)

vegtnplot <- ggplot(veg_combined, aes(x = Animal, y = TN, fill = Animal_treatment)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days", "18" = "18 days", "49-50" = "49-50 days")), strip.position = "bottom") +
  geom_point(position = position_dodge(width = 1), size = 1.3, aes(color = Animal)) +
  geom_errorbar(data = vegtn_means, position = position_dodge(width = 1), linewidth = 0.3, aes(ymin=TN-TN_se, ymax=TN+TN_se)) +
  scale_color_manual(values = c("#A4AC86", "#A68A64"), guide = "none") +
  new_scale_color() +
  guides(fill = "none") +
  geom_point(
    data = vegtn_means,
    size = 3,
    position = position_dodge(width = 1), aes(shape = treatment, color = Animal)
  ) +
  ylab("Total N % in soil") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) +
  labs(shape = "Treatment") +
  scale_shape_manual(values = c(23,24,22), labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("#656D4A", "#656D4A", "#7F4F24", "#7F4F24")) +
  scale_color_manual("Herbivore species", values = c("#656D4A","#7F4F24"),
                     labels = c("Cattle", "Horse"))
vegtnplot
ggsave(filename = "veg_plots/vegtn_new.jpeg", plot = vegtnplot, width = 6, height = 4)

# modelling
run_model <- function(dataset, model) {
  print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$treatment)
  plotResiduals(simuOutput, form = dataset$dung_age)
  test <- emmeans(model, ~ treatment|Animal|dung_age)
  test2 <- emmeans(model, ~ Animal|treatment|dung_age)
  test3 <- emmeans(model, ~ dung_age|Animal|treatment)
  contrast(test3, method = "pairwise") %>% as.data.frame()
  eff_size(test3, sigma = sigma(model), edf = df.residual(model))
}
veg_growth_d1 <- glmmTMB(height_value ~ Animal * treatment * dung_age + (1|base_code), data = veg_growth)
veg_cn_d1 <- glmmTMB(CN_ratio ~ Animal * treatment * dung_age + (1|base_code), data = veg_combined)
veg_tn_d1 <- glmmTMB(TN ~ Animal * treatment * dung_age + (1|base_code), data = veg_combined)

run_model(veg_combined, veg_cn_d1) #change accordingly

# model effects
veg_growth_effects <- allEffects(veg_growth_d1)
veg_cn_effects <- allEffects(veg_cn_d1)
veg_tn_effects <- allEffects(veg_tn_d1)
plot(veg_growth_effects)

# RDA with weights
species_df.hell <- decostand(species_df, "hellinger") #hellinger-transform the species dataset

rda_model <- rda(species_df.hell ~ bulk_density + pH + PO4.P + TN + dung_age + sample_type + Animal, data = soil_df)
summary_rda <- summary(rda_model)
summary_rda

plot(rda_model, scaling = 1)
plot(rda_model, scaling = 2)

(R2 <- RsquareAdj(rda_model)$r.squared)
(R2adj <- RsquareAdj(rda_model)$adj.r.squared)

anova(rda_model, permutations = how(nperm = 999))

vif.cca(rda_model) # watch for things >10

# variable selection
mod0 <- rda(species_df.hell ~ 1, data = soil_df)

step.forward <- ordistep(mod0,
                         scope = formula(rda_model),
                         direction = "forward",
                         permutations = how(nperm = 999))


RsquareAdj(step.forward)
step.forward$anova

step.forward2 <- ordiR2step(mod0,
                            scope = formula(rda_model),
                            direction = "forward",
                            R2scope = TRUE,
                            permutations = how(nperm = 999))

RsquareAdj(step.forward2)
step.forward2$anova

# Variation Partitioning
bulk <- soil_df$bulk_density
ph <- soil_df$pH
tn <- soil_df$TN
PO4.P <- soil_df$PO4.P
soil <- soil_df[c("bulk_density", "pH", "TN", "PO4.P")]
age <- soil_df$dung_age
sample_type <- soil_df$sample_type
Animal <- soil_df$Animal

# Soil
(spe.part.all <- varpart(species_df.hell, PO4.P, ph, tn, bulk))

jpeg(file="veg_plots/Soil_variation_partitioning.jpeg")
plot(spe.part.all, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("PO4.P","pH", "N", "Bulk"))
dev.off()

# Test PO4
anova(rda(species_df.hell, PO4.P, permutation = how(nperm = 9999)))

# Fraction [a], pure PO4
anova(rda(species_df.hell, PO4.P, cbind(ph, tn, bulk)), permutation = how(nperm = 9999))

# Test pH
anova(rda(species_df.hell, ph, permutation = how(nperm = 9999)))

# Fraction [b], pure pH
anova(rda(species_df.hell, ph, cbind(PO4.P, tn, bulk)), permutation = how(nperm = 9999))

# Test TN
anova(rda(species_df.hell, tn, permutation = how(nperm = 9999)))

# Fraction [c], pure TN
anova(rda(species_df.hell, tn, cbind(PO4.P, ph, bulk)), permutation = how(nperm = 9999))

# Test bulk density
anova(rda(species_df.hell, bulk, permutation = how(nperm = 9999)))

# Fraction [d], pure bulk density
anova(rda(species_df.hell, bulk, cbind(PO4.P, ph, tn)), permutation = how(nperm = 9999))

# Soil and other variables
(spe.part.all2 <- varpart(species_df.hell, soil, age, sample_type, Animal))
jpeg(file="veg_plots/variation_partitioning.jpeg")
plot(spe.part.all2, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("Soil","Dung age", "Treatment", "Animal"))
dev.off()

# Test soil
anova(rda(species_df.hell, soil, permutation = how(nperm = 9999)))

# Fraction [a], pure soil
anova(rda(species_df.hell, soil, cbind(age, sample_type, Animal)), permutation = how(nperm = 9999))

# Test dung age
anova(rda(species_df.hell, age, permutation = how(nperm = 9999)))

# Fraction [a] only dung age
anova(rda(species_df.hell, age, cbind(soil, sample_type, Animal)), permutation = how(nperm = 9999))

# Test sample_type
anova(rda(species_df.hell, sample_type, permutation = how(nperm = 9999)))

# Fraction [a] only sample_type
anova(rda(species_df.hell, sample_type, cbind(soil, age, Animal)), permutation = how(nperm = 9999))

# Test Animal
anova(rda(species_df.hell, Animal, permutation = how(nperm = 9999)))

# Fraction [a] only Animal
anova(rda(species_df.hell, Animal, cbind(soil, age, sample_type)), permutation = how(nperm = 9999))

# RDA with the presence absence of forbs and graminoids
rda_model2 <- rda(presence_absence_matrix ~ ., soil_df)
summary_rda2 <- summary(rda_model2)
summary_rda2

plot(rda_model2, scaling = 1)
plot(rda_model2, scaling = 2)

(R2 <- RsquareAdj(rda_model2)$r.squared)
(R2adj <- RsquareAdj(rda_model2)$adj.r.squared)

anova(rda_model2, permutations = how(nperm = 999))

vif.cca(rda_model2) # watch for things >10

# variable selection
mod0 <- rda(presence_absence_matrix ~ 1, data = soil_df)

step.forward <- ordistep(mod0,
                         scope = formula(rda_model2),
                         direction = "forward",
                         permutations = how(nperm = 999))


RsquareAdj(step.forward)
step.forward$anova

step.forward2 <- ordiR2step(mod0,
                            scope = formula(rda_model2),
                            direction = "forward",
                            R2scope = TRUE,
                            permutations = how(nperm = 999))

RsquareAdj(step.forward2)
step.forward2$anova

# Variation Partitioning ####
bulk <- soil_df$bulk_density
ph <- soil_df$pH
tn <- soil_df$TN
PO4.P <- soil_df$PO4.P
soil <- soil_df[c("bulk_density", "pH", "TN", "PO4.P")]
age <- soil_df$dung_age
sample_type <- soil_df$sample_type
Animal <- soil_df$Animal

# Soil
(spe.part.all3 <- varpart(presence_absence_matrix, PO4.P, ph, tn, bulk))

jpeg(file="veg_plots/Soil_variation_partitioning_absence.jpeg")
plot(spe.part.all3, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("PO4.P","pH", "N", "Bulk"))
dev.off()

# Test PO4
anova(rda(presence_absence_matrix, PO4.P, permutation = how(nperm = 9999)))

# Fraction [a], pure PO4
anova(rda(presence_absence_matrix, PO4.P, cbind(ph, tn, bulk)), permutation = how(nperm = 9999))

# Test pH
anova(rda(presence_absence_matrix, ph, permutation = how(nperm = 9999)))

# Fraction [b], pure pH
anova(rda(presence_absence_matrix, ph, cbind(PO4.P, tn, bulk)), permutation = how(nperm = 9999))

# Test TN
anova(rda(presence_absence_matrix, tn, permutation = how(nperm = 9999)))

# Fraction [c], pure TN
anova(rda(presence_absence_matrix, tn, cbind(PO4.P, ph, bulk)), permutation = how(nperm = 9999))

# Test bulk density
anova(rda(presence_absence_matrix, bulk, permutation = how(nperm = 9999)))

# Fraction [d], pure bulk density
anova(rda(presence_absence_matrix, bulk, cbind(PO4.P, ph, tn)), permutation = how(nperm = 9999))

# Soil and other variables
(spe.part.all4 <- varpart(presence_absence_matrix, soil, age, sample_type, Animal))
jpeg(file="veg_plots/variation_partitioning_absence.jpeg")
plot(spe.part.all4, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("Soil","Dung age", "Treatment", "Animal"))
dev.off()

# Test soil
anova(rda(presence_absence_matrix, soil, permutation = how(nperm = 9999)))

# Fraction [a], pure soil
anova(rda(presence_absence_matrix, soil, cbind(age, sample_type, Animal)), permutation = how(nperm = 9999))

# Test dung age
anova(rda(presence_absence_matrix, age, permutation = how(nperm = 9999)))

# Fraction [a] only dung age
anova(rda(presence_absence_matrix, age, cbind(soil, sample_type, Animal)), permutation = how(nperm = 9999))

# Test sample_type
anova(rda(presence_absence_matrix, sample_type, permutation = how(nperm = 9999)))

# Fraction [a] only sample_type
anova(rda(presence_absence_matrix, sample_type, cbind(soil, age, Animal)), permutation = how(nperm = 9999))

# Test Animal
anova(rda(presence_absence_matrix, Animal, permutation = how(nperm = 9999)))

# Fraction [a] only Animal
anova(rda(presence_absence_matrix, Animal, cbind(soil, age, sample_type)), permutation = how(nperm = 9999))

# GHG analysis ----------------------------------------------------------
flux_data <- read.csv("clean_data/flux_data.csv")

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
  test1 <- emmeans(model, ~ treatment|Animal|Days_Since_First)
  test2 <- emmeans(model, ~ Animal|treatment|Days_Since_First)#|SWC_.)
  test3 <- emmeans(model, ~ Days_Since_First|treatment|Animal)#|SWC_.)
  contrast(test2, method = "pairwise") %>% as.data.frame()
  eff_size(test2, sigma = sigma(model), edf = df.residual(model))
}

CO2_PS_d0 <- glmmTMB(CO2_PS_flux ~ Animal * treatment * Days_Since_First + ns(par, df = 1) + (1|base_code), data = flux_data)
CO2_RE_d2 <- glmmTMB(CO2_RE_flux ~ Animal * treatment * poly(Days_Since_First, 3) + (1|base_code), data = flux_data)
N2O_d1 <- glmmTMB(normalized_N2O_flux ~ Animal * treatment * Days_Since_First * S_temp * SWC_. + (1|base_code), data = flux_data)
CH4_d6 <- glmmTMB(normalized_CH4_flux ~ Animal * treatment * poly(Days_Since_First, 2) + poly(SWC_., 2) * poly(Days_Since_First, 2) + (1|base_code), data = flux_data)

run_model(flux_data, N2O_d1) # change accordingly

CO2_PS_effects <- allEffects(CO2_PS_d0)
CO2_RE_effects <- allEffects(CO2_RE_d1)
N2O_effects <- allEffects(N2O_d1)
CH4_effects <- allEffects(CH4_d6)
plot(CO2_PS_effects)

#model plots
pred_GPPdata <- get_model_data(CO2_PS_d0,
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


pred_CH4data <- get_model_data(CH4_d6,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))

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

pred_N2Odata <- get_model_data(N2O_d1,
                               type = "pred",
                               terms = c("Days_Since_First[0:47, by=0.5]", "treatment", "Animal"))
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

# combined analysis -----------------------------------------------------
comb_flux <- read.csv("clean_data/combined_flux_data.csv")
comb_soil <- read.csv("clean_data/combined_soil_data.csv")

run_model <- function(dataset, model) {
  print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  testDispersion(simuOutput)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$treatment)
  #plotResiduals(simuOutput, form = dataset$TN_soil)
  #plotResiduals(simuOutput, form = dataset$PO4.P_soil)
  test1 <- emmeans(model, ~ treatment|Animal)
  #test2 <- emmeans(model, ~ Animal|treatment|TN_soil)
  #test3 <- emmeans(model, ~ PO4.P_soil|treatment|Animal)
  contrast(test1, method = "pairwise") %>% as.data.frame()
}

#soil CN and total soil N model
mod_soilCN_soilN <- glmmTMB(CN_ratio_soil ~ TN_soil + treatment, data = comb_soil)
#CN and N in soil compared to GHG
NGPP <- glmmTMB(CO2_PS_flux ~ TN * treatment, data = comb_flux)
NReco <- glmmTMB(CO2_RE_flux ~ TN * treatment, data = comb_flux)
NCH4 <- glmmTMB(normalized_CH4_flux ~ TN * treatment, data = comb_flux)
NN2O <- glmmTMB(normalized_N2O_flux ~ TN * treatment, data = comb_flux)
CNGPP <- glmmTMB(CO2_PS_flux ~ CN_ratio * treatment, data = comb_flux)
CNReco <- glmmTMB(CO2_RE_flux ~ CN_ratio * treatment, data = comb_flux)
CNCH4 <- glmmTMB(normalized_CH4_flux ~ CN_ratio * treatment, data = comb_flux)
CNN2O <- glmmTMB(normalized_N2O_flux ~ CN_ratio * treatment, data = comb_flux)
# veg growth soil influence model - PAPER
vegsoil_model <- glmmTMB(height_value ~ Animal + TN_soil + PO4.P_soil + treatment, data = comb_soil)

run_model(comb_soil, vegsoil_model) #change accordingly

