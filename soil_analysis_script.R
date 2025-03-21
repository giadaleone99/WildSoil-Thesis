### Soil data analysis script 

library(ggplot2)
library(dplyr)
library(patchwork)
library(ARTool)
library(plotrix)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(effects)

soil_data_raw <- read.csv("data/soil_data_raw.csv")
date_data <- readRDS("flux_data/clean_flux_data.rds")

# get the dates per base code for days since first
date_data <- date_data %>% 
  filter(gastype == "CO2") %>% 
  filter(NEERE == "RE") %>% 
  filter(treatment == "F") %>% 
  filter(Days_Since_First == 0) %>% 
  select(base_code, longdate)

date_data <- date_data %>% 
  mutate(days_since_first = as.numeric(as.Date("2024-08-02") - as.Date(longdate))) %>% 
  mutate(dung_age = case_when(days_since_first %in% c(3, 4) ~ "3-4", days_since_first == 18 ~ "18", days_since_first %in% c(49, 50) ~ "49-50")) %>% 
  select(base_code, dung_age)

# create combined dataset with dung plot age
comb_soil_data <- soil_data_raw %>% 
  mutate(Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),) %>% 
  mutate(Animal = as.factor(Animal))

comb_soil_data$sample_type <- factor(comb_soil_data$sample_type, levels = c("Dung soil", "Fresh", "Control"))

comb_soil_data <- comb_soil_data %>% 
  left_join(date_data, by = "base_code")
  
# create subdatasets for the campaigns
gradient_soil <- soil_data_raw %>% 
  filter(grepl("G.", base_code)) %>% 
  mutate(Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),) %>% 
  mutate(Animal = as.factor(Animal))
  
gradient_soil$sample_type <- factor(gradient_soil$sample_type, levels = c("Dung soil", "Fresh", "Control"))

saveRDS(gradient_soil, "data/gradient_soil_data.rds")


daily_soil <- soil_data_raw %>% 
  filter(grepl("D.", base_code)) %>% 
  mutate(Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),) %>% 
  mutate(Animal = as.factor(Animal))

daily_soil$sample_type <- factor(daily_soil$sample_type, levels = c("Dung soil", "Fresh", "Control"))

saveRDS(daily_soil, "data/daily_soil_data.rds")

summary(gradient_soil)

combined_soil <- full_join(gradient_soil, daily_soil, by = intersect(names(gradient_soil), names(daily_soil)))

combined_soil <- combined_soil %>% 
  mutate(Campaign = as.factor(case_when(
    grepl("G", base_code) ~ "Gradient",
    grepl("D", base_code) ~ "Daily")
  ))


# Ensure sample_type_animal is a factor with the correct order
gradient_soil <- gradient_soil %>%
  mutate(
    sample_type_animal = interaction(sample_type, Animal))

gradient_soil$sample_type_animal <- factor(gradient_soil$sample_type_animal,
                       levels = c("Control.Cow", "Fresh.Cow", "Dung soil.Cow",
                                  "Control.Horse", "Fresh.Horse", "Dung soil.Horse"),ordered = TRUE)

# Create the plot
gradientphplot <- ggplot(gradient_soil, aes(x = Animal, y = pH, fill = sample_type_animal)) + 
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  xlab("\nAnimal") + ylab(NULL) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow control", "Cow beside dung", "Cow below dung", 
                               "Horse control", "Horse beside dung", "Horse below dung")) +
  expand_limits(y = c(4, 7)) +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) +
  ggtitle("B")

# Display the plot
gradientphplot

ggsave(filename = "soil_plots/gradientph.jpeg", plot = gradientphplot, width = 6, height = 4)

daily_soil <- daily_soil %>%
  mutate(
    sample_type_animal = interaction(sample_type, Animal))

daily_soil$sample_type_animal <- factor(daily_soil$sample_type_animal,
                                           levels = c("Control.Cow", "Fresh.Cow", "Dung soil.Cow",
                                                      "Control.Horse", "Fresh.Horse", "Dung soil.Horse"),ordered = TRUE)

dailyphplot <- ggplot(daily_soil, aes(x = Animal, y = pH, fill = sample_type_animal)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") +
  ylab("pH") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow below dung", "Cow beside dung", "Cow control", "Horse below dung", "Horse beside dung", "Horse control"))+
  expand_limits(y = c(4, 7)) +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) +
  ggtitle("A") +
  theme(legend.position = "none")
dailyphplot
ggsave(filename = "soil_plots/dailyph.jpeg", plot = dailyphplot, width = 6, height = 4)

gradientpo4plot <- ggplot(gradient_soil, aes(x = Animal, y = PO4.P, fill = sample_type_animal)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab(NULL) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow control", "Cow beside dung", "Cow below dung", 
                               "Horse control", "Horse beside dung", "Horse below dung")) +
  expand_limits(y = c(50, 370)) +
  scale_y_continuous(breaks = seq(50, 370, by = 50)) +
  ggtitle("B")
gradientpo4plot
ggsave(filename = "soil_plots/gradientpo4.jpeg", plot = gradientpo4plot, width = 6, height = 4)

dailypo4plot <- ggplot(daily_soil, aes(x = Animal, y = PO4.P, fill = sample_type_animal)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab("Plant available P (mg P/kg soil)") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow below dung", "Cow beside dung", "Cow control", "Horse below dung", "Horse beside dung", "Horse control"))+
  expand_limits(y = c(50, 370)) +
  scale_y_continuous(breaks = seq(50, 370, by = 50)) +
  ggtitle("A") +
  theme(legend.position = "none")
dailypo4plot
ggsave(filename = "soil_plots/dailypo4.jpeg", plot = dailypo4plot, width = 6, height = 4)

gradientcnplot <- ggplot(gradient_soil, aes(x = Animal, y = CN_ratio, fill = sample_type_animal)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab(NULL) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow control", "Cow beside dung", "Cow below dung", 
                               "Horse control", "Horse beside dung", "Horse below dung")) +
  expand_limits(y = c(11, 16)) +
  scale_y_continuous(breaks = seq(11, 16, by = 1)) +
  ggtitle("B")
gradientcnplot
ggsave(filename = "soil_plots/gradientcn.jpeg", plot = gradientcnplot, width = 6, height = 4)

dailycnplot <- ggplot(daily_soil, aes(x = Animal, y = CN_ratio, fill = sample_type_animal)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab("C:N ratio") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow below dung", "Cow beside dung", "Cow control", "Horse below dung", "Horse beside dung", "Horse control"))+
  expand_limits(y = c(11, 16)) +
  scale_y_continuous(breaks = seq(11, 16, by = 1)) +
  ggtitle("A") +
  theme(legend.position = "none")
dailycnplot
ggsave(filename = "soil_plots/dailycn.jpeg", plot = dailycnplot, width = 6, height = 4)

combinedphplot <- dailyphplot + gradientphplot
combinedphplot
ggsave(filename = "soil_plots/combinedph.jpeg", plot = combinedphplot, width = 8, height = 4)

combinedpo4plot <- dailypo4plot + gradientpo4plot
combinedpo4plot
ggsave(filename = "soil_plots/combinedpo4.jpeg", plot = combinedpo4plot, width = 8, height = 4)

combinedcnplot <- dailycnplot + gradientcnplot
combinedcnplot
ggsave(filename = "soil_plots/combinedcn.jpeg", plot = combinedcnplot, width = 8, height = 4)

# boxplots with dung age - IMPROVED FOR PAPER

comb_soil_data <- comb_soil_data %>%
  mutate(
    sample_type_animal = interaction(sample_type, Animal))
comb_soil_data$dung_age <- ordered(comb_soil_data$dung_age, levels = c("3-4", "18", "49-50"))

# CN ratio
cnplot <- ggplot(comb_soil_data, aes(x = Animal, y = CN_ratio, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days old", "18" = "18 days old", "49-50" = "49-50 days old"))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  ylab("C:N ratio") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow below dung", "Cow beside dung", "Cow control", "Horse below dung", "Horse beside dung", "Horse control"))+
  expand_limits(y = c(11, 16)) +
  scale_y_continuous(breaks = seq(11, 16, by = 1))
  #ggtitle("C:N ratio")
cnplot
ggsave(filename = "soil_plots/newcn.jpeg", plot = cnplot, width = 6, height = 4)

# pH
phplot <- ggplot(comb_soil_data, aes(x = Animal, y = pH, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days old", "18" = "18 days old", "49-50" = "49-50 days old"))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  ylab("pH") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow below dung", "Cow beside dung", "Cow control", "Horse below dung", "Horse beside dung", "Horse control")) +
  expand_limits(y = c(4, 7)) +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5))
  #ggtitle("pH")
phplot
ggsave(filename = "soil_plots/newph.jpeg", plot = phplot, width = 6, height = 4)

# PO4.P
po4plot <- ggplot(comb_soil_data, aes(x = Animal, y = PO4.P, fill = sample_type_animal)) +
  facet_wrap(~dung_age, labeller = labeller(dung_age = c("3-4" = "3-4 days old", "18" = "18 days old", "49-50" = "49-50 days old"))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  ylab("Plant available P (mg P/kg soil)") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) +
  labs(fill = "Plot type & \nsampling location") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow below dung", "Cow beside dung", "Cow control", "Horse below dung", "Horse beside dung", "Horse control"))+
  expand_limits(y = c(50, 370)) +
  scale_y_continuous(breaks = seq(50, 370, by = 50))
  #ggtitle("Plant available P")
po4plot
ggsave(filename = "soil_plots/newpo4p.jpeg", plot = po4plot, width = 6, height = 4)

## Basic statistics 

# Daily 
Daily_summary_stats <- daily_soil %>%
  group_by(Animal, sample_type) %>%
  summarise(
    pH_mean = mean(pH, na.rm = TRUE),
    pH_median = median(pH, na.rm = TRUE),
    pH_sd = sd(pH, na.rm = TRUE),
    pH_se = std.error(pH),
    PO4.P_mean = mean(PO4.P, na.rm = TRUE),
    PO4.P_median = median(PO4.P, na.rm = TRUE),
    PO4.P_sd = sd(PO4.P, na.rm = TRUE),
    PO4.P_se = std.error(PO4.P),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio)
  )


# Gradient
Gradient_summary_stats <- gradient_soil %>%
  group_by(Animal, sample_type) %>%
  summarise(
    pH_mean = mean(pH, na.rm = TRUE),
    pH_median = median(pH, na.rm = TRUE),
    pH_sd = sd(pH, na.rm = TRUE),
    pH_se = std.error(pH),
    PO4.P_mean = mean(PO4.P, na.rm = TRUE),
    PO4.P_median = median(PO4.P, na.rm = TRUE),
    PO4.P_sd = sd(PO4.P, na.rm = TRUE),
    PO4.P_se = std.error(PO4.P),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio)
  )


# ANOVA
# Daily pH 
daily_pH_ANOVA1 <- aov(pH ~ Animal + sample_type, data = daily_soil)
summary(daily_pH_ANOVA1)
res <- residuals(daily_pH_ANOVA1)
qqnorm(res)
qqline(res)
shapiro.test(res)

TukeyHSD(daily_pH_ANOVA1)
plot(TukeyHSD(daily_pH_ANOVA1))

# Interaction *
daily_pH_ANOVA2 <- aov(pH ~ Animal * sample_type, data = daily_soil)
summary(daily_pH_ANOVA2)
res <- residuals(daily_pH_ANOVA2)
qqnorm(res)
qqline(res)
shapiro.test(res)

TukeyHSD(daily_pH_ANOVA2)
plot(TukeyHSD(daily_pH_ANOVA2))

# Model comparison
anova(daily_pH_ANOVA1, daily_pH_ANOVA2)

# Daily PO4
daily_PO4.P_ANOVA1 <- aov(log_PO4.P ~ Animal + sample_type, data = daily_soil)
summary(daily_PO4.P_ANOVA1)
res <- residuals(daily_PO4.P_ANOVA1)
qqnorm(res)
qqline(res)
shapiro.test(res) # Residuals are not normally distributed to using ART model instead

# ART 
# Fit the ART model
daily_PO4.P_ART_model <- art(PO4.P ~ sample_type * Animal, data = daily_soil)
anova(daily_PO4.P_ART_model)

# Daily CN
daily_CN_ANOVA1 <- aov(CN_ratio ~ Animal + sample_type, data = daily_soil)
summary(daily_CN_ANOVA1)
res <- residuals(daily_CN_ANOVA1)
qqnorm(res)
qqline(res)
shapiro.test(res) # normally distributed

# Daily CN interaction *
daily_CN_ANOVA2 <- aov(CN_ratio ~ Animal * sample_type, data = daily_soil)
summary(daily_CN_ANOVA2)
res <- residuals(daily_CN_ANOVA2)
qqnorm(res)
qqline(res)
shapiro.test(res) # normally distributed

# Model comparison
anova(daily_CN_ANOVA1, daily_CN_ANOVA2)

# Gradient pH
gradient_pH_ANOVA1 <- aov(pH ~ Animal + sample_type, data = gradient_soil)
summary(gradient_pH_ANOVA1)
res <- residuals(gradient_pH_ANOVA1)
qqnorm(res)
qqline(res)
shapiro.test(res)

TukeyHSD(gradient_pH_ANOVA1)
plot(TukeyHSD(gradient_pH_ANOVA1))

# Gradient pH interaction *
gradient_pH_ANOVA2 <- aov(pH ~ Animal * sample_type, data = gradient_soil)
summary(gradient_pH_ANOVA2)
res <- residuals(gradient_pH_ANOVA2)
qqnorm(res)
qqline(res)
shapiro.test(res)

TukeyHSD(gradient_pH_ANOVA2)
plot(TukeyHSD(gradient_pH_ANOVA2))

# model comparison
anova(gradient_pH_ANOVA1, gradient_pH_ANOVA2)

# Gradient PO4.P
gradient_PO4.P_ANOVA1 <- aov(PO4.P ~ Animal + sample_type, data = gradient_soil)
summary(gradient_PO4.P_ANOVA1)
res <- residuals(gradient_PO4.P_ANOVA)
qqnorm(res)
qqline(res)
shapiro.test(res)

TukeyHSD(gradient_PO4.P_ANOVA1)
plot(TukeyHSD(gradient_PO4.P_ANOVA))

#Gradient PO4.P interaction *
gradient_PO4.P_ANOVA2 <- aov(PO4.P ~ Animal * sample_type, data = gradient_soil)
summary(gradient_PO4.P_ANOVA)
res <- residuals(gradient_PO4.P_ANOVA)
qqnorm(res)
qqline(res)
shapiro.test(res)

TukeyHSD(gradient_PO4.P_ANOVA)
plot(TukeyHSD(gradient_PO4.P_ANOVA))

# model comparison
anova(gradient_PO4.P_ANOVA1, gradient_PO4.P_ANOVA2)

# Gradient CN
gradient_CN_ANOVA1 <- aov(CN_ratio ~ Animal + sample_type, data = gradient_soil)
summary(gradient_CN_ANOVA1)
res <- residuals(gradient_CN_ANOVA1)
qqnorm(res)
qqline(res)
shapiro.test(res)

#Gradient CN interaction +
gradient_CN_ANOVA2 <- aov(CN_ratio ~ Animal * sample_type, data = gradient_soil)
summary(gradient_CN_ANOVA)
res <- residuals(gradient_CN_ANOVA)
qqnorm(res)
qqline(res)
shapiro.test(res)

anova(gradient_CN_ANOVA1, gradient_CN_ANOVA2)

## We tried separating Horse and Cow and testing sample_type but no significance

# Comparing gradient and daily campaigns
gradient_stats_fresh <- gradient_soil %>%
  filter(sample_type %in% c("Fresh", "Dung soil")) %>%
  summarise(
    mean_CN = mean(CN_ratio, na.rm = TRUE),
    std_error_CN = std.error(CN_ratio, na.rm = TRUE),
    mean_pH = mean(pH, na.rm = TRUE),
    std_error_pH = std.error(pH, na.rm = TRUE),
    mean_PO4 = mean(PO4.P, na.rm = TRUE),
    std_error_PO4 = std.error(PO4.P, na.rm = TRUE)
  ) %>%
  mutate(source = "gradient_soil")

# For daily_soil: filter for both fresh and dung samples
daily_stats_fresh <- daily_soil %>%
  filter(sample_type %in% c("Fresh", "Dung soil")) %>%
  summarise(
    mean_CN = mean(CN_ratio, na.rm = TRUE),
    std_error_CN = std.error(CN_ratio, na.rm = TRUE),
    mean_pH = mean(pH, na.rm = TRUE),
    std_error_pH = std.error(pH, na.rm = TRUE),
    mean_PO4 = mean(PO4.P, na.rm = TRUE),
    std_error_PO4 = std.error(PO4.P, na.rm = TRUE)
  ) %>%
  mutate(source = "daily_soil")

combined_stats <- bind_rows(gradient_stats_fresh, daily_stats_fresh)

# Lasses's method 

# Function for modelling 
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
  contrast(test3, method = "pairwise") %>% as.data.frame()
}

# Soil parameters
Soil_CNratio_m1 <- glmmTMB(CN_ratio ~ Animal * sample_type * Campaign + (1|base_code), data = combined_soil)

Soil_pH_m1 <- glmmTMB(pH ~ Animal * sample_type * Campaign + (1|base_code), data = combined_soil)

Soil_PO4.P_m1 <- glmmTMB(PO4.P ~ Animal * sample_type * Campaign + (1|base_code), data = combined_soil)

# With bulk density as a fixed effect
Soil_CNratio_b1 <- glmmTMB(CN_ratio ~ Animal * sample_type * Campaign * bulk_density + (1|base_code), data = combined_soil)

Soil_pH_b1 <- glmmTMB(pH ~ Animal * sample_type * Campaign * bulk_density + (1|base_code), data = combined_soil)

Soil_PO4.P_b1 <- glmmTMB(PO4.P ~ Animal * sample_type * Campaign * bulk_density + (1|base_code), data = combined_soil)

#models with dung age - IMPROVED FOR PAPER
Soil_CNratio_d1 <- glmmTMB(CN_ratio ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

Soil_pH_d1 <- glmmTMB(pH ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

Soil_PO4.P_d1 <- glmmTMB(PO4.P ~ Animal * sample_type * dung_age + (1|base_code), data = comb_soil_data)

# Change this accordingly 
run_model(comb_soil_data, Soil_CNratio_d1)

# model effects
Soil_CNratio_effects <- allEffects(Soil_CNratio_m1)
Soil_pH_effects <- allEffects(Soil_pH_m1)
Soil_PO4.P_effects <- allEffects(Soil_PO4.P_m1)
print(Soil_CNratio_effects)
plot(Soil_CNratio_effects)
