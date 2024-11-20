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

veg_raw <- read.csv("data/vegetation_data.csv")
soil_data_raw <- read.csv("data/soil_data_raw.csv")

# get the veg species data and transform it into a matrix
veg_new <- veg_raw %>%
  mutate(
    `Total harvested weight` = dry_weight / harvested_area,
    treatment = case_when(
      grepl("_F$", plot_id) ~ "Fresh",
      grepl("_C$", plot_id) ~ "Control",
      TRUE ~ NA_character_
    ),
    Campaign = case_when(
      grepl("G", plot_id) ~ "Gradient",
      grepl("D", plot_id) ~ "Daily"
    ),
    Animal = case_when(grepl("^C", plot_id) ~ "Cow", grepl("^H", plot_id) ~ "Horse"),
    harvested_area_m2 = harvested_area / 10000
  ) %>% 
  mutate(
    weight = dry_weight/harvested_area_m2)

species_df <- veg_new[c("plot_id", "veg_class", "weight")]

species_df <- species_df %>% 
  pivot_wider(names_from = veg_class, values_from = weight, values_fill = 0)

species_df <- species_df[order(species_df$plot_id), ]

species_df <- species_df[, !names(species_df) %in% "plot_id"]

# get the relevant soil data
soil_data <- soil_data_raw %>% 
  filter(sample_type %in% c("Fresh", "Control")) %>% 
  mutate(Campaign = as.factor(case_when(grepl("G", plotID) ~ "Gradient",grepl("D", plotID) ~ "Daily")),
  Animal = as.factor(case_when(grepl("^C", plotID) ~ "Cow", grepl("^H", plotID) ~ "Horse")),
  sample_type = as.factor(sample_type))
  
soil_data <- soil_data[c("plotID", "bulk_density", "pH", "PO4.P", "CN_ratio", "Campaign", "Animal", "sample_type")]

soil_df <- soil_data[order(soil_data$plotID), ]
rownames(soil_df) <- NULL

soil_df <- soil_df[c("bulk_density", "pH", "PO4.P", "CN_ratio", "Campaign", "Animal", "sample_type")]

# RDA
rda_model <- rda(species_df ~ ., soil_df)
rda_model <- rda(species_df ~ bulk_density + pH + PO4.P + CN_ratio + Campaign + sample_type + Animal, data = soil_df)
summary_rda <- summary(rda_model)
summary_rda
# Checking the R2 values of the separate variables to see which ones are not explaining the data
eigenvalues <- summary_rda$cont$importance[2, ]
total_variance <- sum(eigenvalues)
r2_individual <- eigenvalues / total_variance
print(r2_individual)
names(r2_individual) <- c("bulk_density", "pH", "PO4.P", "CN_ratio", "Campaign", "sample_type", "Animal")
print(r2_individual)

# Check for multicollinearity - not working
vif_model <- lm(species_df ~ bulk_density + pH + PO4.P + CN_ratio + Campaign + sample_type + Animal, data = soil_df)

# total R2 value
(r2 <- RsquareAdj(rda_model)$adj.r.squared)
print(r2)
plot(rda_model, scaling = 3)

# model P value
anova_result <- anova.cca(rda_model)
print(anova_result)



# get P values for the variables
anova_terms <- anova.cca(rda_model, by = "term")
print(anova_terms)

# get P values for something...
anova_axes <- anova.cca(rda_model, by = "axis")
print(anova_axes)
