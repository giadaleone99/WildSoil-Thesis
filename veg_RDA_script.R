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
species_list <- read.csv("data/species_lists.csv")

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

# Create a presence absence matrix for another RDA
sp_presence <- species_list %>% 
  mutate(presence = 1) 

# Get a full list of the species
unique_species <- unique(species_list$species_list)

# Removing the Bryophytes
unique_species <- unique_species[unique_species != "Bryophytes"]

# Create presence absence matrix
presence_absence_matrix <- sp_presence %>%
  dplyr::select(plot_id, species_list, presence) %>%  
  pivot_wider(names_from = species_list, values_from = presence, values_fill = list(presence = 0)) %>% 
  dplyr::select(-Bryophytes)

presence_absence_matrix <- presence_absence_matrix[order(presence_absence_matrix$plot_id), ]
presence_absence_matrix <- presence_absence_matrix[, !names(presence_absence_matrix) %in% "plot_id"]

# Create species weight matrix and deleting Bryophytes
species_df <- species_df %>% 
  pivot_wider(names_from = veg_class, values_from = weight, values_fill = 0) %>% 
  dplyr::select(-c("Bryophytes", "Graminoids")) 

species_df <- species_df[order(species_df$plot_id), ]
species_df <- species_df[, !names(species_df) %in% "plot_id"]

species_df <- species_df %>% 
  round(digits = 4) # Rounding all numbers to 4 decimal places

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


# RDA with weights ----------------------------------------------------------------------------------
species_df.hell <- decostand(species_df, "hellinger") #hellinger-transform the species dataset

#rda_model <- rda(species_df ~ ., soil_df)
rda_model <- rda(species_df.hell ~ bulk_density + pH + PO4.P + CN_ratio + Campaign + sample_type + Animal, data = soil_df)
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

# Variation Partitioning ####
bulk <- soil_df$bulk_density
ph <- soil_df$pH
cn <- soil_df$CN_ratio
PO4.P <- soil_df$PO4.P
soil <- soil_df[c("bulk_density", "pH", "CN_ratio", "PO4.P")]
campaign <- soil_df$Campaign
sample_type <- soil_df$sample_type
Animal <- soil_df$Animal


# Soil
(spe.part.all <- varpart(species_df.hell, PO4.P, ph, cn, bulk))

jpeg(file="veg_plots/Soil_variation_partitioning.jpeg")
plot(spe.part.all, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("PO4.P","pH", "C:N", "Bulk"))
dev.off()

# Test PO4
anova(rda(species_df.hell, PO4.P, permutation = how(nperm = 9999)))

# Fraction [a], pure PO4
anova(rda(species_df.hell, PO4.P, cbind(ph, cn, bulk)), permutation = how(nperm = 9999))

# Test pH
anova(rda(species_df.hell, ph, permutation = how(nperm = 9999)))

# Fraction [b], pure pH
anova(rda(species_df.hell, ph, cbind(PO4.P, cn, bulk)), permutation = how(nperm = 9999))

# Test C:N
anova(rda(species_df.hell, cn, permutation = how(nperm = 9999)))

# Fraction [c], pure CN
anova(rda(species_df.hell, cn, cbind(PO4.P, ph, bulk)), permutation = how(nperm = 9999))

# Test bulk density
anova(rda(species_df.hell, bulk, permutation = how(nperm = 9999)))

# Fraction [d], pure bulk density
anova(rda(species_df.hell, bulk, cbind(PO4.P, ph, cn)), permutation = how(nperm = 9999))

# Soil and other variables
(spe.part.all2 <- varpart(species_df.hell, soil, campaign, sample_type, Animal))
jpeg(file="veg_plots/variation_partitioning.jpeg")
plot(spe.part.all2, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("Soil","Campaign", "Treatment", "Animal"))
dev.off()

# Test soil
anova(rda(species_df.hell, soil, permutation = how(nperm = 9999)))

# Fraction [a], pure soil
anova(rda(species_df.hell, soil, cbind(campaign, sample_type, Animal)), permutation = how(nperm = 9999))

# Test campaign
anova(rda(species_df.hell, campaign, permutation = how(nperm = 9999)))

# Fraction [a] only campaign
anova(rda(species_df.hell, campaign, cbind(soil, sample_type, Animal)), permutation = how(nperm = 9999))

# Test sample_type
anova(rda(species_df.hell, sample_type, permutation = how(nperm = 9999)))

# Fraction [a] only sample_type
anova(rda(species_df.hell, sample_type, cbind(soil, campaign, Animal)), permutation = how(nperm = 9999))

# Test Animal
anova(rda(species_df.hell, Animal, permutation = how(nperm = 9999)))

# Fraction [a] only Animal
anova(rda(species_df.hell, Animal, cbind(soil, campaign, sample_type)), permutation = how(nperm = 9999))

#Not from Lasse ----
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


## Lasses way ---------
# RDA with the presence absence of forbs and graminoids -------------------
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
cn <- soil_df$CN_ratio
PO4.P <- soil_df$PO4.P
soil <- soil_df[c("bulk_density", "pH", "CN_ratio", "PO4.P")]
campaign <- soil_df$Campaign
sample_type <- soil_df$sample_type
Animal <- soil_df$Animal


# Soil
(spe.part.all3 <- varpart(presence_absence_matrix, PO4.P, ph, cn, bulk))

jpeg(file="veg_plots/Soil_variation_partitioning_absence.jpeg")
plot(spe.part.all3, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("PO4.P","pH", "C:N", "Bulk"))
dev.off()

# Test PO4
anova(rda(presence_absence_matrix, PO4.P, permutation = how(nperm = 9999)))

# Fraction [a], pure PO4
anova(rda(presence_absence_matrix, PO4.P, cbind(ph, cn, bulk)), permutation = how(nperm = 9999))

# Test pH
anova(rda(presence_absence_matrix, ph, permutation = how(nperm = 9999)))

# Fraction [b], pure pH
anova(rda(presence_absence_matrix, ph, cbind(PO4.P, cn, bulk)), permutation = how(nperm = 9999))

# Test C:N
anova(rda(presence_absence_matrix, cn, permutation = how(nperm = 9999)))

# Fraction [c], pure CN
anova(rda(presence_absence_matrix, cn, cbind(PO4.P, ph, bulk)), permutation = how(nperm = 9999))

# Test bulk density
anova(rda(presence_absence_matrix, bulk, permutation = how(nperm = 9999)))

# Fraction [d], pure bulk density
anova(rda(presence_absence_matrix, bulk, cbind(PO4.P, ph, cn)), permutation = how(nperm = 9999))

# Soil and other variables
(spe.part.all4 <- varpart(presence_absence_matrix, soil, campaign, sample_type, Animal))
jpeg(file="veg_plots/variation_partitioning_absence.jpeg")
plot(spe.part.all4, digits = 2, bg = c("red", "blue", "green", "yellow"),
     Xnames = c("Soil","Campaign", "Treatment", "Animal"))
dev.off()

# Test soil
anova(rda(presence_absence_matrix, soil, permutation = how(nperm = 9999)))

# Fraction [a], pure soil
anova(rda(presence_absence_matrix, soil, cbind(campaign, sample_type, Animal)), permutation = how(nperm = 9999))

# Test campaign
anova(rda(presence_absence_matrix, campaign, permutation = how(nperm = 9999)))

# Fraction [a] only campaign
anova(rda(presence_absence_matrix, campaign, cbind(soil, sample_type, Animal)), permutation = how(nperm = 9999))

# Test sample_type
anova(rda(presence_absence_matrix, sample_type, permutation = how(nperm = 9999)))

# Fraction [a] only sample_type
anova(rda(presence_absence_matrix, sample_type, cbind(soil, campaign, Animal)), permutation = how(nperm = 9999))

# Test Animal
anova(rda(presence_absence_matrix, Animal, permutation = how(nperm = 9999)))

# Fraction [a] only Animal
anova(rda(presence_absence_matrix, Animal, cbind(soil, campaign, sample_type)), permutation = how(nperm = 9999))



