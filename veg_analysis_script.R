#### Analyzing Vegetation Data from Mols
# Libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(gridExtra)
library(patchwork)
library(cowplot)
library(tidyr)
library(ggpattern)
#library(ARTool)
library(lme4)
library(lmerTest)
library(plotrix)
library(DHARMa)
library(emmeans)
library(car)
library(glmmTMB)
library(sjPlot)
library(vegan)
library(paletteer)

# Import data
vegdung_lab <- read.csv("data/Plant_Dung_CN_pH_elements.csv", sep = ";")
veg_raw <- read.csv("data/vegetation_data.csv")
fieldwork_data_raw <- read.csv("data/Fieldwork_data_final.csv")
species_list <- read.csv("data/species_lists.csv")

#dung_data <- readRDS("data/dung_data.rds")
#gradient_soil_data <- readRDS("data/gradient_soil_data.rds")
#daily_soil_data <- readRDS("data/daily_soil_data.rds")
#species_data <- readRDS("data/species_data.rds")
#vegetation_data <- readRDS("data/vegetation_data.rds")

# Manipulating vegetation data
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
    )
  )

veg_new <- veg_new %>% 
  mutate(veg_category = case_when(
    veg_class == "Bryophytes" ~ "Bryophytes",
    veg_class == "Graminoids" ~ "Graminoids",
    TRUE ~ "Forbs"
    
  ))

veg_summary <- veg_new %>%
  group_by(plot_id) %>%
  dplyr::summarise(
    total_veg_weight = sum(`Total harvested weight`),
    veg_height = unique(veg_height),
    veg_height_2 = unique(veg_height_2)
  ) %>%
  mutate(
    base_code = sub("_[FC]$", "", plot_id),
    Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),
    treatment = case_when(
      grepl("_F$", plot_id) ~ "Fresh",
      grepl("_C$", plot_id) ~ "Control",
      TRUE ~ NA_character_
    ),
    Animal_treatment = paste(Animal, treatment),
    Campaign = case_when(
      grepl("G", base_code) ~ "Gradient",
      grepl("D", base_code) ~ "Daily"
    )
  )

# Fieldwork data manipulation and joining with vegetation data
fieldwork_data <- fieldwork_data_raw %>%
  mutate(plot_id = str_remove(Plot_ID, "_[^_]+$"))

veg_combined <- veg_summary %>%
  left_join(fieldwork_data %>% select(plot_id, dung_area_cm2), by = "plot_id", multiple = "first") %>%
  mutate(
    frame_area_cm2 = 3058.15,
    area_minus_dung = frame_area_cm2 - dung_area_cm2
  ) %>%
  left_join(select(veg_new, plot_id, harvested_area), by = "plot_id") %>%
  distinct() %>%
  mutate(biomass = total_veg_weight * area_minus_dung,
         Animal = as.factor(Animal),
         treatment = as.factor(treatment))

# get CN and plot_id from the lab data sheet and merge with the rest
vegdung_data <- vegdung_lab %>% 
  select(Kode_1, CN_ratio) %>% 
  filter(!Kode_1 %in% c("CG dung", "HG dung", "HOD dung", "COD dung", "HND dung", "CND dung"))
colnames(vegdung_data) <- c("plot_id", "CN_ratio")

veg_combined <- veg_combined %>% 
  left_join(vegdung_data, by = "plot_id")

# create a df for dung
dung_data <- vegdung_lab %>% 
  select(Kode_1, CN_ratio) %>% 
  filter(Kode_1 %in% c("CG dung", "HG dung", "HOD dung", "COD dung", "HND dung", "CND dung"))
colnames(vegdung_data) <- c("plot_id", "CN_ratio")

# mutating the data for the stacked height bar plots
veg_combined2 <- veg_combined
na_index <- is.na(veg_combined2$veg_height) & !is.na(veg_combined2$veg_height_2)
veg_combined2$veg_height[na_index] <- veg_combined2$veg_height_2[na_index]
veg_combined2$veg_height_2[na_index] <- NA


veg_heightdata <- pivot_longer(veg_combined2, cols = c("veg_height", "veg_height_2"), 
                               names_to = "height_type", 
                               values_to = "height_value")
# Filtering for veg height growth and bar plot
veg_growth <- veg_summary %>%
  filter(!is.na(veg_height)) %>%
  mutate(height_value = veg_height_2 - veg_height,
         height_type = "veg_growth") 


# Combine the original DataFrame with the new rows
veg_heightdata <- bind_rows(veg_heightdata, veg_growth) %>% 
  filter(height_type != "veg_height_2")
  


#make subsets for the campaigns
daily_veg_combined <- veg_combined %>% 
  filter(Campaign == "Daily")
daily_veg_growth <- veg_growth %>% 
  filter(Campaign == "Daily")
gradient_veg_combined <- veg_combined %>% 
  filter(Campaign == "Gradient")
gradient_veg_growth <- veg_growth %>% 
  filter(Campaign == "Gradient")

##-----------------------------------------------------------------------------

# Plotting function for weight and height
save_plot <- function(df, y_var, filename, y_label) {
  p <- ggplot(df, aes(x = plot_id, y = !!sym(y_var), fill = treatment)) +
    geom_bar(stat = "identity") +
    xlab("Plot ID") + ylab(y_label) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
  ggsave(filename = filename, plot = p, width = 6, height = 4)
}

veg_filtered <- veg_combined %>% filter(grepl("_F$", plot_id))

save_plot(veg_filtered, "total_veg_weight", "veg_plots/Dung_veg_weight.jpeg", "Total Standardized Veg Weight")
save_plot(veg_filtered, "veg_height_2", "veg_plots/Dung_veg_height.jpeg", "Veg Height")

# plots for veg weight grouped together per campaign
veg_daily <- veg_combined %>% filter(grepl("D.", base_code))
dailyvegweight <- ggplot(veg_daily, aes(x = plot_id, y = total_veg_weight, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free") +
  xlab("\nPlot ID") + ylab("vegetation weight (grams)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetation weight of the short term campaign")
dailyvegweight
ggsave(filename = "veg_plots/dailyvegweight.jpeg", plot = dailyvegweight, width = 6, height = 4)

veg_gradient <- veg_combined %>% filter(grepl("G.", base_code))
gradientvegweight <- ggplot(veg_gradient, aes(x = plot_id, y = total_veg_weight, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free") +
  xlab("\nPlot ID") + ylab("vegetation weight (grams)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetation weight of the long term campaign")
gradientvegweight
ggsave(filename = "veg_plots/gradientvegweight.jpeg", plot = gradientvegweight, width = 6, height = 4)

# stacked bar plot for veg height 1 and 2 per campaign
veg_heightdatadaily <- veg_heightdata %>% filter(grepl("D.", base_code))
veg_heightdatadaily$height_type <- factor(veg_heightdatadaily$height_type, 
                                          levels = c("veg_growth", "veg_height"))

dailyvegheights <- ggplot(veg_heightdatadaily, aes(x = plot_id, y = height_value, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  geom_bar_pattern(aes(pattern = height_type), 
                   position = "stack", 
                   stat = "identity", 
                   pattern_fill = "black",     
                   pattern_density = 0.1, 
                   pattern_spacing = 0.03) +
  facet_wrap("Animal", scales = "free_x") +
  xlab("\nPlot ID") + 
  ylab("Vegetation height (cm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +  
  scale_pattern_manual(name = "Measurement", 
                       values = c("veg_height" = "none", "veg_growth" = "stripe"),
                       labels = c("Height", "Growth"), 
                       guide = guide_legend(override.aes = list(
                         fill = "transparent", 
                         pattern_fill = c("black", "transparent"),
                         pattern = c("none", "stripe"),  # Ensure correct pattern display
                         color = "black"  # Border color for legend boxes
                       ))) + 
  #guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  theme(legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5), expand = c(0, 0)) +  ggtitle("Vegetation height and growth of the short term campaign")

dailyvegheights
ggsave(filename = "veg_plots/veg_growth_daily_stacked.jpeg", plot = dailyvegheights, width = 6, height = 4)


veg_heightdatagradient <- veg_heightdata %>% filter(grepl("G.", base_code))
veg_heightdatagradient$height_type <- factor(veg_heightdatagradient$height_type, 
                                             levels = c("veg_growth", "veg_height"))


gradientvegheights <- ggplot(veg_heightdatagradient, aes(x = plot_id, y = height_value, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  geom_bar_pattern(aes(pattern = height_type), 
                   position = "stack", 
                   stat = "identity", 
                   pattern_fill = "black",     
                   pattern_density = 0.1, 
                   pattern_spacing = 0.03) +
  facet_wrap("Animal", scales = "free_x") +  # Changed to 'free_x' to match the daily plot
  xlab("\nPlot ID") + 
  ylab("Vegetation height (cm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +  
  scale_pattern_manual(name = "Measurement period", 
                       values = c("veg_height" = "none", "veg_growth" = "stripe"),
                       labels = c("Height", "Growth"), 
                       guide = guide_legend(override.aes = list(
                         fill = "transparent", 
                         pattern_fill = c("black", "transparent"),
                         pattern = c("none", "stripe"),  
                         color = "black"  
                       ))) + 
  theme(legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), expand = c(0, 0)) +
  ggtitle("Vegetation height and growth of the long term campaign")


gradientvegheights
ggsave(filename = "veg_plots/veg_growth_gradient_stacked.jpeg", plot = gradientvegheights, width = 6, height = 4)

# Biomass per campaign barcharts
veg_daily <- veg_combined %>% filter(grepl("D.", base_code))
dailyvegbiomass <- ggplot(veg_daily, aes(x = plot_id, y = biomass, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free_x") +
  xlab("\nPlot ID") + ylab("Estimated biomass per plot (g)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetation biomass per plot of the short term campaign")
dailyvegbiomass
ggsave(filename = "veg_plots/dailyvegbiomass.jpeg", plot = dailyvegbiomass, width = 6, height = 4)


veg_gradient <- veg_combined %>% filter(grepl("G.", base_code))
gradientvegbiomass <- ggplot(veg_gradient, aes(x = plot_id, y = biomass, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free_x") +
  xlab("\nPlot ID") + ylab("Estimated biomass per plot (g)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetation biomass per plot of the long term campaign")
gradientvegbiomass
ggsave(filename = "veg_plots/gradientvegbiomass.jpeg", plot = gradientvegbiomass, width = 6, height = 4)


veg_gradient$plot_id <- as.factor(veg_gradient$plot_id)

ggplot(veg_gradient, aes(x = plot_id, y = CN_ratio)) +
  geom_bar(stat = "identity") 

barplot(veg_gradient$biomass)
# Biomass per campaign boxplots
dailyvegbiomassbox <- ggplot(veg_daily, aes(x = Animal, y = biomass, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  xlab("\nAnimal") + ylab("Estimated biomass (g)") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung","Horse control", "Horse dung"))+
  ggtitle("Short term campaign") +
  scale_y_continuous(limits = c(20, 140), breaks = seq(20, 140, by = 20)) +
  theme(legend.position = "none")
dailyvegbiomassbox

gradientvegbiomassbox <- ggplot(veg_gradient, aes(x = Animal, y = biomass, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  xlab("\nAnimal") + ylab(NULL) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung","Horse control", "Horse dung"))+
  scale_y_continuous(limits = c(20, 140), breaks = seq(20, 140, by = 20)) +
  ggtitle("Long term campaign")
gradientvegbiomassbox

vegbiomassbox <- dailyvegbiomassbox + gradientvegbiomassbox
vegbiomassbox
ggsave(filename = "veg_plots/vegbiomassbox.jpeg", plot = vegbiomassbox, width = 6, height = 4)

# Percent stacked barplot for veg weight per species
veg_class_order <- c("Bryophytes", "Graminoids", "Achillea millefolium", "Campanula rotundifolia", "Cerastium fontanum", "Dandelions and false dandelions", "Daucus carota", 
                     "Euphrasia stricta", "Galium verum", "Hypericum perforatum", "Jacobaea vulgaris", "Knautia arvensis", "Pilosella officinarum", "Pimpinella saxifraga", 
                     "Plantago lanceolata", "Ranunculus sp.", "Rosa canina", "Rumex acetosa", "Rumex acetosella", "Stellaria graminea", "Trifolium arvense", 
                     "Trifolium campestre", "Trifolium pratense", "Trifolium repens", "Veronica chamaedrys", "Veronica officinalis", "Vicia sativa")

veg_new <- veg_new %>% mutate(veg_class = factor(veg_class, levels = veg_class_order))

lookup_species <- with(species_data, setNames(species_per_vegclass, paste(plot_id, veg_category, sep = "_")))

veg_new <- veg_new %>%
  dplyr::mutate(
    species_per_vegclass = lookup_species[paste(plot_id, veg_category, sep = "_")]
  )

percentspecies <- ggplot(veg_new, aes(x = plot_id, y = dry_weight, fill = veg_class)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap("Campaign", scales = "free") +
  xlab("\nPlot ID") + ylab("Proportion of total weight") +
  theme_minimal() +
  scale_fill_paletteer_d("Polychrome::palette36") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), legend.position = "bottom") +
  labs(fill = "Vegetation species/groups") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Proportions of vegetation species/groups per campaign")

percentspecies

ggsave(filename = "veg_plots/percentspecies.jpeg", plot = percentspecies, width = 11, height = 8)



# Scatter plot of veg weight vs height
scatterplot <- ggplot(veg_combined, aes(x = biomass, y = veg_height_2, color = treatment)) +
  geom_point() + geom_smooth(method = lm)
print(scatterplot)

gradient_regression <- ggplot(veg_gradient, aes(x = biomass, y = veg_height_2, color = treatment)) +
  geom_point() + geom_smooth(method = lm) +
  theme_minimal() +
  xlab("Estimated biomass per plot (g)") +
  ylab("Vegetation height (cm)") +
  ggtitle("Gradient") +
  scale_color_manual(values = c("Dung" = "black", "Control" = "gray")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        legend.position = "none") 

print(gradient_regression)

daily_regression <- ggplot(veg_daily, aes(x = biomass, y = veg_height_2, color = treatment)) +
  geom_point() + geom_smooth(method = lm) +
  theme_minimal() +
  xlab("Estimated biomass per plot (g)") +
  ylab("Vegetation height (cm)") +
  ggtitle("Daily") +
  labs(colour = "Treatment") +
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray"),
                     labels = c("Control", "Dung")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line())

print(daily_regression)

combined_regression <- gradient_regression + daily_regression 
print(combined_regression)
ggsave(filename = "veg_plots/combined_regression.jpeg", plot = combined_regression, width = 6, height = 4)

# Fitting models to plot linear regressions between biomass and veg_height_2 with SjPlot
# Fit linear models for each dataset
SjPlot_model_gradient <- glmmTMB(veg_height_2 ~ biomass * treatment + (1|base_code), data = veg_gradient)

SjPlot_model_daily <- lm(veg_height_2 ~ biomass * treatment , data = veg_daily)


# Effect plot with predicted values for gradient model
plot_gradient_effect <- plot_model(SjPlot_model_gradient, type = "pred", 
                                   terms = c("biomass", "treatment"), 
                                   title = "Gradient") +
  theme_minimal() +  
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray")) +
  scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  labs(x = "Estimated biomass per plot (g)", y = "Predicted vegetation height (cm)",
       colour = "Treatment") + 
  scale_y_continuous(limits = c(0, 30)) + 
  theme(
    panel.grid.minor = element_blank(),  
    panel.grid.major.x = element_blank(),  
    axis.line = element_line(color = "black",),
    axis.title.y = element_blank(),
  )

plot_gradient_effect

# Effect plot with predicted values for daily model
plot_daily_effect <- plot_model(SjPlot_model_daily, type = "pred", 
                                terms = c("biomass", "treatment"), 
                                title = "Daily") +
  theme_minimal() +  
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray"),
                     labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  labs(x = "Estimated biomass per plot (g)", y = "Predicted vegetation height (cm)",
       colour = "Treatment") + 
  scale_y_continuous(limits = c(0, 15)) + 
  theme(
    panel.grid.minor = element_blank(),  
    panel.grid.major.x = element_blank(),  
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

combined_model_plot <- plot_daily_effect + plot_gradient_effect
combined_model_plot
ggsave("veg_plots/model_biomass_vegheight_combined.jpeg", plot = combined_model_plot, width = 6, height = 4)

scatterplot <- ggplot(veg_combined, aes(x = total_veg_weight, y = veg_height_2, color = treatment)) +
  geom_point(size = 2) +  
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("Control" = "grey", "Fresh" = "black")) +  
  labs(color = "treatment") + 
  labs(
    x = "Total vegetation weight (g)",  
    y = "Vegetation height (cm)",        
    color = "Treatment"                   
  ) +
  theme_minimal()  

print(scatterplot)
ggsave("veg_plots/scatterplot_veg_weight_height.jpeg", plot = scatterplot, width = 6, height = 4)

save_plot(veg_growth, "veg_growth", "veg_plots/veg_growth.jpeg", "Veg Growth")

# Bar plot for biomass
save_plot(veg_combined, "estimated_biomass_plot", "veg_plots/estimated_biomass_plot.jpeg", "Estimated Biomass")

# Custom colors for each combination of treatment and animal
custom_colors <- c(
  "Control.Cow" = "#A4AC86",    # Cow Control
  "Fresh.Cow" = "#656D4A",      # Cow Fresh
  "Control.Horse" = "#A68A64",  # Horse Control
  "Fresh.Horse" = "#7F4F24"     # Horse Fresh
)

custom_labels <- c(
  "Control.Cow" = "Cow Control",
  "Fresh.Cow" = "Cow Fresh",
  "Control.Horse" = "Horse Control",
  "Fresh.Horse" = "Horse Fresh"
)

plot_by_code <- function(df, y_var, folder, y_label) {
  for (code in unique(df$base_code)) {
    df_subset <- df %>% filter(base_code == code)
    
    # Check if df_subset is created
    if (nrow(df_subset) == 0) {
      print(paste("No data for base_code:", code))
      next  # Skip to the next iteration if no data is found
    }
    
    # Create the plot
    p <- ggplot(df_subset, aes(x = plot_id, y = !!sym(y_var), fill = interaction(treatment, Animal))) +
      geom_bar(stat = "identity") +
      ggtitle(paste("Bar Plot for", code)) +
      scale_fill_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels
      labs(fill = "Plot Type",  # Change legend title to "Plot Type"
           x = "Plot ID",
           y = y_label) +  # Set dynamic y-axis label
      theme_minimal()
    
    # Save the plot
    ggsave(filename = paste0("veg_plots/", folder, "/", code, "_plot.jpeg"), plot = p, width = 6, height = 4)
  }
}


plot_by_code(veg_combined, "total_veg_weight", "veg_weight", "Vegetation weight (grams)")
plot_by_code(veg_growth, "veg_growth", "veg_growth", "Vegetation growth (cm)")

# CN plots
gradientcnplot <- ggplot(veg_gradient, aes(x = Animal, y = CN_ratio, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab(NULL) +
  scale_y_continuous(limits = c(12, 38), breaks = seq(12, 38, by = 4)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung","Horse control", "Horse dung"))+
  ggtitle("Long term campaign")
gradientcnplot
ggsave(filename = "veg_plots/gradientcn.jpeg", plot = gradientcnplot, width = 6, height = 4)

dailycnplot <- ggplot(veg_daily, aes(x = Animal, y = CN_ratio, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab("CN ratio") +
  scale_y_continuous(limits = c(12, 38), breaks = seq(12, 38, by = 4)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  ggtitle("Short term campaign") +
  theme(legend.position = "none")
dailycnplot
ggsave(filename = "veg_plots/dailycn.jpeg", plot = dailycnplot, width = 6, height = 4)

combinedcnplot <- dailycnplot + gradientcnplot
combinedcnplot
ggsave(filename = "veg_plots/combinedcn.jpeg", plot = combinedcnplot, width = 8, height = 4)

### STATISTICS ------------------------------------------------------------------------

#residuals are not normally distributed so log transform the data
veg_growth <- veg_growth %>% 
  mutate(log_veg_growth = log(2+veg_growth))
veg_combined <- veg_combined %>% 
  mutate(log_total_veg_weight = log(total_veg_weight),
         log_estimated_biomass_plot = log(biomass))
gradient_veg_growth <- gradient_veg_growth %>% 
  mutate(log_veg_growth = log(2+veg_growth))
gradient_veg_combined <- gradient_veg_combined %>% 
  mutate(log_total_veg_weight = 1/(total_veg_weight),
         log_estimated_biomass_plot = log(biomass),
         log_CN_ratio = log(CN_ratio))

# T-tests and ANOVAs
run_anova <- function(df, y_var) {
  res.aov <- aov(as.formula(paste(y_var, "~ Animal + treatment")), data = df)
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  return(res.aov)
}

#check if an anova with * is better than with +. not the case for any of the dailies and gradients
compare_anova <- function(df, y_var) {
  res.aov1 <- aov(as.formula(paste(y_var, "~ Animal + treatment")), data = df)
  res.aov2 <- aov(as.formula(paste(y_var, "~ Animal * treatment")), data = df)
  print(anova(res.aov1, res.aov2))
}

veg_aov <- run_anova(veg_combined, "log_total_veg_weight")
height_aov <- run_anova(veg_growth, "log_veg_growth")
biomass_aov <- run_anova(veg_combined, "log_estimated_biomass_plot")

daily_veg_aov <- run_anova(daily_veg_combined, "total_veg_weight")
daily_height_aov <- run_anova(daily_veg_growth, "veg_growth")
daily_biomass_aov <- run_anova(daily_veg_combined, "estimated_biomass_plot")
daily_CN_ratio_aov <- run_anova(daily_veg_combined, "CN_ratio")

gradient_veg_aov <- run_anova(gradient_veg_combined, "log_total_veg_weight")
gradient_height_aov <- run_anova(gradient_veg_growth, "log_veg_growth")
gradient_biomass_aov <- run_anova(gradient_veg_combined, "log_estimated_biomass_plot")
gradient_CN_ratio_aov <- run_anova(gradient_veg_combined, "log_CN_ratio")

# Non-parametric tests due to non-normality of CN data, transforming the data did not help
CN_ratio_art <- art(CN_ratio ~ Animal * treatment, data = veg_combined)
anova(CN_ratio_art)
summary(CN_ratio_art) # good enough, they are all 0.00

# Plot residuals
plot_residuals <- function(aov_model) {
  plot(aov_model, 2)
  shapiro.test(residuals(aov_model))
}

plot_residuals(gradient_CN_ratio_aov)
compare_anova(veg_combined_gradient, "log_CN_ratio")

#summary stats
daily_veg_combined <- daily_veg_combined %>% left_join(daily_veg_growth, by=intersect(names(daily_veg_combined), names(daily_veg_growth)))
daily_veg_combined <- daily_veg_combined %>% left_join(final_species_data, by=intersect(names(daily_veg_combined), names(final_species_data)), multiple = "last")
daily_summary_stats <- daily_veg_combined %>%
  group_by(Animal, treatment) %>%
  dplyr::summarise(
    total_veg_weight_mean = mean(total_veg_weight, na.rm = TRUE),
    total_veg_weight_median = median(total_veg_weight, na.rm = TRUE),
    total_veg_weight_sd = sd(total_veg_weight, na.rm = TRUE),
    total_veg_weight_se = std.error(total_veg_weight),
    veg_growth_mean = mean(veg_growth, na.rm = TRUE),
    veg_growth_median = median(veg_growth, na.rm = TRUE),
    veg_growth_sd = sd(veg_growth, na.rm = TRUE),
    veg_growth_se = std.error(veg_growth),
    estimated_biomass_plot_mean = mean(biomass, na.rm = TRUE),
    estimated_biomass_plot_median = median(biomass, na.rm = TRUE),
    estimated_biomass_plot_sd = sd(biomass, na.rm = TRUE),
    estimated_biomass_plot_se = std.error(biomass),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    species_count_mean = mean(species_count, na.rm = TRUE),
    species_count_median = median(species_count, na.rm = TRUE),
    species_count_sd = sd(species_count, na.rm = TRUE),
    species_count_se = std.error(species_count)
  )

gradient_veg_combined <- gradient_veg_combined %>% left_join(gradient_veg_growth, by=intersect(names(gradient_veg_combined), names(gradient_veg_growth)))
gradient_veg_combined <- gradient_veg_combined %>% left_join(final_species_data, by=intersect(names(gradient_veg_combined), names(final_species_data)), multiple = "last")
gradient_summary_stats <- gradient_veg_combined %>%
  group_by(Animal, treatment) %>%
  dplyr::summarise(
    total_veg_weight_mean = mean(total_veg_weight, na.rm = TRUE),
    total_veg_weight_median = median(total_veg_weight, na.rm = TRUE),
    total_veg_weight_sd = sd(total_veg_weight, na.rm = TRUE),
    total_veg_weight_se = std.error(total_veg_weight),
    veg_growth_mean = mean(veg_growth, na.rm = TRUE),
    veg_growth_median = median(veg_growth, na.rm = TRUE),
    veg_growth_sd = sd(veg_growth, na.rm = TRUE),
    veg_growth_se = std.error(veg_growth),
    estimated_biomass_plot_mean = mean(biomass, na.rm = TRUE),
    estimated_biomass_plot_median = median(biomass, na.rm = TRUE),
    estimated_biomass_plot_sd = sd(biomass, na.rm = TRUE),
    estimated_biomass_plot_se = std.error(biomass),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    species_count_mean = mean(species_count, na.rm = TRUE),
    species_count_median = median(species_count, na.rm = TRUE),
    species_count_sd = sd(species_count, na.rm = TRUE),
    species_count_se = std.error(species_count)
  )

# Species data analysis # HERE!
species_summary <- species_list %>%
  group_by(plot_id) %>%
  dplyr::summarise(species_count = n_distinct(species_list)) %>% 
  ungroup()

species_data <- species_list %>%
  group_by(plot_id, veg_class) %>%
  dplyr::summarise(species_per_vegclass = n_distinct(species_list)) %>% 
  ungroup() %>% 
  mutate(veg_category = veg_class) %>% 
  mutate(veg_cat = case_when(veg_category == "Bryophytes" ~ "Bryophytes",
                             veg_category == "Graminoids" ~ "Graminoids", 
                             TRUE ~ "Forbs")) %>% 
  select(!veg_category) %>% 
  mutate(veg_category = veg_cat)

veg_forbs <- species_data %>% 
  filter(veg_cat == "Forbs")

veg_non_forbs <- species_data %>% filter(veg_cat != "Forbs")

veg_merged <- bind_rows(veg_forbs, veg_non_forbs)

final_species_data <- left_join(species_data, species_summary, by = c("plot_id", "veg_cat"))

final_species_data <- final_species_data %>%
  mutate(
    Animal = case_when(
      grepl("^C", plot_id) ~ "Cow",
      grepl("^H", plot_id) ~ "Horse",
      TRUE ~ NA_character_  # Assign NA if neither
    ),
    Campaign = case_when(
      grepl("G", plot_id) ~ "Gradient",
      grepl("D", plot_id) ~ "Daily",
      TRUE ~ NA_character_  # Assign NA if neither
    ),
    treatment = case_when(
      grepl("_F$", plot_id) ~ "Fresh",  # Ends with "_F"
      grepl("_C$", plot_id) ~ "Control",  # Ends with "_C"
      TRUE ~ NA_character_  # Assign NA if neither
    ))

# create df with species    
comb_data <- final_species_data %>%
  left_join(species_list, by = c("plot_id", "veg_class"))

# summarize to create species_list
species_summary <- comb_data %>%
  group_by(plot_id, veg_class) %>%
  dplyr::summarize(
    species_list = list(unique(species_list)),
    .groups = "drop"
  ) %>% 
  ungroup()

# Join back to retain the original data frame structure, now with species added as a list per row
final_species_data <- final_species_data %>%
  left_join(species_summary, by = c("plot_id", "veg_class"))

# add veg weight per veg clas
final_species_data <- final_species_data %>%
  left_join(veg_new, by = c("plot_id", "veg_class"))

#join dfs to add weight per species
forb_weights <- veg_new %>%
  group_by(plot_id) %>%
  dplyr::summarise(
    forb_weight = sum(dry_weight[!(veg_class %in% c("Bryophytes", "Graminoids"))], na.rm = TRUE),
    .groups = 'drop'
  )

veg_weight <- veg_new %>%
  left_join(forb_weights, by = "plot_id") %>%
  mutate(
    weight_per_class = case_when(
      veg_class %in% c("Bryophytes", "Graminoids") ~ dry_weight,
      TRUE ~ forb_weight
    )
  ) %>%
  select(-forb_weight)

veg_weight <- veg_weight %>% 
  mutate(
    veg_class = ifelse(veg_class %in% c("Bryophytes", "Graminoids"), veg_class, "Forbs")
  ) 

forb_unique <- veg_weight %>% 
  filter(veg_class == "Forbs") %>% 
  distinct(plot_id, .keep_all = TRUE)

veg_weight <- veg_weight %>% 
  filter(veg_class != "Forbs") %>% 
  bind_rows(forb_unique) %>% 
  select(-dry_weight)

veg_weight <- veg_weight %>% 
  mutate(harvested_area_m2 = harvested_area / 10000) %>% 
  mutate(adjweight_per_class = weight_per_class/harvested_area_m2) %>% 
  mutate(Animal = case_when(grepl("^C", plot_id) ~ "Cow", grepl("^H", plot_id) ~ "Horse")) %>% 
  mutate(Animal = factor(Animal))


# T-tests and ANOVAs for the species data are no good, the data distribution is unsuitable
species_aov <- run_anova(final_species_data, "species_count")
species_model <- lmer(species_per_vegclass ~ Animal * treatment + (1|veg_class), data = final_species_data)
summary(species_model)
shapiro.test(residuals(species_model))

saveRDS(final_species_data, "data/species_data.rds")
saveRDS(dung_data, "data/dung_data.rds")
saveRDS(veg_combined, "data/vegetation_data.rds")



# Create the stacked bar plots per campaign

# Gradient
gradient_stacked_weights <- ggplot(data = veg_weight %>% filter(grepl("G", plot_id)), 
                                   aes(x = plot_id, y = adjweight_per_class, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free") +
  geom_bar_pattern(aes(pattern = veg_class), 
                   position = "stack", 
                   stat = "identity", 
                   pattern_fill = "black",     
                   pattern_density = 0.1, 
                   pattern_spacing = 0.03) +
  ggtitle("Gradient weight per vegetation class") +
  xlab("\nPlot ID") +
  ylab(expression("Adjusted Weight (g/m"^2*")")) +
  theme_minimal() +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +  # No patterns in fill legend
  scale_pattern_manual(name = "Vegetation class",
                       values = c("Bryophytes" = "stripe", "Graminoids" = "crosshatch", "Forbs" = "circle"), # Customize as needed
                       guide = guide_legend(override.aes = list(fill = "transparent", 
                                                                color = "black"))) + # Custom pattern legend
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line())

print(gradient_stacked_weights)

ggsave(filename = "veg_plots/gradientvegweightperclass.jpeg", plot = gradient_stacked_weights, width = 6, height = 4)

## Daily
daily_stacked_weights <- ggplot(data = veg_weight %>% filter(grepl("D", plot_id)), 
                                   aes(x = plot_id, y = adjweight_per_class, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free") +
  geom_bar_pattern(aes(pattern = veg_class), 
                   position = "stack", 
                   stat = "identity", 
                   pattern_fill = "black",     
                   pattern_density = 0.1, 
                   pattern_spacing = 0.03) +
  ggtitle("Daily weight per vegetation class") +
  xlab("\nPlot ID") +
  ylab(expression("Adjusted Weight (g/m"^2*")")) +
  theme_minimal() +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +  # No patterns in fill legend
  scale_pattern_manual(name = "Vegetation class",
                       values = c("Bryophytes" = "stripe", "Graminoids" = "crosshatch", "Forbs" = "circle"), # Customize as needed
                       guide = guide_legend(override.aes = list(fill = "transparent", 
                                                                color = "black"))) + # Custom pattern legend
  scale_y_continuous(limits = c(0, 400), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line())

print(daily_stacked_weights)

ggsave(filename = "veg_plots/dailyvegweightperclass.jpeg", plot = daily_stacked_weights, width = 6, height = 4)


# modelling Lasses way xD
veg_height_m1 <- lm(veg_height_2 ~ Animal * treatment, data = veg_combined)
veg_height_m2 <- lmer(veg_height_2 ~ Animal * treatment + (1|base_code), data = veg_combined)
veg_height_m3 <- glmmTMB(veg_height_2 ~ Animal * treatment, data = veg_combined)
veg_height_m4 <- glmmTMB(veg_height_2 ~ Animal * treatment + (1|base_code), data = veg_combined)
veg_height_m5 <- glmmTMB(veg_height_2 ~ Animal * treatment, family = poisson, data = veg_combined)
summary(veg_height_daily_m1)
Anova(veg_height_daily_m1)
# Model validation
simulationOutput <- simulateResiduals(fittedModel = veg_height_gradient_m1, n = 1000)
testDispersion(simulationOutput)

plot(simulationOutput)
plotResiduals(simulationOutput, form = gradient_veg_combined$Animal)

plotResiduals(simulationOutput, form = gradient_veg_combined$treatment)
# Post-hoc-test
test <- emmeans(veg_height_gradient_m1, ~ treatment|Animal)
contrast(test, method = "pairwise") %>% as.data.frame()

# make it into a function
run_model <- function(dataset, model) {
  #print(summary(model))
  print(Anova(model))
  simuOutput <- simulateResiduals(fittedModel = model, n = 1000)
  plot(simuOutput)
  plotResiduals(simuOutput, form = dataset$Animal)
  plotResiduals(simuOutput, form = dataset$treatment)
  plotResiduals(simuOutput, form = dataset$Campaign)
  test <- emmeans(model, ~ Campaign|treatment|Animal)
  contrast(test, method = "pairwise") %>% as.data.frame()
}

#dailies
veg_weight_daily_m1 <- glmmTMB(total_veg_weight ~ Animal * treatment + (1|base_code), data = daily_veg_combined)
veg_height_daily_m1 <- glmmTMB(veg_height_2 ~ Animal * treatment + (1|base_code), data = daily_veg_combined)
veg_growth_daily_m1 <- glmmTMB(veg_growth ~ Animal * treatment + (1|base_code), data = daily_veg_growth)
veg_biomass_daily_m1 <- glmmTMB(biomass ~ Animal * treatment + (1|base_code), data = daily_veg_combined)
veg_cn_daily_m1 <- glmmTMB(CN_ratio ~ Animal * treatment + (1|base_code), data = daily_veg_combined)

#gradients
veg_weight_gradient_m1 <- glmmTMB(total_veg_weight ~ Animal * treatment + (1|base_code), data = gradient_veg_combined)
veg_height_gradient_m1 <- glmmTMB(veg_height_2 ~ Animal * treatment + (1|base_code), data = gradient_veg_combined)
veg_growth_gradient_m1 <- glmmTMB(veg_growth ~ Animal * treatment + (1|base_code), data = gradient_veg_growth)
veg_biomass_gradient_m1 <- glmmTMB(biomass ~ Animal * treatment + (1|base_code), data = gradient_veg_combined)
veg_cn_gradient_m1 <- glmmTMB(CN_ratio ~ Animal * treatment + (1|base_code), data = gradient_veg_combined)

#both campaigns together
veg_weight_m1 <- glmmTMB(total_veg_weight ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)
veg_height_m1 <- glmmTMB(veg_height_2 ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)
veg_growth_m1 <- glmmTMB(height_value ~ Animal * treatment * Campaign + (1|base_code), data = veg_growth)
veg_biomass_m1 <- glmmTMB(biomass ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)
veg_cn_m1 <- glmmTMB(CN_ratio ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)

run_model(veg_combined, veg_cn_m1)

#dung data
dung_data <- dung_data %>%
  mutate(
    campaign = case_when(
      Kode_1 %in% c("CG dung", "HG dung") ~ "gradient",
      Kode_1 %in% c("HOD dung", "COD dung", "HND dung", "CND dung") ~ "daily"
    ),
    Animal = case_when(
      Kode_1 %in% c("CG dung", "COD dung", "CND dung") ~ "cow",
      Kode_1 %in% c("HG dung", "HOD dung", "HND dung") ~ "horse"
    )
  )

dung_cn_m1 <- glmmTMB(CN_ratio ~ Animal * campaign, data = dung_data)
summary(dung_cn_m1)
Anova(dung_cn_m1)
simulationOutput <- simulateResiduals(fittedModel = dung_cn_m1, n = 1000)
testDispersion(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = dung_data$Animal)
plotResiduals(simulationOutput, form = dung_data$campaign)
test <- emmeans(dung_cn_m1, ~ campaign|Animal)
contrast(test, method = "pairwise") %>% as.data.frame()

# Remove trash and save for combined analysis
clean_veg_data <- veg_combined %>% 
  subset(select = -c(Animal, Animal_treatment, Campaign, frame_area_cm2))

saveRDS(clean_veg_data, file = "data/clean_veg_data.rds")
