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
library(effects)

# Importeffects# Import data
vegdung_lab <- read.csv("data/Plant_Dung_CN_pH_elements.csv", sep = ";", check.names = FALSE)
veg_raw <- read.csv("data/vegetation_data.csv")
fieldwork_data_raw <- read.csv("data/Fieldwork_data_final.csv")
species_list <- read.csv("data/species_lists.csv")
date_data <- readRDS("flux_data/clean_flux_data.rds")

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
  left_join(fieldwork_data %>% dplyr::select(plot_id, dung_area_cm2), by = "plot_id", multiple = "first") %>%
  mutate(
    frame_area_cm2 = 3058.15,
    area_minus_dung = frame_area_cm2 - dung_area_cm2
  ) %>%
  left_join(dplyr::select(veg_new, plot_id, harvested_area), by = "plot_id") %>%
  distinct() %>%
  mutate(biomass = total_veg_weight * area_minus_dung,
         total_veg_weight_gm2 = total_veg_weight * 10000,
         Animal = as.factor(Animal),
         treatment = as.factor(treatment))

# get CN and plot_id from the lab data sheet and merge with the rest
veglab_data <- vegdung_lab %>% 
  dplyr::select(Kode_1, CN_ratio, B, Na, Mg, Al, P, S, K, Ca, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, As, Se, Sr, Mo, Cd, Ba, Tl, Pb) %>% 
  filter(!Kode_1 %in% c("CG dung", "HG dung", "HOD dung", "COD dung", "HND dung", "CND dung"))
colnames(veglab_data)[1] <- c("plot_id")

veg_combined <- veg_combined %>% 
  left_join(veglab_data, by = "plot_id")

# create a df for dung
dung_data <- vegdung_lab %>% 
  dplyr::rename(Nperc = "N%") %>% 
  dplyr::select(Kode_1, CN_ratio, Nperc, pH, B, Na, Mg, Al, P, S, K, Ca, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, As, Se, Sr, Mo, Cd, Ba, Tl, Pb) %>% 
  filter(Kode_1 %in% c("CG dung", "HG dung", "HOD dung", "COD dung", "HND dung", "CND dung"))
colnames(dung_data)[1] <- c("plot_id")

vegdung_data <- vegdung_lab %>% 
  dplyr::select(Kode_1, CN_ratio, "N%", B, Na, Mg, Al, P, S, K, Ca, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, As, Se, Sr, Mo, Cd, Ba, Tl, Pb) %>% 
  filter(!grepl("_C$", Kode_1)) %>%
  dplyr::rename(Nperc = "N%") %>% 
  mutate(cordung = case_when(
    grepl("dung", Kode_1) ~ Kode_1,
    grepl("CG", Kode_1) ~ "CG dung",
    grepl("HG", Kode_1) ~ "HG dung",
    Kode_1 %in% c("CD1_F", "CD2_F") ~ "COD dung",
    Kode_1 %in% c("HD1_F", "HD2_F") ~ "HOD dung",
    Kode_1 %in% c("CD3_F", "CD4_F") ~ "CND dung",
    Kode_1 %in% c("HD3_F", "HD4_F", "HD5_F") ~ "HND dung"
  )) %>% 
  mutate(Animal = case_when(grepl("^C", Kode_1) ~ "Cow", grepl("^H", Kode_1) ~ "Horse"))

colnames(vegdung_data)[1] <- c("plot_id")

#dung data
dung_data <- dung_data %>%
  mutate(
    campaign = case_when(
      plot_id %in% c("CG dung", "HG dung") ~ "gradient",
      plot_id %in% c("HOD dung", "COD dung", "HND dung", "CND dung") ~ "daily"
    ),
    Animal = case_when(
      plot_id %in% c("CG dung", "COD dung", "CND dung") ~ "cow",
      plot_id %in% c("HG dung", "HOD dung", "HND dung") ~ "horse"
    )
  )

# get the dates per base code for days since first and add to veg combined
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

veg_combined <- veg_combined %>% 
  left_join(date_data, by = "base_code")

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

veg_growth <- veg_growth %>% 
  left_join(date_data, by = "base_code")


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

veg_daily <- veg_combined %>% filter(grepl("D.", base_code))
veg_gradient <- veg_combined %>% filter(grepl("G.", base_code))


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
dailyvegbiomassbox <- ggplot(veg_daily, aes(x = Animal, y = total_veg_weight_gm2, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  xlab("\nAnimal") + ylab(expression("Adjusted biomass (g/m"^2*")")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung","Horse control", "Horse dung"))+
  ggtitle("A") +
  scale_y_continuous(limits = c(80, 450), breaks = seq(100, 400, by = 100)) +
  theme(legend.position = "none")
dailyvegbiomassbox

gradientvegbiomassbox <- ggplot(veg_gradient, aes(x = Animal, y = total_veg_weight_gm2, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1)) +
  xlab("\nAnimal") + ylab(NULL) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung","Horse control", "Horse dung"))+
  scale_y_continuous(limits = c(80, 450), breaks = seq(100, 400, by = 100)) +
  ggtitle("B")
gradientvegbiomassbox

vegbiomassbox <- dailyvegbiomassbox + gradientvegbiomassbox
vegbiomassbox
ggsave(filename = "veg_plots/vegbiomassbox.jpeg", plot = vegbiomassbox, width = 6, height = 4)

# Percent stacked barplot for veg weight per species
veg_class_order <- c("Achillea millefolium", "Campanula rotundifolia", "Cerastium fontanum", "Daucus carota", 
                     "Euphrasia stricta", "Galium verum", "Hypericum perforatum", "Jacobaea vulgaris", "Knautia arvensis", "Pilosella officinarum", "Pimpinella saxifraga", 
                     "Plantago lanceolata", "Ranunculus sp.", "Rosa canina", "Rumex acetosa", "Rumex acetosella", "Stellaria graminea", "Dandelions and false dandelions", "Trifolium arvense", 
                     "Trifolium campestre", "Trifolium pratense", "Trifolium repens", "Veronica chamaedrys", "Veronica officinalis", "Vicia sativa")

veg_percent <- veg_new %>% 
  mutate(veg_class = factor(veg_class, levels = veg_class_order)) %>% 
  filter(veg_category == "Forbs")

veg_percent <- veg_percent %>% 
  mutate(across("veg_class", str_replace, "Dandelions and false dandelions", "Taraxacum sp. and false dandelions"))

species_data <- species_list %>%
  group_by(plot_id, veg_class) %>%
  dplyr::summarise(species_per_vegclass = n_distinct(species_list)) %>% 
  ungroup() %>% 
  mutate(veg_category = veg_class) %>% 
  mutate(veg_cat = case_when(veg_category == "Bryophytes" ~ "Bryophytes",
                             veg_category == "Graminoids" ~ "Graminoids", 
                             TRUE ~ "Forbs")) %>% 
  dplyr::select(!veg_category) %>% 
  mutate(veg_category = veg_cat)

lookup_species <- with(species_data, setNames(species_per_vegclass, paste(plot_id, veg_category, sep = "_")))

veg_percent <- veg_percent %>%
  dplyr::mutate(
    species_per_vegclass = lookup_species[paste(plot_id, veg_category, sep = "_")]
  )

percentspecies <- ggplot(veg_percent, aes(x = plot_id, y = dry_weight, fill = factor(veg_class)))  +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap("Campaign", scales = "free", labeller = as_labeller(c(Daily="A", Gradient="B"))) +
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
  ggtitle("Proportions of forb species per campaign")

percentspecies

ggsave(filename = "veg_plots/percentspecies.jpeg", plot = percentspecies, width = 10, height = 8)


# Scatter plot of veg weight vs height
scatterplot <- ggplot(veg_combined, aes(x = biomass, y = veg_height_2, color = treatment)) +
  geom_point() + geom_smooth(method = lm)
print(scatterplot)

# subset by animal
veg_horse <- veg_combined %>% filter(Animal == "Horse")
veg_cow <- veg_combined %>% filter(Animal == "Cow")

# test assumptions of the lms
horse_reg <- lm(total_veg_weight_gm2 ~ veg_height_2, data = veg_horse)
cow_reg <- lm(total_veg_weight_gm2 ~ veg_height_2, data = veg_cow)

plot(daily_reg)

gradient_regression <- ggplot(veg_cow, aes(x = veg_height_2, y = total_veg_weight_gm2, color = treatment)) +
  geom_point() + geom_smooth(method = lm) +
  theme_minimal() +
  stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  xlab("Estimated biomass per plot (g)") +
  ylab("Vegetation height (cm)") +
  ggtitle("A") +
  scale_color_manual(values = c("Dung" = "black", "Control" = "gray")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        legend.position = "none") 

print(gradient_regression)

daily_regression <- ggplot(veg_horse, aes(x = veg_height_2, y = total_veg_weight_gm2, color = treatment)) +
  geom_point() + geom_smooth(method = lm) +
  theme_minimal() +
  stat_poly_eq(use_label(c("eq", "R2", "p"))) +
  xlab("Estimated biomass per plot (g)") +
  ylab("Vegetation height (cm)") +
  ggtitle("Horse") +
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

SjPlot_model_gradient <- glmmTMB(total_veg_weight_gm2 ~ veg_height_2   * treatment, data = veg_gradient)
SjPlot_model_horse <- glmmTMB(total_veg_weight_gm2 ~ veg_height_2   * treatment, data = veg_horse)

SjPlot_model_daily <- glmmTMB(total_veg_weight_gm2 ~ veg_height_2    * treatment, data = veg_daily)
SjPlot_model_cow <- glmmTMB(total_veg_weight_gm2 ~ veg_height_2    * treatment, data = veg_cow)

#Validate modelsAnova(veg_height_daily_m1)
# Model validation
simulationOutput <- simulateResiduals(fittedModel = SjPlot_model_horse, n = 1000)
testDispersion(simulationOutput)

plot(simulationOutput)
plotResiduals(simulationOutput, form = SjPlot_model_horse$treatment)
#plotResiduals(simulationOutput, form = gradient_veg_combined$treatment)
# Post-hoc-test
test <- emmeans(SjPlot_model_horse, ~ treatment)
contrast(test, method = "pairwise") %>% as.data.frame()

# Effect plot with predicted values for gradient model
plot_horse_effect <- plot_model(SjPlot_model_horse, type = "pred", 
                                   terms = c("veg_height_2", "treatment"), 
                                   title = "B") +
  theme_minimal() +  
  scale_color_manual(values = c("Fresh" = "#7F4F24", "Control" = "#A68A64"), labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  xlab("Vegetation height (cm)") +
  ylab(expression("Adjusted biomass (g/m"^2*")")) +
  labs(colour = "Treatment") + 
  scale_y_continuous(limits = c(100, 800)) + 
  # scale_x_continuous(breaks = seq(0, 18, by = 2)) +
  theme(
    panel.grid.major.x = element_blank(),  
    axis.line = element_line(color = "black",),
    axis.title.y = element_blank()) +
  geom_point(aes(x = veg_height_2, y = total_veg_weight_gm2, color = treatment),
             data = veg_horse,
             inherit.aes = FALSE)

plot_horse_effect

# Effect plot with predicted values for daily model 
plot_cow_effect <- plot_model(SjPlot_model_cow, type = "pred", show.values = TRUE,
                                terms = c("veg_height_2", "treatment"), 
                                title = "A") +
  theme_minimal() +  
  scale_color_manual(values = c("Fresh" = "#656D4A", "Control" = "#A4AC86"),
                     labels = c("Control", "Dung")) +
  scale_fill_manual(values = c("Fresh" = "#696969", "Control" = "#696969")) +
  xlab("Vegetation height (cm)") +
  ylab(expression("Adjusted biomass (g/m"^2*")")) +
  labs(colour = "Treatment") +  
  # scale_y_continuous(limits = c(0, 650)) + 
  scale_x_continuous(breaks = seq(0, 12, by = 4)) +
  theme( 
    panel.grid.major.x = element_blank(),  
    axis.line = element_line(color = "black"),
   legend.position = "none") +
  geom_point(aes(x = veg_height_2, y = total_veg_weight_gm2, color = treatment),
                                        data = veg_daily,
                                        inherit.aes = FALSE)
   
plot_cow_effect

combined_model_plot <- plot_cow_effect + plot_horse_effect
combined_model_plot
ggsave("veg_plots/model_biomass_vegheight_animal.jpeg", plot = combined_model_plot, width = 7, height = 4)

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
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung","Horse control", "Horse dung"))+
  ggtitle("B")
gradientcnplot
ggsave(filename = "veg_plots/gradientcn.jpeg", plot = gradientcnplot, width = 6, height = 4)

dailycnplot <- ggplot(veg_daily, aes(x = Animal, y = CN_ratio, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab("C:N ratio") +
  scale_y_continuous(limits = c(12, 38), breaks = seq(12, 38, by = 4)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  ggtitle("A") +
  theme(legend.position = "none")
dailycnplot
ggsave(filename = "veg_plots/dailycn.jpeg", plot = dailycnplot, width = 6, height = 4)

combinedcnplot <- dailycnplot + gradientcnplot
combinedcnplot
ggsave(filename = "veg_plots/combinedcn.jpeg", plot = combinedcnplot, width = 6, height = 4)

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
    veg_growth_mean = mean(height_value, na.rm = TRUE),
    veg_growth_median = median(height_value, na.rm = TRUE),
    veg_growth_sd = sd(height_value, na.rm = TRUE),
    veg_growth_se = std.error(height_value),
    estimated_harvested_biomass_plot_mean = mean(total_veg_weight_gm2, na.rm = TRUE),
    estimated_harvested_biomass_plot_median = median(total_veg_weight_gm2, na.rm = TRUE),
    estimated_harvested_biomass_plot_sd = sd(total_veg_weight_gm2, na.rm = TRUE),
    estimated_harvested_biomass_plot_se = std.error(total_veg_weight_gm2),
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
    veg_growth_mean = mean(height_value, na.rm = TRUE),
    veg_growth_median = median(height_value, na.rm = TRUE),
    veg_growth_sd = sd(height_value, na.rm = TRUE),
    veg_growth_se = std.error(height_value),
    estimated_harvested_biomass_plot_mean = mean(total_veg_weight_gm2, na.rm = TRUE),
    estimated_harvested_biomass_plot_median = median(total_veg_weight_gm2, na.rm = TRUE),
    estimated_harvested_biomass_plot_sd = sd(total_veg_weight_gm2, na.rm = TRUE),
    estimated_harvested_biomass_plot_se = std.error(total_veg_weight_gm2),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    species_count_mean = mean(species_count, na.rm = TRUE),
    species_count_median = median(species_count, na.rm = TRUE),
    species_count_sd = sd(species_count, na.rm = TRUE),
    species_count_se = std.error(species_count)
  )

dung_summary_stats <- dung_data %>%
  group_by(Animal) %>% 
  dplyr::summarise(
    Nperc_mean = mean(Nperc, na.rm = TRUE),
    Nperc_se = std.error(Nperc),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    pH_mean = mean(pH, na.rm = TRUE),
    pH_median = median(pH, na.rm = TRUE),
    pH_sd = sd(pH, na.rm = TRUE),
    pH_se = std.error(pH),
    P_mean = mean(P, na.rm = TRUE),
    P_se = std.error(P),
    K_mean = mean(K, na.rm = TRUE),
    K_se = std.error(K),
    Ca_mean = mean(Ca, na.rm = TRUE),
    Ca_se = std.error(Ca),
    S_mean = mean(S, na.rm = TRUE),
    S_se = std.error(S),
    Mg_mean = mean(Mg, na.rm = TRUE),
    Mg_se = std.error(Mg)
  )

# Species data analysis
species_summary <- species_list %>%
  group_by(plot_id) %>%
  dplyr::summarise(species_count = n_distinct(species_list)) %>% 
  ungroup()


veg_forbs <- species_data %>% 
  filter(veg_cat == "Forbs")

veg_non_forbs <- species_data %>% filter(veg_cat != "Forbs")

veg_merged <- bind_rows(veg_forbs, veg_non_forbs)

final_species_data <- left_join(species_data, species_summary, by = c("plot_id"))

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
  dplyr::select(-forb_weight)

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
  dplyr::select(-dry_weight)

veg_weight <- veg_weight %>% 
  mutate(harvested_area_m2 = harvested_area / 10000) %>% 
  mutate(adjweight_per_class = weight_per_class/harvested_area_m2) %>% 
  mutate(Animal = case_when(grepl("^C", plot_id) ~ "Cow", grepl("^H", plot_id) ~ "Horse")) %>% 
  mutate(Animal = factor(Animal))

# Calculating forb:graminoid ratio here
forb_graminoid_ratio <- veg_weight %>% 
  filter(veg_class %in% c("Forbs", "Graminoids")) %>% # Keep only Forbs and Graminoids
  group_by(plot_id, veg_class) %>%
  summarise(total_weight = sum(adjweight_per_class, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = veg_class, values_from = total_weight, values_fill = 0) %>%
  mutate(forb_graminoid_ratio = Forbs / Graminoids) %>% 
  mutate(Animal = case_when(grepl("^C", plot_id) ~ "Cow", grepl("^H", plot_id) ~ "Horse"),
         treatment = case_when(
           grepl("_F$", plot_id) ~ "Fresh",
           grepl("_C$", plot_id) ~ "Control",
           TRUE ~ NA_character_
         ),
        Campaign = case_when(
          grepl("G", plot_id) ~ "Gradient",
          grepl("D", plot_id) ~ "Daily"
        )) %>% 
  ungroup()
  
# plot
ggplot(filter(forb_graminoid_ratio, Campaign == "Gradient"), aes( y = forb_graminoid_ratio, fill = interaction(Animal, treatment))) +
  geom_boxplot() +
  labs(
    title = "Forb-to-Graminoid Ratio by plot",
    x = "Plot ID",
    y = "Forb:Graminoid ratio"
  ) +
  theme_minimal()

ggplot(filter(forb_graminoid_ratio, Campaign == "Gradient"), aes(x = plot_id, y = forb_graminoid_ratio, fill = interaction(Animal, treatment))) +
  geom_bar(stat = "identity") +
  labs(
    title = "Forb-to-Graminoid Ratio by plot",
    x = "Plot ID",
    y = "Forb:Graminoid ratio"
  ) +
  theme_minimal()

# table
kable(forb_graminoid_ratio, col.names = c("Plot ID", "Forbs", "Graminoids", "Forb-to-Graminoid Ratio"))

# stats
gradient_forb_graminoid <- forb_graminoid_ratio %>% filter(Campaign == "Gradient")
forb_graminoid_model <- glmmTMB(forb_graminoid_ratio ~ treatment * Animal + (1|plot_id), data = gradient_forb_graminoid)
Anova(forb_graminoid_model)
simulationOutput <- simulateResiduals(fittedModel = forb_graminoid_model, n = 1000)
testDispersion(simulationOutput)

plot(simulationOutput)
forb_graminoid_effect <- allEffects(forb_graminoid_model)
plot(forb_graminoid_effect)
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
  facet_wrap("Animal", scales = "free_x") +
  geom_bar_pattern(aes(pattern = veg_class), 
                   position = "stack", 
                   stat = "identity", 
                   pattern_fill = "black",     
                   pattern_density = 0.1, 
                   pattern_spacing = 0.03) +
  ggtitle("Biomass per vegetation class of the long term campaign") +
  xlab("\nPlot ID") +
  ylab(expression("Adjusted biomass (g/m"^2*")")) +
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
  facet_wrap("Animal", scales = "free_x") +
  geom_bar_pattern(aes(pattern = veg_class), 
                   position = "stack", 
                   stat = "identity", 
                   pattern_fill = "black",     
                   pattern_density = 0.1, 
                   pattern_spacing = 0.03) +
  ggtitle("Biomass per vegetation class of the short term campaign") +
  xlab("\nPlot ID") +
  ylab(expression("Adjusted biomass (g/m"^2*")")) +
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
#veg_weight_m1 <- glmmTMB(total_veg_weight ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)
#veg_height_m1 <- glmmTMB(veg_height_2 ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)
veg_growth_m1 <- glmmTMB(height_value ~ Animal * treatment * Campaign + (1|base_code), data = veg_growth)
veg_biomass_m1 <- glmmTMB(total_veg_weight_gm2 ~ treatment * Animal * Campaign + (1|base_code), data = veg_combined)
veg_cn_m1 <- glmmTMB(CN_ratio ~ Animal * treatment * Campaign + (1|base_code), data = veg_combined)

#models with dung age - IMPROVED FOR PAPER
veg_growth_d1 <- glmmTMB(height_value ~ Animal * treatment * dung_age + (1|base_code), data = veg_growth)
veg_cn_d1 <- glmmTMB(CN_ratio ~ Animal * treatment * dung_age + (1|base_code), data = veg_combined)

run_model(veg_growth, veg_growth_d1)

#test model effects
veg_height_effects <- allEffects(veg_height_m1)
veg_growth_effects <- allEffects(veg_growth_m1)
veg_biomass_effects <- allEffects(veg_biomass_m1)
veg_cn_effects <- allEffects(veg_cn_m1)
plot(veg_biomass_effects)


dung_cn_m1 <- glmmTMB(CN_ratio ~ Animal * campaign + (1|plot_id), data = dung_data)
dung_ph_m1 <- glmmTMB(pH ~ Animal * campaign + (1|plot_id), data = dung_data)
summary(dung_ph_m1)
Anova(dung_ph_m1)
simulationOutput <- simulateResiduals(fittedModel = dung_ph_m1, n = 1000)
testDispersion(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = dung_data$Animal)
plotResiduals(simulationOutput, form = dung_data$campaign)
test <- emmeans(dung_cn_m1, ~ campaign|Animal)
contrast(test, method = "pairwise") %>% as.data.frame()

# Remove trash and save for combined analysis
clean_veg_data <- veg_combined %>% 
  subset(dplyr::dplyr::select = -c(Animal, Animal_treatment, Campaign, frame_area_cm2))

#saveRDS(clean_veg_data, file = "data/clean_veg_data.rds")

# Vegetation elemental properties ------------------------------------------------------

veg_phosphorus <- ggplot(veg_combined, aes(x = plot_id, y = P, fill = interaction(treatment, Animal))) +
  geom_bar(stat = "identity") +
  facet_wrap("Animal", scales = "free") +
  xlab("\nPlot ID") + ylab("P (mg/kg)") +
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
  ggtitle("P  graph")
veg_phosphorus

phosphorus <- ggplot(veg_combined, aes(x = plot_id, y = P, colour = interaction(treatment, Animal))) +
  geom_point(stat = "identity") +
  xlab("\nPlot ID") + ylab("P (mg/kg)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_colour_manual(values = c("Control.Cow" = "#A4AC86", 
                               "Fresh.Cow" = "#656D4A", 
                               "Control.Horse" = "#A68A64", 
                               "Fresh.Horse" = "#7F4F24"),
                    labels = c("Cow control", "Cow dung", "Horse control", "Horse dung"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("P  graph")
phosphorus

dung_phosphorus <- ggplot(dung_data, aes(x = plot_id, y = P)) +
  geom_bar(stat = "identity") +
  xlab("\nPlot ID") + ylab("P (mg/kg)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Dung P  graph")
dung_phosphorus

# Create the nested bar chart
level_order <- c("CG dung", "CG1_F", "CG2_F", "CG3_F", "CG4_F", "CG5_F", "COD dung", "CD1_F", "CD2_F", "CND dung", "CD3_F", "CD4_F", "HG dung", "HG1_F", "HG2_F", "HG3_F", "HG4_F", "HG5_F", "HOD dung", "HD1_F", "HD2_F", "HND dung", "HD3_F", "HD4_F", "HD5_F") 
ggplot(vegdung_data, aes(x = plot_id, y = Ca, fill = cordung)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_x_discrete(limits = level_order) +
  #facet_wrap("Animal", scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(veg_growth, file = "plant_data/veg_growth_data.rds")


## Plotting relationship between CN ratio in the dung compared to the CN ratio of the vegetation
vegCN <- vegdung_lab %>% 
  filter(Kode_3 != "Dung") %>% 
  dplyr::select(Kode_1, CN_ratio, "N%") %>% 
  filter(grepl(".F", Kode_1)) %>% 
  mutate(dung_name =  case_when(
    grepl("^CG", Kode_1) ~ "CG",
    grepl("^HG", Kode_1) ~ "HG",
    Kode_1 %in% c("CD1_F", "CD2_F") ~ "COD",
    Kode_1 %in% c("CD3_F", "CD4_F") ~ "CND",
    Kode_1 %in% c("HD1_F", "HD2_F") ~ "HOD",
    Kode_1 %in% c("HD3_F", "HD4_F", "HD5_F") ~ "HND")) %>% 
  rename(base_code = Kode_1,
         N_perc = "N%") %>% 
  group_by(dung_name) %>% 
  summarise(
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    Nperc_mean = mean(N_perc, na.rm = TRUE),
    Nperc_se = std.error(N_perc)) %>% 
  ungroup() 

dung_data <- dung_data %>% 
  mutate(dung_name =
    case_when(
    plot_id %in% "CG dung" ~ "CG",
    plot_id %in% "HG dung" ~ "HG",
    plot_id %in% "HOD dung" ~ "HOD",
    plot_id %in% "COD dung" ~ "COD",
    plot_id %in% "HND dung" ~ "HND",
    plot_id %in% "CND dung" ~ "CND",
  )) %>% 
  dplyr::select(-5:-28)

vegdungCN <- left_join(vegCN, dung_data, by = "dung_name" )


Npercentage <- ggplot(vegdungCN, aes(x = Nperc, y = Nperc_mean, colour = Animal)) +
  geom_point(stat = "identity") +
  xlab("Dung N%") + ylab("Plant N%") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(colour = "Animal") +
  scale_colour_manual(values = c( 
                                 "cow" = "#656D4A", 
                                 "horse" = "#7F4F24"),
                      labels = c("Cow dung", "Horse dung"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Dung and plant N%")
Npercentage
_