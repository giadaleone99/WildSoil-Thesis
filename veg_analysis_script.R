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

# Import data
vegdung_lab <- read.csv("data/vegdung_lab_data.csv")
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
    )
  )

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
                      labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh"))+
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
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetation weight of the long term campaign")
gradientvegweight
ggsave(filename = "veg_plots/gradientvegweight.jpeg", plot = gradientvegweight, width = 6, height = 4)

# stacked bar plot for veg height 1 and 2 per campaign
veg_heightdatadaily <- veg_heightdata %>% filter(grepl("D.", base_code))
veg_heightdatadaily$height_type <- factor(veg_heightdatadaily$height_type, 
                                          levels = c("veg_height_2", "veg_height"))

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
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +  
  scale_pattern_manual(name = "Measurement period", 
                       values = c("veg_height" = "none", "veg_height_2" = "stripe"),
                       labels = c("1", "2"), 
                       guide = guide_legend(override.aes = list(
                         fill = "transparent", 
                         pattern_fill = c("black", "transparent"),
                         pattern = c("none", "stripe"),  # Ensure correct pattern display
                         color = "black"  # Border color for legend boxes
                       ))) + 
                       #guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  theme(legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), expand = c(0, 0)) +  ggtitle("Vegetation height and growth of the short term campaign")

dailyvegheights
ggsave(filename = "veg_plots/veg_growth_daily_stacked.jpeg", plot = dailyvegheights, width = 6, height = 4)


veg_heightdatagradient <- veg_heightdata %>% filter(grepl("G.", base_code))
veg_heightdatagradient$height_type <- factor(veg_heightdatagradient$height_type, 
                                          levels = c("veg_height_2", "veg_height"))


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
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh"),
                    guide = guide_legend(override.aes = list(pattern = "none"))) +  
  scale_pattern_manual(name = "Measurement period", 
                       values = c("veg_height" = "none", "veg_height_2" = "stripe"),
                       labels = c("1", "2"), 
                       guide = guide_legend(override.aes = list(
                         fill = "transparent", 
                         pattern_fill = c("black", "transparent"),
                         pattern = c("none", "stripe"),  
                         color = "black"  
                       ))) + 
  theme(legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10), expand = c(0, 0)) +
  ggtitle("Vegetation height and growth of the long term campaign")


gradientvegheights
ggsave(filename = "veg_plots/veg_growth_gradient_stacked.jpeg", plot = gradientvegheights, width = 6, height = 4)

# Biomass per campaign
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
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh"))+
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
                    labels = c("Cow control", "Cow fresh", "Horse control", "Horse fresh"))+
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetation biomass per plot of the long term campaign")
gradientvegbiomass
ggsave(filename = "veg_plots/gradientvegbiomass.jpeg", plot = gradientvegbiomass, width = 6, height = 4)


veg_gradient$plot_id <- as.factor(veg_gradient$plot_id)

ggplot(veg_gradient, aes(x = plot_id, y = CN_ratio)) +
  geom_bar(stat = "identity") 

barplot(veg_gradient$biomass)


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
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray")) +
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
  scale_color_manual(values = c("Fresh" = "black", "Control" = "gray")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line())

print(daily_regression)

combined_regression <- gradient_regression + daily_regression 
# print(combined_plot)
print(combined_regression)
ggsave(filename = "veg_plots/combined_regression.jpeg", plot = combined_regression, width = 6, height = 4)



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

# Filtering for veg height growth and bar plot
veg_growth <- veg_summary %>%
  filter(!is.na(veg_height)) %>%
  mutate(veg_growth = veg_height_2 - veg_height)

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
  xlab("\nAnimal") + ylab("CN ratio") +
  scale_y_continuous(breaks = seq(12, 38, by = 4)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Treatment") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow fresh", "Cow control","Horse fresh", "Horse control"))+
  ggtitle("Long term campaign") +
  theme(legend.position = "none")
gradientcnplot
ggsave(filename = "veg_plots/gradientcn.jpeg", plot = gradientcnplot, width = 6, height = 4)

dailycnplot <- ggplot(veg_daily, aes(x = Animal, y = CN_ratio, fill = interaction(treatment, Animal))) +
  geom_boxplot(position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1))+
  xlab("\nAnimal") + ylab(NULL) +
  scale_y_continuous(breaks = seq(12, 38, by = 4)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text.x = element_text(size = 12)) +
  labs(fill = "Treatment") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64"),
                    labels = c("Cow fresh", "Cow control", "Horse fresh", "Horse control"))+
  ggtitle("Short term campaign")
dailycnplot
ggsave(filename = "veg_plots/dailycn.jpeg", plot = dailycnplot, width = 6, height = 4)
combinedcnplot <- gradientcnplot + dailycnplot
combinedcnplot
ggsave(filename = "veg_plots/combinedcn.jpeg", plot = combinedcnplot, width = 8, height = 4)

### STATISTICS ------------------------------------------------------------------------
#make subsets for the campaigns
daily_veg_combined <- veg_combined %>% 
  filter(Campaign == "Daily")
daily_veg_growth <- veg_growth %>% 
  filter(Campaign == "Daily")

gradient_veg_combined <- veg_combined %>% 
  filter(Campaign == "Gradient")
gradient_veg_growth <- veg_growth %>% 
  filter(Campaign == "Gradient")

#residuals are not normally distributed so log transform the data
veg_growth <- veg_growth %>% 
  mutate(log_veg_growth = log(2+veg_growth))
veg_combined <- veg_combined %>% 
  mutate(log_total_veg_weight = log(total_veg_weight),
         log_estimated_biomass_plot = log(estimated_biomass_plot))
gradient_veg_growth <- gradient_veg_growth %>% 
  mutate(log_veg_growth = log(2+veg_growth))
gradient_veg_combined <- gradient_veg_combined %>% 
  mutate(log_total_veg_weight = 1/(total_veg_weight),
         log_estimated_biomass_plot = log(estimated_biomass_plot),
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
    estimated_biomass_plot_mean = mean(estimated_biomass_plot, na.rm = TRUE),
    estimated_biomass_plot_median = median(estimated_biomass_plot, na.rm = TRUE),
    estimated_biomass_plot_sd = sd(estimated_biomass_plot, na.rm = TRUE),
    estimated_biomass_plot_se = std.error(estimated_biomass_plot),
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
    estimated_biomass_plot_mean = mean(estimated_biomass_plot, na.rm = TRUE),
    estimated_biomass_plot_median = median(estimated_biomass_plot, na.rm = TRUE),
    estimated_biomass_plot_sd = sd(estimated_biomass_plot, na.rm = TRUE),
    estimated_biomass_plot_se = std.error(estimated_biomass_plot),
    CN_mean = mean(CN_ratio, na.rm = TRUE),
    CN_median = median(CN_ratio, na.rm = TRUE),
    CN_sd = sd(CN_ratio, na.rm = TRUE),
    CN_se = std.error(CN_ratio),
    species_count_mean = mean(species_count, na.rm = TRUE),
    species_count_median = median(species_count, na.rm = TRUE),
    species_count_sd = sd(species_count, na.rm = TRUE),
    species_count_se = std.error(species_count)
  )

# Species data analysis
species_summary <- species_list %>%
  group_by(plot_id) %>%
  dplyr::summarise(species_count = n_distinct(species_list))

species_data <- species_list %>%
  group_by(plot_id, veg_class) %>%
  dplyr::summarise(species_per_vegclass = n_distinct(species_list))

final_species_data <- left_join(species_data, species_summary, by = "plot_id")

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
  )

# Join back to retain the original data frame structure, now with species added as a list per row
final_species_data <- final_species_data %>%
  left_join(species_summary, by = c("plot_id", "veg_class"))

# add veg weight per veg clas
final_species_data <- final_species_data %>%
  left_join(veg_new, by = c("plot_id", "veg_class"))

#join dfs to add weight per species
veg_weight <- veg_new %>% 
  group_by(plot_id) %>% 
  mutate(
    weight_per_class = case_when(
      veg_class == "Bryophytes" ~ dry_weight,
      veg_class == "Graminoids" ~ dry_weight,
      TRUE ~ sum(dry_weight[veg_class != "Bryophytes" & veg_class != "Graminoids"])
    )
  )
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

# T-tests and ANOVAs for the species data are no good, the data distribution is unsuitable
species_aov <- run_anova(final_species_data, "species_count")
species_model <- lmer(species_per_vegclass ~ Animal * treatment + (1|veg_class), data = final_species_data)
summary(species_model)
shapiro.test(residuals(species_model))

saveRDS(final_species_data, "data/species_data.rds")
saveRDS(dung_data, "data/dung_data.rds")
saveRDS(veg_combined, "data/vegetation_data.rds")

# Create the stacked bar plots per campaign
gradient_stacked_weights <- ggplot(data = veg_weight %>% filter(grepl("G", plot_id)), 
                     aes(x = plot_id, y = weight_per_class, fill = veg_class)) +
  geom_bar(stat = "identity") +
  ggtitle("Gradient weight per vegetation class") +
  xlab("Plot ID") +
  ylab("Weight (g)") +
  scale_fill_brewer(palette = "Set3", name = "Vegetation type") +  # Set the legend title
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))
        #legend.position = "none")  # Hide legend for this plot
print(gradient_stacked_weights)
daily_stacked_weights <- ggplot(data = veg_weight %>% filter(grepl("D", plot_id)), 
                                   aes(x = plot_id, y = weight_per_class, fill = veg_class)) +
  geom_bar(stat = "identity") +
  ggtitle("Daily weight per vegetation class") +
  xlab("Plot ID") +
  ylab("Weight (g)") +
  scale_fill_brewer(palette = "Set3", name = "Vegetation type") +  # Set the legend title
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))
print(daily_stacked_weights)


# modelling Lasses way
veg_height_m1 <- lm(veg_height_2 ~ Animal * treatment, data = veg_combined)
veg_height_m2 <- lmer(veg_height_2 ~ Animal * treatment + (1|base_code), data = veg_combined)
veg_height_m3 <- glmmTMB(veg_height_2 ~ Animal * treatment, data = veg_combined)
veg_height_m4 <- glmmTMB(veg_height_2 ~ Animal * treatment + (1|base_code), data = veg_combined)
veg_height_m5 <- glmmTMB(veg_height_2 ~ Animal * treatment, family = poisson, data = veg_combined)
summary(veg_height_daily_m1)
Anova(veg_height_daily_m1)
# Model validation
simulationOutput <- simulateResiduals(fittedModel = veg_height_daily_m1, n = 1000)
testDispersion(simulationOutput)

plot(simulationOutput)
plotResiduals(simulationOutput, form = daily_veg_combined$Animal)

plotResiduals(simulationOutput, form = daily_veg_combined$treatment)

# Post-hoc-test
test <- emmeans(veg_height_daily_m1, ~ treatment|Animal)
contrast(test, method = "pairwise") %>% as.data.frame()

#per campaign
veg_height_gradient_m1 <- glmmTMB(veg_height_2 ~ Animal * treatment + (1|base_code), data = gradient_veg_combined)
veg_height_daily_m1 <- glmmTMB(veg_height_2 ~ Animal * treatment + (1|base_code), data = daily_veg_combined)
