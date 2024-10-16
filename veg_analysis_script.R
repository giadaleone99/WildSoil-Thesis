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

# Import data
veg_raw <- read.csv("data/vegetation_data.csv")
fieldwork_data_raw <- read.csv("data/Fieldwork_data_final.csv")
species_list <- read.csv("data/species_lists.csv")

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
    Animal_treatment = paste(Animal, treatment)
  )

# Fieldwork data manipulation and joining with vegetation data
fieldwork_data <- fieldwork_data_raw %>%
  mutate(plot_id = str_remove(Plot_ID, "_[^_]+$"))

veg_combined <- veg_summary %>%
  left_join(fieldwork_data %>% select(plot_id, dung_area_cm2), by = "plot_id") %>%
  mutate(
    frame_area_cm2 = 3058.15,
    area_minus_dung = frame_area_cm2 - dung_area_cm2
  ) %>%
  left_join(select(veg_new, plot_id, harvested_area), by = "plot_id") %>%
  distinct() %>%
  mutate(estimated_biomass_plot = total_veg_weight * area_minus_dung)

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



# Scatter plot of veg weight vs height
scatterplot <- ggplot(veg_combined, aes(x = total_veg_weight, y = veg_height_2, color = treatment)) +
  geom_point() + geom_smooth(method = lm)
ggsave("veg_plots/scatterplot_veg_weight_height.jpeg", plot = scatterplot, width = 6, height = 4)

# Filtering for veg height growth and bar plot
veg_growth <- veg_summary %>%
  filter(!is.na(veg_height)) %>%
  mutate(veg_growth = veg_height_2 - veg_height)

save_plot(veg_growth, "veg_growth", "veg_plots/veg_growth.jpeg", "Veg Growth")

# Bar plot for biomass
save_plot(veg_combined, "estimated_biomass_plot", "veg_plots/estimated_biomass_plot.jpeg", "Estimated Biomass")



# Loop for plotting based on base_code
plot_by_code <- function(df, y_var, folder, y_label) {
  for (code in unique(df$base_code)) {
    df_subset <- df %>% filter(base_code == code)
    p <- ggplot(df_subset, aes(x = plot_id, y = !!sym(y_var), fill = plot_id)) +
      geom_bar(stat = "identity") +
      ggtitle(paste("Bar Plot for", code)) +
      theme_minimal()
    ggsave(filename = paste0("veg_plots/", folder, "/", code, "_plot.jpeg"), plot = p, width = 6, height = 4)
  }
}

plot_by_code(veg_combined, "total_veg_weight", "veg_weight", "Veg Weight")
plot_by_code(veg_growth, "veg_growth", "veg_height", "Veg Growth")

# T-tests and ANOVAs
run_anova <- function(df, y_var) {
  res.aov <- aov(as.formula(paste(y_var, "~ Animal_treatment")), data = df)
  print(summary(res.aov))
  print(TukeyHSD(res.aov))
  return(res.aov)
}

veg_aov <- run_anova(veg_combined, "total_veg_weight")
height_aov <- run_anova(veg_growth, "veg_growth")

# Non-parametric tests due to non-normality
kruskal.test(veg_growth ~ Animal_treatment, data = veg_growth)
bartlett.test(veg_growth ~ Animal_treatment, data = veg_growth)

# Plot residuals
plot_residuals <- function(aov_model) {
  plot(aov_model, 2)
  shapiro.test(residuals(aov_model))
}

plot_residuals(height_aov)

# Species data analysis
species_summary <- species_list %>%
  group_by(plot_id) %>%
  summarise(species_count = n_distinct(species_list))

species_data <- species_list %>%
  group_by(plot_id, veg_class) %>%
  summarise(species_per_vegclass = n_distinct(species_list))

final_species_data <- left_join(species_data, species_summary, by = "plot_id")

final_species_data <- final_species_data %>% 
  mutate(treatment = case_when(
    grepl("_F$", plot_id) ~ "Fresh",  # Ends with "_F"
    grepl("_C$", plot_id) ~ "Control",  # Ends with "_C"
    TRUE ~ NA_character_  # Assign NA if neither
  ))

species_plot <- ggplot(final_species_data, aes(x = plot_id, y = species_per_vegclass, fill = veg_class)) +
  geom_bar(stat = "identity") +  # Use stat = "identity" to plot the species counts
  xlab("Plot ID") +
  ylab("Species Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

# Create the Fresh stacked bar plot
fresh_plot <- ggplot(data = final_species_data %>% filter(treatment == "Fresh"), 
                     aes(x = plot_id, y = species_per_vegclass, fill = veg_class)) +
  geom_bar(stat = "identity") +
  ggtitle("Species Count - Fresh Treatment") +
  xlab("Plot ID") +
  ylab("Species Count") +
  scale_fill_brewer(palette = "Set3", name = "Vegetation type") +  # Set the legend title
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  # Hide legend for this plot

# Create the Control stacked bar plot without the legend
control_plot <- ggplot(data = final_species_data %>% filter(treatment == "Control"), 
                       aes(x = plot_id, y = species_per_vegclass, fill = veg_class)) +
  geom_bar(stat = "identity") +
  ggtitle("Species Count - Control Treatment") +
  xlab("Plot ID") +
  ylab("Species Count") +
  scale_fill_brewer(palette = "Set3", name = "Vegetation type") +  # Set the legend title
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

control_plot
# Print the combined plot

combined_plot <- fresh_plot + control_plot 
# print(combined_plot)
print(combined_plot)


# Save the final plot to the specified directory
ggsave(plot = combined_plot, filename = "veg_plots/combined_species_count_plot.jpeg",
       width = 10, 
       height = 6, 
       dpi = 300)
