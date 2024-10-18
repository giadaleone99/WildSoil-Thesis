### Soil data analysis script 

library(ggplot2)
library(patchwork)

soil_data_raw <- read.csv("data/soil_data_raw.csv")

# create subdatasets for the campaigns

gradient_soil <- soil_data_raw %>% 
  filter(grepl("G.", base_code)) %>% 
  mutate(Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),)
  
gradient_soil$sample_type <- factor(gradient_soil$sample_type, levels = c("Dung soil", "Fresh", "Control"))

daily_soil <- soil_data_raw %>% 
  filter(grepl("D.", base_code)) %>% 
  mutate(Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),)

daily_soil$sample_type <- factor(daily_soil$sample_type, levels = c("Dung soil", "Fresh", "Control"))

summary(gradient_soil)

# create some plots 

gradientphplot <- ggplot(gradient_soil, aes(x = Animal, y = pH, fill = interaction(sample_type, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  xlab("\nAnimal") + ylab("pH") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung soil", "Cow fresh", "Cow control", "Horse dung soil", "Horse fresh", "Horse control"))+
  expand_limits(y = 6) +
  ggtitle("Soil pH of the long term campaign") +
  theme(legend.position = "none")
gradientphplot
#ggsave(filename = "soil_plots/gradientph.jpeg", plot = gradientphplot, width = 6, height = 4)

dailyphplot <- ggplot(daily_soil, aes(x = Animal, y = pH, fill = interaction(sample_type, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  xlab("\nAnimal") + ylab("pH") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung soil", "Cow fresh", "Cow control", "Horse dung soil", "Horse fresh", "Horse control"))+
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) +
  ggtitle("Soil pH of the short term campaign")
dailyphplot
ggsave(filename = "soil_plots/dailyph.jpeg", plot = dailyphplot, width = 6, height = 4)

gradientpo4plot <- ggplot(gradient_soil, aes(x = Animal, y = PO4.P, fill = interaction(sample_type, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  xlab("\nAnimal") + ylab("PO4 P (mg/kg)") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung soil", "Cow fresh", "Cow control", "Horse dung soil", "Horse fresh", "Horse control"))+
  scale_y_continuous(breaks = seq(50, 300, by = 50)) +
  ggtitle("Plant available P of the long term campaign") +
  theme(legend.position = "none")
gradientpo4plot
#ggsave(filename = "soil_plots/gradientpo4.jpeg", plot = gradientpo4plot, width = 6, height = 4)

dailypo4plot <- ggplot(daily_soil, aes(x = Animal, y = PO4.P, fill = interaction(sample_type, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  xlab("\nAnimal") + ylab("PO4 P (mg/kg)") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung soil", "Cow fresh", "Cow control", "Horse dung soil", "Horse fresh", "Horse control"))+
  scale_y_continuous(breaks = seq(50, 370, by = 50)) +
  ggtitle("Plant avaiable P of the short term campaign")
dailypo4plot
ggsave(filename = "soil_plots/dailypo4.jpeg", plot = dailypo4plot, width = 6, height = 4)

gradientcnplot <- ggplot(gradient_soil, aes(x = Animal, y = CN_ratio, fill = interaction(sample_type, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  xlab("\nAnimal") + ylab("CN ratio") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung soil", "Cow fresh", "Cow control", "Horse dung soil", "Horse fresh", "Horse control"))+
  scale_y_continuous(breaks = seq(11, 16, by = 1)) +
  ggtitle("CN ratio of the long term campaign") +
  theme(legend.position = "none")
gradientcnplot
#ggsave(filename = "soil_plots/gradientcn.jpeg", plot = gradientcnplot, width = 6, height = 4)

dailycnplot <- ggplot(daily_soil, aes(x = Animal, y = CN_ratio, fill = interaction(sample_type, Animal))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8))+
  xlab("\nAnimal") + ylab("CN ratio") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line()) +
  labs(fill = "Plot type") +
  scale_fill_manual(values = c("Fresh.Cow" = "#656D4A",
                               "Control.Cow" = "#A4AC86",
                               "Dung soil.Cow"= "#333D29", 
                               "Fresh.Horse" = "#7F4F24",
                               "Control.Horse" = "#A68A64", 
                               "Dung soil.Horse" = "#582F0E"),
                    labels = c("Cow dung soil", "Cow fresh", "Cow control", "Horse dung soil", "Horse fresh", "Horse control"))+
  scale_y_continuous(breaks = seq(11, 15, by = 1)) +
  ggtitle("CN ratio of the short term campaign")
dailycnplot
ggsave(filename = "soil_plots/dailycn.jpeg", plot = dailycnplot, width = 6, height = 4)

combinedphplot <- gradientphplot + dailyphplot
combinedphplot
ggsave(filename = "soil_plots/combinedph.jpeg", plot = combinedphplot, width = 6, height = 4)

combinedpo4plot <- gradientpo4plot + dailypo4plot
combinedpo4plot
ggsave(filename = "soil_plots/combinedpo4.jpeg", plot = combinedpo4plot, width = 6, height = 4)

combinedcnplot <- gradientcnplot + dailycnplot
combinedcnplot
ggsave(filename = "soil_plots/combinedcn.jpeg", plot = combinedcnplot, width = 6, height = 4)
