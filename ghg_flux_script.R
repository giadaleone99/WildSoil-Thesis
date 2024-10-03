# script to process ghg flux data

library(stringr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(lubridate)

# Import data 
flux_data_raw <- read.csv("flux_data/combined_data.csv")

flux_data <- flux_data_raw %>%
  mutate(
    base_code = NA,
    treatment = NA,
    NEERE = NA,
    date = NA
  )

# Split the uniqueID into its components
split_values <- str_split_fixed(flux_data_raw$UniqueID, "_", 4)

# Use a loop to fill in the new columns with the split components
for (i in 1:nrow(flux_data_raw)) {
  flux_data$base_code[i] <- split_values[i, 1]  # First component goes into base_code
  flux_data$treatment[i] <- split_values[i, 2]  # Second component goes into treatment
  flux_data$NEERE[i] <- split_values[i, 3]  
  flux_data$date[i] <- split_values[i, 4] 
}

flux_data <- flux_data %>% relocate(date, .after=UniqueID)
flux_data <- flux_data %>% relocate(NEERE, .after=UniqueID)
flux_data <- flux_data %>% relocate(treatment, .after=UniqueID)
flux_data <- flux_data %>% relocate(base_code, .after=UniqueID)

flux_data <- flux_data %>% 
  mutate(
    Animal = case_when(
      grepl("^C", base_code) ~ "Cow",
      grepl("^H", base_code) ~ "Horse"
    )
  )

flux_data <- flux_data %>% 
  mutate(
    Campaign = case_when(
      grepl("G", base_code) ~ "Gradient",
      grepl("D", base_code) ~ "Daily",
      grepl("^PIT", base_code) ~ "PIT"
    )
  )

flux_data <- flux_data %>% 
  mutate(plotID = paste(base_code, treatment, sep = "_"),
         longdate = dmy(paste0(date, " 2024")))

flux_data <- flux_data %>% 
  mutate(plotNEERE = paste(plotID, NEERE, sep = "_"))

        

# Calculating photosynthesis
flux_data <- flux_data %>% 
  mutate(photosynthesis = NA)
  
for(i in 1:(nrow(flux_data) - 1)) {  # Loop until n-1 to avoid going out of bounds
  if (flux_data$NEERE[i] == "NEE" && flux_data$NEERE[i+1] == "RE") {
    if (flux_data$gastype[i] == "CO2") {
      # Subtract the best.flux values from NEE and RE
      flux_data$photosynthesis[i] <-  flux_data$best.flux[i+1] - flux_data$best.flux[i]
    }
  }
}


gradients <-  flux_data %>%  filter(Campaign == "Gradient", Animal == "Cow", gastype == "CH4")
dailies <- flux_data %>%  filter(Campaign == "Daily", Animal == "Horse", gastype == "CH4")

scatterplot <- ggplot(gradients, aes(x = UniqueID, y = best.flux))+
  geom_point(aes(color = treatment)) +
  geom_smooth(method = lm) + 
  xlab("plot ID") +
  ylab("gas type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

scatterplot

scatterplot1 <- ggplot(dailies, aes(x = UniqueID, y = best.flux))+
  geom_point(aes(color = treatment)) +
  geom_smooth(method = lm) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
 
scatterplot1

zoomdaily <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul"), Animal == "Horse")

scatterplot2 <- ggplot(zoomdaily, aes(x = UniqueID, y = best.flux))+
  geom_point(aes(color = treatment)) +
  geom_smooth(method = lm) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
scatterplot2

zoomdaily <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul"), Animal %in% c("Horse", "Cow"))

ggplot(zoomdaily, aes(x = date, y = best.flux, color = treatment)) +
  geom_point() +
  facet_wrap(~ Animal, scales = "free_y") +
  labs(x = "Date", y = "CO2 Flux")

# Horse daily plots CO2
HD_CO2 <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Horse")

HD_CO2_plot <- ggplot(HD_CO2, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free_x")+
  labs(title = "Daily Horse CO2",
       x = "Date", y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

ggsave(filename = "plots/HD_CO2.jpeg", plot = HD_CO2_plot, width = 10, height = 8)

# Horse daily plots CH4
HD_CH4 <- flux_data %>% filter(Campaign == "Daily", gastype == "CH4", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Horse")

HD_CH4_plot <- ggplot(HD_CH4, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Daily Horse CH4",
       x = "Date", y = expression("nmol CH4 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

ggsave(filename = "plots/HD_CH4.jpeg", plot = HD_CH4_plot, width = 10, height = 8)

# Horse Daily N2O plots
HD_N2O <- flux_data %>% filter(Campaign == "Daily", gastype == "N2O", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Horse")

HD_N2O_plot <- ggplot(HD_N2O, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free_x")+
  labs(title = "Daily Horse N2O",
       x = "Date", y = expression("nmol N2O m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

ggsave(filename = "plots/HD_N2O.jpeg", plot = HD_N2O_plot, width = 10, height = 8)

# Cow daily plots CO2
CD_CO2 <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Cow")

CD_CO2_plot <- ggplot(CD_CO2, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free_x")+
  labs(title = "Daily Cow CO2",
       x = "Date", y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = "bottom",
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

ggsave(filename = "plots/CD_CO2.jpeg", plot = CD_CO2_plot, width = 10, height = 8)

# Cow daily plots CH4
CD_CH4 <- flux_data %>% filter(Campaign == "Daily", gastype == "CH4", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Cow")

CD_CH4_plot <- ggplot(CD_CH4, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Daily Cow CH4",
       x = "Date", y = expression("nmol CH4 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = "bottom",
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

ggsave(filename = "plots/CD_CH4.jpeg", plot = CD_CH4_plot, width = 10, height = 8)

# Cow daily plots N2O
CD_N2O <- flux_data %>% filter(Campaign == "Daily", gastype == "N2O", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Cow")

CD_N2O_plot <- ggplot(CD_N2O, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Daily Cow N2O",
       x = "Date", y = expression("nmol N2O m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  ) +
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = "bottom",
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

CD_N2O_plot

ggsave(filename = "plots/CD_N2O.jpeg", plot = CD_N2O_plot, width = 10, height = 8)

# GRADIENTS PLOTS
# Horse Gradient CO2 plots
HG_CO2 <- flux_data %>% filter(Campaign == "Gradient", gastype == "CO2", date %in% c("13jun", "14jun",  "16jul", "17jul", "18jul", "29jul", "30jul"), Animal == "Horse")

HG_CO2_plot <- ggplot(HG_CO2, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free_x")+
  labs(title = "Gradient Horse CO2",
       x = "Date", y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

HG_CO2_plot

ggsave(filename = "plots/HG_CO2.jpeg", plot = HG_CO2_plot, width = 10, height = 8)

# Horse Gradient CH4 plot
HG_CH4 <- flux_data %>% filter(Campaign == "Gradient", gastype == "CH4", date %in% c("13jun", "14jun",  "16jul", "17jul", "18jul", "29jul", "30jul"), Animal == "Horse")

HG_CH4_plot <- ggplot(HG_CH4, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Gradient Horse CH4",
       x = "Date", y = expression("nmol CH4 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

HG_CH4_plot

ggsave(filename = "plots/HG_CH4.jpeg", plot = HG_CH4_plot, width = 10, height = 8)

# Horse Gradient N2O plot
HG_N2O <- flux_data %>% filter(Campaign == "Gradient", gastype == "N2O", date %in% c("13jun", "14jun",  "16jul", "17jul", "18jul", "29jul", "30jul"), Animal == "Horse")

HG_N2O_plot <- ggplot(HG_N2O, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Gradient Horse N2O",
       x = "Date", y = expression("nmol N2O m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

HG_N2O_plot

ggsave(filename = "plots/HG_N2O.jpeg", plot = HG_N2O_plot, width = 10, height = 8)

# Cow Gradient CO2 plots
CG_CO2 <- flux_data %>% filter(Campaign == "Gradient", gastype == "CO2", date %in% c("13jun", "14jun",  "16jul", "17jul", "18jul", "29jul", "30jul"), Animal == "Cow")

CG_CO2_plot <- ggplot(CG_CO2, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free_x")+
  labs(title = "Gradient Cow CO2",
       x = "Date", y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

CG_CO2_plot

ggsave(filename = "plots/CG_CO2.jpeg", plot = CG_CO2_plot, width = 10, height = 8)

# Cow Gradient CH4 plots
CG_CH4 <- flux_data %>% filter(Campaign == "Gradient", gastype == "CH4", date %in% c("13jun", "14jun",  "16jul", "17jul", "18jul", "29jul", "30jul"), Animal == "Cow")

CG_CH4_plot <- ggplot(CG_CH4, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Gradient Cow CH4",
       x = "Date", y = expression("nmol CH4 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )

CG_CH4_plot

ggsave(filename = "plots/CG_CH4.jpeg", plot = CG_CH4_plot, width = 10, height = 8)

# Cow Gradient N2O plots
CG_N2O <- flux_data %>% filter(Campaign == "Gradient", gastype == "N2O", date %in% c("13jun", "14jun",  "16jul", "17jul", "18jul", "29jul", "30jul"), Animal == "Cow")

CG_N2O_plot <- ggplot(CG_N2O, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
  geom_point(size = 2) +
  geom_line()+
  facet_wrap(~base_code, scales = "free")+
  labs(title = "Gradient Cow N2O",
       x = "Date", y = expression("nmol N2O m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  labs(
    color = "Light/Dark",
    shape = "Treatment"
  )+
  scale_color_manual(values = c("NEE" = "gray", "RE" = "black"))+

  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.minor = element_line(color = "gray88", linewidth = 0.05),
    panel.grid.major = element_blank(), 
    panel.background = element_rect(fill = "white")
  )

CG_N2O_plot

ggsave(filename = "plots/CG_N2O.jpeg", plot = CG_N2O_plot, width = 10, height = 8)

photosynthesis_data <- flux_data %>% filter(gastype == "CO2", Campaign == "Daily", NEERE == "NEE")
plot(photosynthesis_data$date, photosynthesis_data$photosynthesis)

photosynthesis_plot <- ggplot(photosynthesis_data, aes(x = Animal, y = photosynthesis, color = treatment)) +
  geom_boxplot() +
  labs(y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1}),
  title = "Cow vs Horse Photosynthesis") +
  theme_minimal()

ggsave(filename = "plots/HorseCow_photosynthesis.jpeg", plot = photosynthesis_plot, width = 10, height = 8)

# Dailys COW