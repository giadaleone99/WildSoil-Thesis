# script to process ghg flux data

library(stringr)
library(dplyr)
library(ggplot2)
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
HD <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Horse")

ggplot(HD, aes(x = longdate, y = best.flux, color = NEERE, group = plotNEERE, shape = treatment)) +
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
  theme(
    legend.position = c(0.75, 0.25),
    scale_shape_manual(values = c(16, 17)),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )
  
HD_CH4 <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Horse")

ggplot(HD, aes(x = longdate, y = best.flux, color = treatment, group = plotID)) +
  geom_point() +
  geom_line()+
  facet_wrap(~plotID, scales = "free_x")+
  labs(title = "Daily Horse CO2",
       x = "Date", y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )





# Cow daily plots
CD <- flux_data %>% filter(Campaign == "Daily", gastype == "CO2", date %in% c("15jul",  "16jul", "17jul", "18jul", "19jul", "29jul", "30jul", "31jul", "01aug", "02aug"), Animal == "Horse", NEERE == "NEE")

ggplot(HD, aes(x = longdate, y = best.flux, color = treatment, group = plotID)) +
  geom_point() +
  geom_line()+
  facet_wrap(~plotID, scales = "free_x")+
  labs(title = "Daily Horse CO2",
       x = "Date", y = expression(mu * "mol CO2 m"^{-2} * " s"^{-1})) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black"),   
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),               # No major grid lines
    panel.background = element_rect(fill = "white")   # White background
  )



