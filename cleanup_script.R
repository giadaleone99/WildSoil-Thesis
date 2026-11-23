### Script with all data processing steps, resulting in clean data files 

library(dplyr)
#library(patchwork)
#library(car)

# soil data -------------------------------------------------------------
soil_data_raw <- read.csv("data/soil_data_raw.csv")
date_data <- readRDS("flux_data/clean_flux_data.rds")

# get the dates per base code for dung age
date_data <- date_data %>% 
  filter(gastype == "CO2") %>% 
  filter(NEERE == "RE") %>% 
  filter(treatment == "F") %>% 
  filter(Days_Since_First == 0) %>% 
  dplyr::select(base_code, longdate)

date_data <- date_data %>% 
  mutate(days_since_first = as.numeric(as.Date("2024-08-02") - as.Date(longdate))) %>% 
  mutate(dung_age = case_when(days_since_first %in% c(3, 4) ~ "3-4", days_since_first == 18 ~ "18", days_since_first %in% c(49, 50) ~ "49-50")) %>% 
  dplyr::select(base_code, dung_age)

# create combined dataset with dung plot age
comb_soil_data <- soil_data_raw %>% 
  mutate(Animal = case_when(grepl("^C", base_code) ~ "Cow", grepl("^H", base_code) ~ "Horse"),) %>% 
  mutate(Animal = as.factor(Animal))

comb_soil_data$sample_type <- factor(comb_soil_data$sample_type, levels = c("Dung soil", "Fresh", "Control"))

comb_soil_data <- comb_soil_data %>% 
  left_join(date_data, by = "base_code")

write.csv(comb_soil_data, file="clean_data/soil_data.csv")

# vegetation data -------------------------------------------------------

# Import data
vegdung_lab <- read.csv("data/Plant_Dung_CN_pH_elements.csv", sep = ";", check.names = FALSE)
veg_raw <- read.csv("data/vegetation_data.csv")
fieldwork_data_raw <- read.csv("data/Fieldwork_data_final.csv")
species_list <- read.csv("data/species_lists.csv")
date_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")

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
  dplyr::select(Kode_1, 'N%', CN_ratio, B, Na, Mg, Al, P, S, K, Ca, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, As, Se, Sr, Mo, Cd, Ba, Tl, Pb) %>% 
  filter(!Kode_1 %in% c("CG dung", "HG dung", "HOD dung", "COD dung", "HND dung", "CND dung"))
colnames(veglab_data)[1] <- c("plot_id")

veg_combined <- veg_combined %>% 
  left_join(veglab_data, by = "plot_id") %>% 
  dplyr::rename("TN" = "N%")

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
  dplyr::select(base_code, longdate)

date_data <- date_data %>% 
  mutate(days_since_first = as.numeric(as.Date("2024-08-02") - as.Date(longdate))) %>% 
  mutate(dung_age = case_when(days_since_first %in% c(3, 4) ~ "3-4", days_since_first == 18 ~ "18", days_since_first %in% c(49, 50) ~ "49-50")) %>% 
  dplyr::select(base_code, dung_age)

veg_combined <- veg_combined %>% 
  left_join(date_data, by = "base_code")

# Filtering for veg height growth and bar plot
veg_growth <- veg_summary %>%
  filter(!is.na(veg_height)) %>%
  mutate(height_value = veg_height_2 - veg_height,
         height_type = "veg_growth") 

veg_growth <- veg_growth %>% 
  left_join(date_data, by = "base_code")

# data for stacked species plot
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
  ) %>% 
  mutate(
    base_code = substr(plot_id, 1, 3)
  ) %>% 
  left_join(date_data, by = "base_code")

veg_percent <- veg_percent %>%
  mutate(dung_age = as.factor(dung_age))

# vegetation RDA data

# get the veg species data and transform it into a matrix
veg_now <- veg_raw %>%
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
    weight = dry_weight/harvested_area_m2) %>% 
  mutate(
    base_code = substr(plot_id, 1, 3)
  )

veg_now <- veg_now %>% 
  left_join(date_data, by = "base_code")

species_df <- veg_now[c("plot_id", "veg_class", "weight")]

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
         Animal = as.factor(case_when(grepl("^C", plotID) ~ "Cattle", grepl("^H", plotID) ~ "Horse")),
         sample_type = as.factor(sample_type))

soil_data <- soil_data %>% 
  left_join(date_data, by = "base_code")

soil_data <- soil_data %>% 
  mutate(dung_age = as.factor(dung_age))

soil_data <- soil_data[c("plotID", "bulk_density", "pH", "PO4.P", "TN", "Campaign", "Animal", "sample_type", "dung_age")]

soil_df <- soil_data[order(soil_data$plotID), ]
rownames(soil_df) <- NULL

soil_df <- soil_df[c("bulk_density", "pH", "PO4.P", "TN", "Campaign", "Animal", "sample_type", "dung_age")]

write.csv(veg_percent, file="clean_data/veg_species.csv")
write.csv(veg_combined, file="clean_data/veg_data.csv")
write.csv(veg_growth, file="clean_data/veg_growth.csv")
write.csv(presence_absence_matrix, file="clean_data/RDA_presence.csv")
write.csv(species_df, file="clean_data/RDA_species.csv")
write.csv(soil_df, file="clean_data/RDA_soil.csv")

# GHG data --------------------------------------------------------------
allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
par_data <- read.csv("data/PAR_file.csv", sep = ",")
start_data <- read.csv("data/Fieldwork_data_final.csv") %>% 
  dplyr::select(UniqueID, Start_timehhmmss)

# add the GHG data to dfs
CO2_PS_data <- allflux_data %>% 
  filter(gastype == "CO2" & NEERE == "NEE") %>% 
  dplyr::rename(CO2_PS_flux = best.flux)
CO2_RE_data <- allflux_data %>% 
  filter(gastype == "CO2" & NEERE == "RE") %>%
  dplyr::rename(CO2_RE_flux = best.flux)
CH4_data <- allflux_data %>% 
  filter(gastype == "CH4") %>% 
  dplyr::rename(CH4_flux = best.flux)
N2O_data <- allflux_data %>% 
  filter(gastype == "N2O") %>% 
  dplyr::rename(N2O_flux = best.flux)

# PAR stuff
CO2_PS_data <- CO2_PS_data %>% 
  left_join(start_data, by = "UniqueID") 
CO2_PS_data <- CO2_PS_data %>% 
  mutate(minute = str_sub(Start_timehhmmss, 1, 5),
         par = NA)
for (i in 1:nrow(CO2_PS_data)) {
  if (substring(CO2_PS_data$minute[i], 1, 1) == "0") {
    CO2_PS_data$minute[i] = substring(CO2_PS_data$minute[i], 2)
  }
}

# Join gas data per date
flux_data <- CO2_PS_data %>% 
  left_join(CO2_RE_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-21:-25, -27:-35)
flux_data <- flux_data %>% 
  left_join(CH4_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-22:-26, -28:-36)
flux_data <- flux_data %>% 
  left_join(N2O_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-1:-5, -7:-9, -12:-17, -23, -29)
flux_data <- flux_data[, c(9, 10, 2, 11, 12, 14, 15, 16, 3, 20, 21, 1, 7, 8, 13, 17, 18, 19, 4, 5, 6)]

# more PAR stuff
par_data <- par_data %>% 
  mutate(date = str_sub(min.time, 1, 10),
         min = str_sub(min.time, 12, 16))
par_data <- par_data %>% 
  mutate(longdate = as.Date(par_data$date, format = "%d-%m-%Y"))
for (j in 1:nrow(flux_data)) {
  for (i in 1:nrow(par_data)) {
    if (par_data$longdate[i] == flux_data$longdate[j] && par_data$min[i] == flux_data$minute[j]) {
      flux_data$par[j] <- mean(c(par_data$PAR_umol_m2_s[i], par_data$PAR_umol_m2_s[i+1], par_data$PAR_umol_m2_s[i+2], par_data$PAR_umol_m2_s[i+3], par_data$PAR_umol_m2_s[i+4]))
    }
  }
} 

# add bulk density to df
soil_data <- soil_data_raw %>% 
  filter(sample_type %in% c("Fresh", "Control"))
bulkdensity <- soil_data %>% 
  dplyr::select(plotID, bulk_density)
flux_data <- flux_data %>% 
  left_join(bulkdensity, by = "plotID")

# create data transformed columns for CH4 and N2O fluxes
flux_data <- flux_data %>% 
  mutate(ranked_CH4_flux = rank(CH4_flux),
         normalized_CH4_flux = bestNormalize(CH4_flux)$x.t,
         addition_CH4_flux = CH4_flux + 2.904) %>% 
  mutate(ranked_N2O_flux = rank(N2O_flux),
         normalized_N2O_flux = bestNormalize(N2O_flux)$x.t,
         addition_N2O_flux = N2O_flux + 1.471)

write.csv(flux_data, file="clean_data/flux_data.csv")

# combined data ---------------------------------------------------------

# Import data
allflux_data <- readRDS("flux_data/clean_flux_data.rds")
soil_data_raw <- read.csv("data/soil_data_raw.csv")
par_data <- read.csv("data/PAR_file.csv", sep = ",")
start_data <- read.csv("data/Fieldwork_data_final.csv") %>% 
  dplyr::select(UniqueID, Start_timehhmmss)
veg_growth_data <- readRDS("plant_data/veg_growth_data.rds") %>% 
  dplyr::rename("plotID"= "plot_id") %>% 
  dplyr::select(-2:-9)
dung_lab <- read.csv("data/Plant_Dung_CN_pH_elements.csv", sep = ";", check.names = FALSE) %>% 
  filter(Kode_3 %in% "Dung") %>% 
  rename("base_code" = "Kode_2",
         "sample_type" = "Kode_3",
         "plotID" = "Kode_1") %>% 
  dplyr::select(-4,-5,-9)
veg_combined_data <- readRDS("plant_data/veg_combined_data.rds")%>% 
  dplyr::rename("plotID" = "plot_id",
                "plant_CN" = "CN_ratio") %>% 
  dplyr::select(-5)

# get the dates per base code for days since first
date_data <- allflux_data %>% 
  filter(gastype == "CO2") %>% 
  filter(NEERE == "RE") %>% 
  filter(treatment == "F") %>% 
  filter(Days_Since_First == 0) %>% 
  dplyr::select(base_code, longdate)

date_data <- date_data %>% 
  mutate(days_since_first = as.numeric(as.Date("2024-08-02") - as.Date(longdate))) %>% 
  mutate(dung_age = case_when(days_since_first %in% c(3, 4) ~ "3-4", days_since_first == 18 ~ "18", days_since_first %in% c(49, 50) ~ "49-50")) %>% 
  dplyr::select(base_code, dung_age)

# Combine data
CO2_PS_data <- allflux_data %>% 
  filter(gastype == "CO2" & NEERE == "NEE") %>% 
  dplyr::rename(CO2_PS_flux = best.flux)
CO2_RE_data <- allflux_data %>% 
  filter(gastype == "CO2" & NEERE == "RE") %>%
  dplyr::rename(CO2_RE_flux = best.flux)
CH4_data <- allflux_data %>% 
  filter(gastype == "CH4") %>% 
  dplyr::rename(CH4_flux = best.flux)
N2O_data <- allflux_data %>% 
  filter(gastype == "N2O") %>% 
  dplyr::rename(N2O_flux = best.flux)

# PAR stuff
CO2_PS_data <- CO2_PS_data %>% 
  left_join(start_data, by = "UniqueID") 
CO2_PS_data <- CO2_PS_data %>% 
  mutate(minute = str_sub(Start_timehhmmss, 1, 5),
         par = NA)

for (i in 1:nrow(CO2_PS_data)) {
  if (substring(CO2_PS_data$minute[i], 1, 1) == "0") {
    CO2_PS_data$minute[i] = substring(CO2_PS_data$minute[i], 2)
  }
}

# Join gas data per date
flux_data <- CO2_PS_data %>% 
  left_join(CO2_RE_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-21:-25, -27:-35)
flux_data <- flux_data %>% 
  left_join(CH4_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-22:-26, -28:-36)
flux_data <- flux_data %>% 
  left_join(N2O_data, by = c("plotID", "longdate")) %>% 
  dplyr::select(-1:-5, -7:-9, -12:-17, -23, -29)
flux_data <- flux_data[, c(9, 10, 2, 11, 12, 14, 15, 16, 3, 20, 21, 1, 7, 8, 13, 17, 18, 19, 4, 5, 6)]

# more PAR stuff
par_data <- par_data %>% 
  mutate(date = str_sub(min.time, 1, 10),
         min = str_sub(min.time, 12, 16))
par_data <- par_data %>% 
  mutate(longdate = as.Date(par_data$date, format = "%d-%m-%Y"))

for (j in 1:nrow(flux_data)) {
  for (i in 1:nrow(par_data)) {
    if (par_data$longdate[i] == flux_data$longdate[j] && par_data$min[i] == flux_data$minute[j]) {
      flux_data$par[j] <- mean(c(par_data$PAR_umol_m2_s[i], par_data$PAR_umol_m2_s[i+1], par_data$PAR_umol_m2_s[i+2], par_data$PAR_umol_m2_s[i+3], par_data$PAR_umol_m2_s[i+4]))
    }
  }
} 

# create a df summing the gas data on the diff days for each plot
plot_flux_data <- flux_data %>% 
  group_by(plotID) %>% 
  dplyr::summarise(CO2_PS_sum = sum(CO2_PS_flux), 
                   CO2_RE_sum = sum(CO2_RE_flux), 
                   CH4_sum = sum(CH4_flux), 
                   N2O_sum = sum(N2O_flux), .groups = 'drop')

# join with veg growth
plot_flux_data <- plot_flux_data %>% 
  left_join(veg_growth_data, by = "plotID")

soil_data <- soil_data_raw %>% 
  filter(sample_type %in% c("Fresh", "Control"))

plot_flux_data <- plot_flux_data %>% 
  left_join(veg_combined_data, by = "plotID")
plot_flux_data <- plot_flux_data %>% 
  full_join(soil_data, by = "plotID")

soil_data_test <- soil_data %>% 
  left_join(plot_flux_data, by = "plotID", suffix = c("_soil", "_veg")) %>% 
  mutate(Animal = case_when(
    grepl("^C", plotID) ~ "Cow",
    grepl("^H", plotID) ~ "Horse")) %>% 
  dplyr::rename("base_code" = "base_code_soil")

soil_data_test <- soil_data_test %>% 
  left_join(date_data, by = "base_code")

# add bulk density to gas data
bulkdensity <- soil_data %>% 
  select(plotID, bulk_density)

flux_data <- flux_data %>% 
  left_join(bulkdensity, by = "plotID")

# create data transformed columns for CH4 and N2O fluxes
flux_data <- flux_data %>% 
  mutate(ranked_CH4_flux = rank(CH4_flux),
         normalized_CH4_flux = bestNormalize(CH4_flux)$x.t,
         addition_CH4_flux = CH4_flux + 2.904) %>% 
  mutate(ranked_N2O_flux = rank(N2O_flux),
         normalized_N2O_flux = bestNormalize(N2O_flux)$x.t,
         addition_N2O_flux = N2O_flux + 1.471)

# create df with only the last GHG measurements per plot and soil variables, does not have dung soil rn
last_flux_data <- flux_data %>% 
  filter(period == 3) %>% 
  filter(Days_Since_First %in% c(3, 4, 16, 46, 47)) %>% 
  filter(UniqueID != 'CD3_C_RE_01aug' & UniqueID != 'CD4_C_RE_01aug' & UniqueID != 'HD3_C_RE_01aug' & UniqueID != 'HD4_C_RE_01aug' & UniqueID != 'CD3_F_RE_01aug' & UniqueID != 'CD4_F_RE_01aug' & UniqueID != 'HD3_F_RE_01aug' & UniqueID != 'HD4_F_RE_01aug')

last_combined <- last_flux_data %>% 
  left_join(soil_data, by = 'plotID')

write.csv(last_combined, file="clean_data/combined_flux_data.csv")
write.csv(soil_data_test, file="clean_data/combined_soil_data.csv")

