# Trying out ordination with Davis Data

#install.packages(c("vegan"))
#devtools::install_github("katemmiller/forestMIDN")

library(tidyverse)
library(forestMIDN)
library(vegan)
# Load data ---------------------------------------------------------------

QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("qd_ch", "seeds", "qd_sp", "saps", "soil_d", "trees59", "trees", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
list2env(data, envir = .GlobalEnv)

# Setting up matrix for ordination  ---------------------------------------
# Will do species comp for just 2020/2021 first: Proportion of basal area for each species at each site
# Problem that the sample areas are different 2020/21 vs 1959

#combining 2020 trees and sapling (>2.4cm) to match 1959 tree definition
saps2 <- saps %>% drop_na() %>% mutate(stem=1) %>% 
                                select(Site, Transect, Species, DBH, stem) # dropped quadrants with no species
treesL <- filter(trees, Status == "L") %>% droplevels() %>% 
                                           mutate(stem=1) %>% 
                                           select(Site, Transect, Species, DBH, stem)
names(saps2);names(treesL)
Ltreesap <- rbind(treesL, saps2)
Ltreesap2 <- Ltreesap %>% filter(DBH > 2.4) # Davis included trees/saps over 1in DBH
Ltreesap2$BA_cm2 <- round(pi * ((Ltreesap2$DBH/2)^2), 4)# add column for basal area

#2020 plot area is 2000 m2; conversion factor to hectare is 5; converting from kg to megagram
Ltreesap_ha <- Ltreesap2 %>% group_by(Site, Species) %>% 
                             summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000)

#dropping frequency, only looking at basal area
Ltreesap_ba_wide <- Ltreesap_ha %>% select(Site, Species, BA_m2ha) %>% 
                                    arrange(Site) %>% 
                                    pivot_wider(names_from = Species,
                                                values_from = BA_m2ha,
                                                values_fill = 0)

Ltreesap_ba_wide2 <- Ltreesap_ba_wide %>% mutate(Tot_BA = sum(c_across())) %>% 
                                          transmute(
                                            ABBA_prop = round((ABBA/Tot_BA), digits = 7),
                                            ACER_prop = round((ACER/Tot_BA), digits = 7),
                                            ACPE_prop = round((ACPE/Tot_BA), digits = 7),
                                            ACRU_prop = round((ACRU/Tot_BA), digits = 7),
                                            AMELANCHIER_prop = round((AMELANCHIER/Tot_BA), digits = 7),
                                            BEAL_prop = round((BEAL/Tot_BA), digits = 7),
                                            BECO_prop = round((BECO/Tot_BA), digits = 7),
                                            BEPA_prop = round((BEPA/Tot_BA), digits = 7),
                                            PIGL_prop = round((PIGL/Tot_BA), digits = 7),
                                            PIRU_prop = round((PIRU/Tot_BA), digits = 7),
                                            PIST_prop = round((PIST/Tot_BA), digits = 7),
                                            PRPE_prop = round((PRPE/Tot_BA), digits = 7),
                                            SODE_prop = round((SODE/Tot_BA), digits = 7),
                                            THOC_prop = round((THOC/Tot_BA), digits = 7),
                                            TSCA_prop = round((TSCA/Tot_BA), digits = 7)) %>% 
                                            column_to_rownames(var = "Site")
head(Ltreesap_ba_wide2)

#removed species that were <1% at any site
Ltreesap_ba_wide3 <- Ltreesap_ba_wide2 %>% map_dfc(sum)
Ltreesap_ba_wide4 <- Ltreesap_ba_wide2 %>% select(-c(ACPE_prop, AMELANCHIER_prop, PRPE_prop, SODE_prop, BECO_prop, ACER_prop))

BA_comp_NMDS1 <- metaMDS(Ltreesap_ba_wide2)
BA_comp_NMDS2 <- metaMDS(Ltreesap_ba_wide4)#removed rare species
plot(BA_comp_NMDS1)
plot(BA_comp_NMDS2)
#ordiplot(BA_comp_NMDS, type = "text")
orditorp(BA_comp_NMDS2, display = "sites")
#orditorp(BA_comp_NMDS, display = "species")
