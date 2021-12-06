# Trying out ordination with Davis Data

install.packages(c("Rcpp"))
#devtools::install_github("katemmiller/NETN_forest_summaries")
#devtools::install_github("gavinsimpson/ggvegan")
library(tidyverse)
library(forestMIDN)
library(vegan)
library(ggvegan)
library(readxl)
library(ggrepel)
library(Rcpp)
# Load data ---------------------------------------------------------------

QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("cwd", "qd_ch", "seeds", "qd_sp", "saps", "soil_d", "trees59", "trees59_calc", "trees", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
list2env(data, envir = .GlobalEnv)

food <- read_xlsx("C:/Users/cseirup/Documents/Personal/grad school/Stats/Food_Preference_2021_Complete.xlsx")


# Practice food preference ordination -------------------------------------
food1 <- food %>% column_to_rownames(var = "Food Item")

food_MDS <- metaMDS(food1)

plot(food_MDS)
#ordiplot(BA_comp_NMDS, type = "text")
orditorp(food_MDS, display = "species")
#orditorp(food_MDS, display = "sites")


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
                             summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000) %>%
                             mutate(Site = recode(Site, BC = "BC20", BH = "BH20", BM = "BM20", 
                                                        OP = "OP20", PM = "PM20", WP = "WP20",
                                                        IB = "IB20"))

#setting up 1959 tree data
trees59_ha <- trees59_calc %>% group_by(Site, Species) %>% 
                               summarise(num_stems_ha = sum(num_stem) * 20, BA_m2ha = sum(BA_cm2)/500) %>%
                               mutate(Site = recode(Site, BC = "BC59", BH = "BH59", BM = "BM59", 
                                                          OP = "OP59", PM = "PM59", WP = "WP59",
                                                          IB = "IB59"))

#combine 1959 and 2020/21 data for ordination
ord_trees <- rbind(Ltreesap_ha, trees59_ha)

#dropping frequency, only looking at basal area
ord_trees_wide <- ord_trees %>% select(Site, Species, BA_m2ha) %>% 
                                    arrange(Site) %>% 
                                    pivot_wider(names_from = Species,
                                                values_from = BA_m2ha,
                                                values_fill = 0)

ord_trees_wide2 <- ord_trees_wide %>% mutate(Tot_BA = sum(c_across())) %>% 
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
names(ord_trees_wide2)
ord_trees_wide3 <- ord_trees_wide2 %>% map_dfc(sum) %>% 
                                           pivot_longer(cols = c(ABBA_prop:TSCA_prop), names_to = "Species",
                                                                 values_to = "tot_prop")
print(arrange(ord_trees_wide3, desc(tot_prop)))

#removing least common species (don't want a lot of zeros)                                    
ord_trees_wide4 <- ord_trees_wide2 %>% select(-c(SODE_prop, PRPE_prop, ACER_prop, BECO_prop, 
                                                     AMELANCHIER_prop, ACPE_prop, BEAL_prop)) 

names(ord_trees_wide4)

BA_comp_NMDS1 <- metaMDS(ord_trees_wide2)
BA_comp_NMDS2 <- metaMDS(ord_trees_wide4)#removed rare species
plot(BA_comp_NMDS1)
plot(BA_comp_NMDS2) 
#ordiplot(BA_comp_NMDS2, type = "text")
orditorp(BA_comp_NMDS2, display = "sites")
#orditorp(BA_comp_NMDS, display = "species")

#using ggvegan plotting
NMDS2_gg<- fortify(BA_comp_NMDS2)
names(NMDS2_gg)
NMDS2_gg2<- NMDS2_gg %>% filter(Score == "sites") 


ord_plot <- NMDS2_gg2 %>%  ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
                     geom_point()+
                     geom_text()+
                     theme_FVM()
ord_plot
                   
                
