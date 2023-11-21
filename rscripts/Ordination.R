# Trying out ordination with Davis Data

#install.packages(c("ggvegan"))
#devtools::install_github("katemmiller/NETN_forest_summaries")
devtools::install_github("gavinsimpson/ggvegan")
library(tidyverse)
library(forestMIDN)
library(vegan)
library(ggvegan)
library(readxl)
library(ggrepel)
#library(Rcpp)
# Load data ---------------------------------------------------------------

QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("cwd", "events", "qd_ch", "seeds", "qd_sp", "saps", "soil_d", "trees59", "trees59qmd", "trees", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
list2env(data, envir = .GlobalEnv)

food <- read_xlsx("C:/Users/cseirup/Documents/Personal/grad school/Stats/Food_Preference_2021_Complete.xlsx")


# Practice food preference ordination -------------------------------------
food1 <- food %>% column_to_rownames(var = "Food Item")

food_MDS <- metaMDS(food1)

plot(food_MDS)
#ordiplot(BA_comp_NMDS, type = "text")
orditorp(food_MDS, display = "species")
orditorp(food_MDS, display = "sites")


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
trees59_ha <- trees59qmd %>% group_by(Site, Species) %>% 
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

#using ggvegan package + ggplot to get a decent plot
NMDS2_gg<- fortify(BA_comp_NMDS2) #transforms ordination results form ggplot can use
names(NMDS2_gg)
NMDS2_gg2<- NMDS2_gg %>% filter(Score == "sites") 


ord_plot <- NMDS2_gg2 %>%  ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
                     geom_point()+
                     geom_text()+
                     theme_FVM()
ord_plot

#ggrepel package for optimizing labels
options(ggrepel.max.overlaps = Inf)#Prints all labels even if overlap

ord_plot2 <- NMDS2_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
                             geom_point()+
                             coord_cartesian(clip = "off") +
                             geom_text_repel(xlim = c(-Inf, NA),
                                             ylim = c(-Inf, Inf))+ #labels can overlap plot borders
                             theme_FVM()

ord_plot2
                

# Trial ordinations from Ch1 Exploratory ----------------------------------

# Species composition ordinations -----------------------------------------
#use tree_sp_ha created earlier: all live stems >2.5cm DBH summed by species, site, and sample event. 
#Includes duplicated 1959 WP data named as WP2
#Make matrix with species on one axis and site+sampling event on the other
ord_tree_sp_ha <- Comb_tree_event %>% dplyr::select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
  group_by(Site, Species, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, SampleYear, 
                Species, num_stems_ha, BA_m2ha) 

names(ord_tree_sp_ha)
ord_trees <- ord_tree_sp_ha %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
  mutate(OrdSite = case_when(SampleEventNum == 1 & Site == 'BC' ~ 'BC59',
                             SampleEventNum == 1 & Site == 'BM' ~ 'BM59',
                             SampleEventNum == 1 & Site == 'BH' ~ 'BH59',
                             SampleEventNum == 1 & Site == 'IB' ~ 'IB59',
                             SampleEventNum == 1 & Site == 'OP' ~ 'OP59',
                             SampleEventNum == 1 & Site == 'PM' ~ 'PM59',
                             SampleEventNum == 1 & Site == 'WP2' ~ 'WP59',
                             SampleEventNum == 2 & Site == 'BC' ~ 'BC20',
                             SampleEventNum == 2 & Site == 'BM' ~ 'BM20',
                             SampleEventNum == 2 & Site == 'BH' ~ 'BH21',
                             SampleEventNum == 2 & Site == 'IB' ~ 'IB21',
                             SampleEventNum == 2 & Site == 'OP' ~ 'OP20',
                             SampleEventNum == 2 & Site == 'PM' ~ 'PM20',
                             SampleEventNum == 2 & Site == 'WP' ~ 'WP21',
                             SampleEventNum == 2 & Site == 'WP2' ~ 'WP22')) 

#Calculating relative basal area, frequency, and importance value for each site by species
ord_trees2 <- ord_trees %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                         Tot_BA = sum(BA_m2ha))
ord_trees3 <- left_join(ord_trees, ord_trees2, by = c("Site", "SampleEventNum"))

ord_trees4 <- ord_trees3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2) %>% ungroup() %>% 
  droplevels()

str(ord_trees4)
table(ord_trees4$Site)
# Relative Basal Area -----------------------------------------------------

ord_BA_wide <- ord_trees4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Rel_BA) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Rel_BA,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

BA_sp_rank <- ord_BA_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_BA")
print(arrange(BA_sp_rank, desc(Rel_BA)))

#removing least common species (don't want a lot of zeros): includes top 8 species                                    
ord_BA_wide2 <- ord_BA_wide %>% dplyr::select(-c(SODE, PRPE, BECO, AMELANCHIER, ACPE, BEAL, ACSP)) 

sort(names(ord_BA_wide2))

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA

#All species included: **Not working because removed extra WP site and now there are more species than sites**
# BA_pcaALL <- princomp(ord_BA_wide, cor=T)
# summary(BA_pcaALL)#gives proportion of variance explained by each PC; 2 gives you 41%, 3 gives you 56%
# plot(BA_pcaALL)
# biplot(BA_pcaALL)
# loadings(BA_pcaALL)
# 
# #using ggfortify package + ggplot to get a decent plot
# BA_pca_ALL_df <- fortify(BA_pcaALL) %>% dplyr::select(Comp.1, Comp.2) %>%  rownames_to_column(var = "OrdSite")
# 
# #ggrepel package for optimizing labels
# options(ggrepel.max.overlaps = Inf)#Prints all labels even if overlap
# 
# BA_ord_plotALL <- BA_pca_ALL_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
#   geom_point()+
#   coord_cartesian(clip = "off") +
#   geom_text_repel(xlim = c(-Inf, NA),
#                   ylim = c(-Inf, Inf))+ #labels can overlap plot borders
#   theme_FHM()
# 
# BA_ord_plotALL

#Top 8 species
BA_pca <- princomp(ord_BA_wide2, cor=T)
summary(BA_pca)#gives proportion of variance explained by each PC; 2 gives you 59%, 3 gives you 74%
plot(BA_pca)
biplot(BA_pca)
loadings(BA_pca)

#using ggfortify package + ggplot to get a decent plot
BA_pca_df <- fortify(BA_pca) %>% dplyr::select(Comp.1, Comp.2) %>%  rownames_to_column( var = "OrdSite")

#ggrepel package for optimizing labels

BA_ord_plot <- BA_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
  theme_FHM()

BA_ord_plot

# Relative Frequency ------------------------------------------------------
ord_F_wide <- ord_trees4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Rel_freq) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Rel_freq,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

F_sp_rank <- ord_F_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_F")
print(arrange(F_sp_rank, desc(Rel_F)))

#removing least common species (don't want a lot of zeros): includes top 8 species                                    
#ord_BA_wide2 <- ord_BA_wide %>% dplyr::select(-c(SODE, PRPE, BECO, AMELANCHIER, ACPE, BEAL, ACSP)) 

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA
# F_pca <- princomp(ord_F_wide, cor=T) # all species included
# summary(F_pca)#gives proportion of variance explained by each PC; 2 gives you 47%, 3 gives you 62%
# plot(F_pca)
# biplot(F_pca)
# loadings(F_pca)
# 
# #using ggfortify package + ggplot to get a decent plot
# F_pca_df <- fortify(F_pca) %>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")
# 
# F_ord_plot <- F_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
#   geom_point()+
#   coord_cartesian(clip = "off") +
#   geom_text_repel(xlim = c(-Inf, NA),
#                   ylim = c(-Inf, Inf))+ #labels can overlap plot borders
#   theme_FHM()
# 
# F_ord_plot

# Importance Value ------------------------------------------------------
ord_IV_wide <- ord_trees4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

IV_sp_rank <- ord_IV_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_IV")
print(arrange(IV_sp_rank, desc(Rel_IV)))

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA
# IV_pca <- princomp(ord_IV_wide, cor=T) # all species included
# summary(IV_pca)#gives proportion of variance explained by each PC; 2 gives you 45%, 3 gives you 59%
# plot(IV_pca)
# biplot(IV_pca)
# loadings(IV_pca)
# 
# #using ggfortify package + ggplot to get a decent plot
# IV_pca_df <- fortify(IV_pca) %>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")
# 
# IV_ord_plot <- IV_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
#   geom_point()+
#   coord_cartesian(clip = "off") +
#   geom_text_repel(xlim = c(-Inf, NA),
#                   ylim = c(-Inf, Inf))+ #labels can overlap plot borders
#   theme_FHM()
# 
# IV_ord_plot

# Species Composition: Pooled Betula papyrifera and Betula cordifolia (changed from original all Betula) ----------------------------------

#Comb_tree_event_BET <- Comb_tree_event %>% add_column(Betula = "")
#Comb_tree_event_BET$Betula <- str_starts(Comb_tree_event_BET$Species, "BE")

#Comb_tree_event_BET2 <- Comb_tree_event_BET %>% mutate(Species = case_when(Betula == TRUE ~ "BETULA", 
#Betula == FALSE ~ Species))

table(Comb_tree_event$Species)
Comb_tree_event_BET <- Comb_tree_event %>% mutate(Species = case_when(Species == 'BEPA' ~ "BETULA",
                                                                      Species == 'BECO' ~ "BETULA",
                                                                      TRUE ~ Species))
table(Comb_tree_event_BET$Species)

#BEPA and BECO pooled, only saplings >2.5 to 10 cm
sap_sp_ha_BET <- Comb_tree_event_BET %>% filter(size_class =='d2.5_9.9') %>% 
  dplyr::select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
  group_by(Site, Species, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, SampleYear, 
                Species, num_stems_ha, BA_m2ha) 

ord_sap_BET <- sap_sp_ha_BET %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
  mutate(OrdSite = case_when(SampleEventNum == 1 & Site == 'BC' ~ 'BC59',
                             SampleEventNum == 1 & Site == 'BM' ~ 'BM59',
                             SampleEventNum == 1 & Site == 'BH' ~ 'BH59',
                             SampleEventNum == 1 & Site == 'IB' ~ 'IB59',
                             SampleEventNum == 1 & Site == 'OP' ~ 'OP59',
                             SampleEventNum == 1 & Site == 'PM' ~ 'PM59',
                             SampleEventNum == 1 & Site == 'WP' ~ 'WP59',
                             SampleEventNum == 1 & Site == 'WP2' ~ 'WPXX',
                             SampleEventNum == 2 & Site == 'BC' ~ 'BC20',
                             SampleEventNum == 2 & Site == 'BM' ~ 'BM20',
                             SampleEventNum == 2 & Site == 'BH' ~ 'BH21',
                             SampleEventNum == 2 & Site == 'IB' ~ 'IB21',
                             SampleEventNum == 2 & Site == 'OP' ~ 'OP20',
                             SampleEventNum == 2 & Site == 'PM' ~ 'PM20',
                             SampleEventNum == 2 & Site == 'WP' ~ 'WP21',
                             SampleEventNum == 2 & Site == 'WP2' ~ 'WP22')) %>% 
  filter(OrdSite != 'WPXX') #remove 1959 WP duplicate

#Calculating relative basal area, frequency, and importance value for each site by species
ord_sap_BET2 <- ord_sap_BET %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                             Tot_BA = sum(BA_m2ha))
ord_sap_BET3 <- left_join(ord_sap_BET, ord_sap_BET2, by = c("Site", "SampleEventNum"))

ord_sap_BET4 <- ord_sap_BET3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2)

# Pooled BEPA/BECO saplings, importance value  NMDS-----------------------------------------------------

ord_IV_sap_wide <- ord_sap_BET4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

IV_sap_sp_rank <- ord_IV_sap_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Imp_val")
print(arrange(IV_sap_sp_rank, desc(Imp_val)))

# all species
mMDS_IVsap <- metaMDS(ord_IV_sap_wide, distance="bray", k=2)
ordiplot(mMDS_IVsap)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVBETsap_gg <- fortify(mMDS_IVsap) #transforms ordination results form ggplot can use

mMDS_IVBETsap_gg2 <- mMDS_IVBETsap_gg %>% filter(Score == "sites")

IVBETsap_mMDS_plot <- mMDS_IVBETsap_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVBETsap_mMDS_plot

#top 8
#removing least common species (don't want a lot of zeros): includes top 8 species                                    
ord_IV_sap_wide2 <- ord_IV_sap_wide %>% dplyr::select(-c(AMELANCHIER, ACSP, SODE, PRPE)) 

mMDS_IVsap2 <- metaMDS(ord_IV_sap_wide2, distance="bray", k=2)
ordiplot(mMDS_IVsap2)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVBETsap_gg2 <- fortify(mMDS_IVsap2) #transforms ordination results form ggplot can use

mMDS_IVBETsap_gg3 <- mMDS_IVBETsap_gg2 %>% filter(Score == "sites")

IVBETsap_mMDS_plot2 <- mMDS_IVBETsap_gg3 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVBETsap_mMDS_plot2

#Why is IB59 so different?
df1 <- trees59qmd2.5 %>% filter(Site == 'IB')




# Pooled BEPA/BECO all >2.5, importance value NMDS ------------------------
sp_ha_BET <- Comb_tree_event_BET %>% dplyr::select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
  group_by(Site, Species, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, SampleYear, 
                Species, num_stems_ha, BA_m2ha) 

ord_all_BET <- sp_ha_BET %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
  mutate(OrdSite = case_when(SampleEventNum == 1 & Site == 'BC' ~ 'BC59',
                             SampleEventNum == 1 & Site == 'BM' ~ 'BM59',
                             SampleEventNum == 1 & Site == 'BH' ~ 'BH59',
                             SampleEventNum == 1 & Site == 'IB' ~ 'IB59',
                             SampleEventNum == 1 & Site == 'OP' ~ 'OP59',
                             SampleEventNum == 1 & Site == 'PM' ~ 'PM59',
                             SampleEventNum == 1 & Site == 'WP' ~ 'WP59',
                             SampleEventNum == 1 & Site == 'WP2' ~ 'WPXX',
                             SampleEventNum == 2 & Site == 'BC' ~ 'BC20',
                             SampleEventNum == 2 & Site == 'BM' ~ 'BM20',
                             SampleEventNum == 2 & Site == 'BH' ~ 'BH21',
                             SampleEventNum == 2 & Site == 'IB' ~ 'IB21',
                             SampleEventNum == 2 & Site == 'OP' ~ 'OP20',
                             SampleEventNum == 2 & Site == 'PM' ~ 'PM20',
                             SampleEventNum == 2 & Site == 'WP' ~ 'WP21',
                             SampleEventNum == 2 & Site == 'WP2' ~ 'WP22')) %>% 
  filter(OrdSite != 'WPXX') #remove 1959 WP duplicate

#Calculating relative basal area, frequency, and importance value for each site by species
ord_all_BET2 <- ord_all_BET %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                             Tot_BA = sum(BA_m2ha))
ord_all_BET3 <- left_join(ord_all_BET, ord_all_BET2, by = c("Site", "SampleEventNum"))

ord_all_BET4 <- ord_all_BET3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2)

ord_IV_all_wide <- ord_all_BET4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

IV_all_sp_rank <- ord_IV_all_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Imp_val")
print(arrange(IV_all_sp_rank, desc(Imp_val)))

# all species
mMDS_IVall <- metaMDS(ord_IV_all_wide, distance="bray", k=2)
ordiplot(mMDS_IVall)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVBETall_gg <- fortify(mMDS_IVall) #transforms ordination results form ggplot can use

mMDS_IVBETall_gg2 <- mMDS_IVBETall_gg %>% filter(Score == "sites")

IVBETall_mMDS_plot <- mMDS_IVBETall_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVBETall_mMDS_plot

#top 8
#removing least common species (don't want a lot of zeros): includes top 8 species                                    
IV_all_sp_rank <- ord_IV_all_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Imp_val")
print(arrange(IV_all_sp_rank, desc(Imp_val)))
ord_IV_all_wide2 <- ord_IV_all_wide %>% dplyr::select(-c(AMELANCHIER, PIST, BEAL, ACSP, SODE, PRPE)) 

mMDS_IVall2 <- metaMDS(ord_IV_all_wide2, distance="bray", k=2)
ordiplot(mMDS_IVall2)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVBETall_gg2 <- fortify(mMDS_IVall2) #transforms ordination results form ggplot can use

mMDS_IVBETall_gg3 <- mMDS_IVBETall_gg2 %>% filter(Score == "sites")

IVBETall_mMDS_plot2 <- mMDS_IVBETall_gg3 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVBETall_mMDS_plot2

#Why is IB59 so different?
df1 <- trees59qmd2.5 %>% filter(Site == 'IB')

# Species Composition: >10cm DBH only -------------------------------------

tree10_sp_ha <- Comb_tree_event %>% filter(size_class != 'd2.5_9.9') %>% dplyr::select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
  group_by(Site, Species, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, SampleYear, 
                Species, num_stems_ha, BA_m2ha) 

ord_10trees <- tree10_sp_ha %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
  mutate(OrdSite = case_when(SampleEventNum == 1 & Site == 'BC' ~ 'BC59',
                             SampleEventNum == 1 & Site == 'BM' ~ 'BM59',
                             SampleEventNum == 1 & Site == 'BH' ~ 'BH59',
                             SampleEventNum == 1 & Site == 'IB' ~ 'IB59',
                             SampleEventNum == 1 & Site == 'OP' ~ 'OP59',
                             SampleEventNum == 1 & Site == 'PM' ~ 'PM59',
                             SampleEventNum == 1 & Site == 'WP' ~ 'WP59',
                             SampleEventNum == 1 & Site == 'WP2' ~ 'WPXX',
                             SampleEventNum == 2 & Site == 'BC' ~ 'BC20',
                             SampleEventNum == 2 & Site == 'BM' ~ 'BM20',
                             SampleEventNum == 2 & Site == 'BH' ~ 'BH21',
                             SampleEventNum == 2 & Site == 'IB' ~ 'IB21',
                             SampleEventNum == 2 & Site == 'OP' ~ 'OP20',
                             SampleEventNum == 2 & Site == 'PM' ~ 'PM20',
                             SampleEventNum == 2 & Site == 'WP' ~ 'WP21',
                             SampleEventNum == 2 & Site == 'WP2' ~ 'WP22')) %>% 
  filter(OrdSite != 'WPXX') #remove 1959 WP duplicate

#Calculating relative basal area, frequency, and importance value for each site by species
ord_10trees2 <- ord_10trees %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                             Tot_BA = sum(BA_m2ha))
ord_10trees3 <- left_join(ord_10trees, ord_10trees2, by = c("Site", "SampleEventNum"))

ord_10trees4 <- ord_10trees3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2)

#write.csv(ord_10trees4, './tables/Ord_data_10cm.csv', row.names = FALSE)

# Importance Value for all species------------------------------------------------------
ord_10IV_wide <- ord_10trees4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

IV10_sp_rank <- ord_10IV_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_IV")
print(arrange(IV10_sp_rank, desc(Rel_IV)))

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA
IV10_pca <- princomp(ord_10IV_wide, cor=T) # all species included
summary(IV10_pca)#gives proportion of variance explained by each PC; 2 gives you 44%, 3 gives you 61%
plot(IV10_pca)
biplot(IV10_pca)
loadings(IV10_pca)

#using ggfortify package + ggplot to get a decent plot
IV10_pca_df <- fortify(IV10_pca) %>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")

IV10_ord_plot <- IV10_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
  theme_FHM()

IV10_ord_plot

# Importance Value: >10cm  w/ NETN plots ----------------------------------
### getting errors with this now (20230826, turning off so .rmd will run)

#ord_ANP4 <- ord_ANP3.5 %>%  filter(PlotCode %in% ANPsprucePlots)

#ord_ANP5 <- ord_ANP4 %>% rename(OrdSite =  PlotCode) %>% ungroup() %>% dplyr::select(-ScientificName)

#ord_10trees5 <- ord_10trees4 %>% ungroup() %>% dplyr::select(-c(Site, SiteName, SampleEventNum))

#names(ord_10trees5)
#names(ord_ANP5)

#ordNETN <- rbind(ord_10trees5, ord_ANP5)

#ordNETN_wide <- ordNETN %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
#arrange(OrdSite) %>% 
#pivot_wider(names_from = Species,
#values_from = Imp_val,
#values_fill = 0) %>% 
#column_to_rownames(var = "OrdSite")

#IVNETN_sp_rank <- ordNETN_wide %>% map_dfc(sum) %>% 
#pivot_longer(cols = everything(), names_to = "Species",
#values_to = "Rel_IV")
#print(arrange(IVNETN_sp_rank, desc(Rel_IV)))

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA
#NETN_pca <- princomp(ordNETN_wide, cor=T) # all species included
#summary(NETN_pca)#gives proportion of variance explained by each PC; 2 gives you 19%, 3 gives you 28%
#plot(NETN_pca)
#biplot(NETN_pca)
#loadings(NETN_pca)

#using ggfortify package + ggplot to get a decent plot
#NETN_pca_df <- fortify(NETN_pca) %>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")

#NETN_pca_plot <- NETN_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
#geom_point()+
#coord_cartesian(clip = "off") +
#geom_text_repel(xlim = c(NA, NA),#labels can't overlap plot borders,
#ylim = c(NA, NA)) +
#scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
#expand = expansion(mult = 0.5)) + 
# theme_FHM()

#NETN_pca_plot



# NMDS: Basal Area top 8 species >10 cm  -----------------------------------------------------------
#Shawn wants ordination per time period first but I get too little data error.

#Basal area for 2020s visits
ord_10BA_wide <- ord_10trees4 %>% ungroup() %>% 
  dplyr::select(OrdSite, Species, Rel_BA) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Rel_BA,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")



BA_10sp_rank <- ord_10BA_wide %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_BA")
print(arrange(BA_10sp_rank, desc(Rel_BA))) # Only 11 species

#removing least common species (don't want a lot of zeros): includes top 7 species                                    
ord_10BA_wide2 <- ord_10BA_wide %>% dplyr::select(-c(AMELANCHIER, ACPE, BECO, BEAL, PIST)) 
sort(names(ord_10BA_wide2))

#metaMDS: Shawn says commonly used for vegetation 
mMDS_BA10 <- metaMDS(ord_10BA_wide2, distance="bray", k=2)# error stress nearly zero
ordiplot(mMDS_BA10)
text(mMDS_BA10)


# Relative Basal Area >2.5cm DBH, splitting by event-----------------------------------------------------
#Basal Area MNDS 1959 sites, top 8 species
ord_BA_wideOLD <- ord_trees4 %>% ungroup() %>% filter(SampleEventNum == 1) %>% dplyr::select(OrdSite, Species, Rel_BA) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Rel_BA,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

BA_sp_rankOLD <- ord_BA_wideOLD %>% map_dfc(sum) %>% 
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_BA")
print(arrange(BA_sp_rankOLD, desc(Rel_BA)))

#removing least common species (don't want a lot of zeros): includes top 8 species                                    
ord_BA_wideOLD2 <- ord_BA_wideOLD %>% dplyr::select(-c(BEAL, THOC, ACSP)) 

mMDS_BAold <- metaMDS(ord_BA_wideOLD2, distance="bray", k=2)# error stress nearly zero
ordiplot(mMDS_BAold)# doesn't work

##
#Basal Area MNDS all sites, top 8 species
##
mMDS_BA <- metaMDS(ord_BA_wide2, distance="bray", k=2)#ord_BA_wide2 made earlier for pca
ordiplot(mMDS_BA)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_BA_gg <- fortify(mMDS_BA) #transforms ordination results form ggplot can use

BA_mMDS_gg2 <- mMDS_BA_gg %>% filter(Score == "sites")

BA_mMDS_plot <- BA_mMDS_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
BA_mMDS_plot

# Diameter Distribution Ordination ----------------------------------------
#will use all individual >2.5

sz_ha <- Comb_tree_event %>% dplyr::select(-c(NumSubplots, SubplotArea, Module, DBH_QMD)) %>% 
  group_by(Site, size_class, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, SampleYear, 
                size_class, num_stems_ha, BA_m2ha) 
table(sz_ha$Site)

ord_sz <- sz_ha %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
  mutate(OrdSite = case_when(SampleEventNum == 1 & Site == 'BC' ~ 'BC59',
                             SampleEventNum == 1 & Site == 'BM' ~ 'BM59',
                             SampleEventNum == 1 & Site == 'BH' ~ 'BH59',
                             SampleEventNum == 1 & Site == 'IB' ~ 'IB59',
                             SampleEventNum == 1 & Site == 'OP' ~ 'OP59',
                             SampleEventNum == 1 & Site == 'PM' ~ 'PM59',
                             SampleEventNum == 1 & Site == 'WP2' ~ 'WM59',
                             SampleEventNum == 2 & Site == 'BC' ~ 'BC20',
                             SampleEventNum == 2 & Site == 'BM' ~ 'BM20',
                             SampleEventNum == 2 & Site == 'BH' ~ 'BH21',
                             SampleEventNum == 2 & Site == 'IB' ~ 'IB21',
                             SampleEventNum == 2 & Site == 'OP' ~ 'OP20',
                             SampleEventNum == 2 & Site == 'PM' ~ 'PM20',
                             SampleEventNum == 2 & Site == 'WP2' ~ 'WM22')) 

#Calculating relative basal area, frequency, and importance value for each site by size class
ord_sz2 <- ord_sz %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                   Tot_BA = sum(BA_m2ha))
ord_sz3 <- left_join(ord_sz, ord_sz2, by = c("Site", "SampleEventNum"))

ord_sz4 <- ord_sz3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2)


# Diameter dist. Relative Basal Area NMDS-----------------------------------------------------
#not super meaningful: commenting out

# ord_BA_sz_wide <- ord_sz_BET4 %>% ungroup() %>% dplyr::select(OrdSite, size_class, Rel_BA) %>% 
#   arrange(OrdSite) %>% 
#   pivot_wider(names_from = size_class,
#               values_from = Rel_BA,
#               values_fill = 0) %>% 
#   column_to_rownames(var = "OrdSite")
# 
# mMDS_BAsz <- metaMDS(ord_BA_sz_wide, distance="bray", k=2)
# ordiplot(mMDS_BAsz)
# 
# #Plotting
# library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
# mMDS_BAsz_gg <- fortify(mMDS_BAsz) #transforms ordination results form ggplot can use
# 
# mMDS_BAsz_gg2 <- mMDS_BAsz_gg %>% filter(Score == "sites") %>% 
#                                   mutate(SiteName = str_sub(Label, 1,2)) %>% 
#                                   mutate(SampleEventNum = str_sub(Label, -2)) %>% 
#                                   mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
#                                                                     SampleEventNum <50 ~ 2)) 
# 
# mMDS_BAsz_gg2$SampleEventNum <- as.factor(mMDS_BAsz_gg2$SampleEventNum)
# mMDS_BAsz_gg2$SiteName <- as.factor(mMDS_BAsz_gg2$SiteName)
# str(mMDS_BAsz_gg2)
#   
# 
# BAsz_mMDS_plot <- mMDS_BAsz_gg2 %>% arrange(desc(Label)) %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
#   geom_point(aes(color = SampleEventNum, size = 2), show.legend = FALSE)+
#   coord_cartesian(clip = "off") +
#   geom_path(aes(group = SiteName), color="black",
#             arrow = arrow(type = "closed",
#                           length=unit(0.075, "inches")))+
#   geom_text_repel(xlim = c(NA, NA),
#                   ylim = c(NA, NA))+ 
#   scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
#                      expand = expansion(mult = 0.5)) + #labels can overlap plot borders
#   theme_FHM()
# BAsz_mMDS_plot

# Diameter dist. Relative Density Area NMDS-----------------------------------------------------

ord_den_sz_wide <- ord_sz4 %>% ungroup() %>% dplyr::select(OrdSite, size_class, Rel_freq) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = size_class,
              values_from = Rel_freq,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

mMDS_densz <- metaMDS(ord_den_sz_wide, distance="bray", k=2)
ordiplot(mMDS_densz)
mMDS_densz$stress

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_densz_gg <- fortify(mMDS_densz) #transforms ordination results form ggplot can use

mMDS_densz_gg2 <- mMDS_densz_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_densz_gg2$SampleEventNum <- as.factor(mMDS_densz_gg2$SampleEventNum)
mMDS_densz_gg2$SiteName <- as.factor(mMDS_densz_gg2$SiteName)
str(mMDS_densz_gg2)

densz_mMDS_plot <- mMDS_densz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
  geom_point(size = 4)+
  coord_cartesian(clip = "off") +
  
  geom_path(aes(group = SiteName), linewidth = .7, color="black", lineend = "butt",
            arrow = arrow(type = "closed",
                          length=unit(0.1, "inches"), 
                          ends = "first")) + 
  geom_label_repel(xlim = c(NA, NA),
                   ylim = c(NA, NA), show.legend = FALSE)+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) +#labels can overlap plot borders
  theme(legend.position = c(0.03,.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
densz_mMDS_plot

ggsave("./figures/Ord_den_sizeclass.jpg", densz_mMDS_plot, height = 5, width = 9, dpi = 400)

# Diameter dist. importance value  NMDS-----------------------------------------------------

ord_IV_sz_wide <- ord_sz4 %>% ungroup() %>% dplyr::select(OrdSite, size_class, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = size_class,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

mMDS_IVsz <- metaMDS(ord_IV_sz_wide, distance="bray", k=2)
ordiplot(mMDS_IVsz)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVsz_gg <- fortify(mMDS_IVsz) #transforms ordination results form ggplot can use

mMDS_IVsz_gg2 <- mMDS_IVsz_gg %>% filter(Score == "sites")

IVsz_mMDS_plot <- mMDS_IVsz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVsz_mMDS_plot

# Diameter dist. importance value  NMDS, attempting only 2020's visits-----------------------------------------------------

ord_IV_sz_wide2 <- ord_sz4 %>% filter(SampleEventNum == 2) %>%  ungroup() %>% dplyr::select(OrdSite, size_class, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = size_class,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

mMDS_IVsz2 <- metaMDS(ord_IV_sz_wide2, distance="bray", k=2)
ordiplot(mMDS_IVsz2)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVsz_gg2 <- fortify(mMDS_IVsz2) #transforms ordination results form ggplot can use

mMDS_IVsz_gg4 <- mMDS_IVsz_gg2 %>% filter(Score == "sites")

IVsz_mMDS_plot2 <- mMDS_IVsz_gg4 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVsz_mMDS_plot2

# Diameter dist. importance value  NMDS, attempting only 1959's visits-----------------------------------------------------

ord_IV_sz_wide1 <- ord_sz4 %>% filter(SampleEventNum == 1) %>%  ungroup() %>% dplyr::select(OrdSite, size_class, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = size_class,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

mMDS_IVsz1 <- metaMDS(ord_IV_sz_wide1, distance="bray", k=2)
ordiplot(mMDS_IVsz1)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_IVsz_gg1 <- fortify(mMDS_IVsz1) #transforms ordination results form ggplot can use

mMDS_IVsz_gg3 <- mMDS_IVsz_gg1 %>% filter(Score == "sites")

IVsz_mMDS_plot1 <- mMDS_IVsz_gg3 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
IVsz_mMDS_plot1

# Alternate Ordination Methods --------------------------------------------


#PCoA based on Bray-Curtis distances
#BA_PCoA<-cmdscale(vegdist(ord_trees_wide4),k=2)
#library(MASS)
#eqscplot(BA_PCoA,type="n")
#text(BA_PCoA,rownames(ord_trees_wide4))

#metaMDS
##BA_mMDS<-metaMDS(ord_trees_wide4,distance="bray",k=2)
#ordiplot(BA_mMDS)
#text(BA_mMDS)
#First plotting attempt with mMDS, but sounds like I should be using PCA
#BA_mMDS_gg <- fortify(BA_mMDS) #transforms ordination results form ggplot can use

#BA_mMDS_gg2 <- BA_mMDS_gg %>% filter(Score == "sites")


#ord_plot <- BA_mMDS_gg2 %>%  ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
#geom_point()+
# geom_text()+
#theme_FHM()
#ord_plot

#ggrepel package for optimizing labels
#options(ggrepel.max.overlaps = Inf)#Prints all labels even if overlap

#ord_plot2 <- NMDS_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
# geom_point()+
# coord_cartesian(clip = "off") +
#  geom_text_repel(xlim = c(-Inf, NA),
#  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
# theme_FHM()

#ord_plot2

