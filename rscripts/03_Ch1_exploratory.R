
#install.packages("forcats")
#devtools::install_github("katemmiller/forestNETN")
library(ggalt)
library(forestNETN)
library(knitr)
library(kableExtra)
library(ggfortify)
library(ggrepel)
library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)

source("C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj/rscripts/davis_functions.R")

# Load data ---------------------------------------------------------------
QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("cwd", "BrLi59", "events", "qd_ch", "seeds", "qd_sp", "rwlong", "saps", "seeds59", "soil_d",
               "cores", "trees59", "trees59qmd", "trees","coreIDtrees", "ibutDAYall", "ibutMONall","ibutDAY2020", "ibutMON2020", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
#list2env(data, envir = .GlobalEnv)

# Common plot labels ------------------------------------------------------
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 
                'PM' = "Pemetic Mtn", 'BH' = 'Bass Harbor Head', 
                'IB' = 'Ironbound Island', 'WP' = '1St W. Mtn Plot',
                'WP2' = 'Western Mtn')

sp_names <- c('PIRU' = 'Picea rubens', 'PIGL' = 'Picea glauca', 'TSCA' = 'Tsuga canadensis', 
              'ABBA' = "Abies balsamea", 'ACRU' = 'Acer rubrum', 'BEPA' = 'Betula papyrifera',
              'PIST' = 'Pinus strobus', 'BEAL' = 'Betula alleghaniensis', 'THOC' = 'Thuja occidentalis',
              'BECO' = 'Betula cordifolia', 'ACPE' = 'Acer pensylvanicum', 'AMELANCHIER' = 'Amelanchier spp', 
              'ACSP' = 'Acer spicatum', 'SODE' = 'Sorbus decora', "PRPE" = 'Prunus pensylvanica')

# Stem map ----------------------------------------------------------------
list2env(data["trees"], envir = .GlobalEnv)
trees <- trees %>% filter(Site != "WP") #Remove first Western Mtn plot that was replaced with WP2
trees$Site <- ordered(trees$Site,
                      levels = c("BH", "OP", "BM", "PM", "BC",   "IB", "WP2"))

stem_mapV <- trees %>% 
  ggplot(aes(x = cY, y = cX))+
  geom_hline(yintercept = 5, linetype = 2)+
  geom_vline(xintercept = 100, linetype = 2)+ 
  geom_point(aes(color = Status, size = DBH))+
  #coord_equal()+# makes x and y the same scale
  scale_color_manual(name = "Status", labels = c("Dead", "Live"), 
                     values = c("D" = '#808080', "L" = '#31a354'))+
  scale_size_continuous(range = c(.5, 3.5))+
  facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
  labs(x = "Easting (meters)", y = "Northing (meters)")+ 
  theme(axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), #facet wrap text size
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "#696969", size = 0.01),
        axis.ticks = element_line(color = "#696969", size = 0.5),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5)),
        aspect.ratio = 10/3,
        legend.position = "right")+
  
  coord_flip()+
  theme_FHM()
stem_mapV

# Sapling Heat Maps -------------------------------------------------------
list2env(data["saps"], envir = .GlobalEnv)
saps <- saps %>% filter(Site != "WP") #Remove first Western Mtn plot that was replaced with WP2
saps$Site <- ordered(saps$Site,
                     levels = c("BH", "OP", "BM", "PM", "BC",   "IB", "WP2"))

sapsb <- saps %>% mutate(X = case_when(Transect == 1|Transect == 3 ~ 2.5, 
                                       Transect == 2|Transect == 4 ~ 7.5)) %>%  
  mutate(Y = case_when(Transect == 3|Transect == 4 ~ Quadrat+100,
                       Transect == 1|Transect == 2 ~ Quadrat+0)) %>%
  drop_na() %>% 
  mutate(num_stems = 1) #drop_na deletes quadrats with no saplings recorded

#summarize by site and remove WP (oddly shaped plot, would have to be handled separately)
sapsc <- sapsb %>% filter(DBH >= 2.5) %>% group_by(Site, Transect, Quadrat)  %>% 
                                       summarise(sum_stems = sum(num_stems),
                                                X = first(X),
                                                Y = first(Y))

sap_map <- ggplot(sapsc, aes(X,Y)) +
  geom_raster(aes(fill = sum_stems), hjust = 1, vjust = 0.5)+
  facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
  scale_fill_gradient(low = "white", high = "#006600")+
  theme(axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), #facet wrap text size
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "#696969", size = 0.01),
        axis.ticks = element_line(color = "#696969", size = 0.5),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5)),
        #aspect.ratio = 10/3,
        legend.position = "right")+
  theme_FHM()

sap_map

# Combined tree stem map + sap heat map -----------------------------------
sapsd <- sapsc %>% rename(cX = X, cY = Y) %>% dplyr::select(-Quadrat) #renameing columns to match trees
treesB <- trees %>% dplyr::select(Site, Transect, Status, DBH, cX, cY) 

treesap_map <- ggplot(treesB, aes(x = cY, y = cX))+
  geom_raster(data = sapsd, aes(fill = sum_stems), hjust = 1, vjust = 0.5)+
  geom_hline(yintercept = 5, linetype = 2)+
  geom_vline(xintercept = 100, linetype = 2)+ 
  scale_fill_gradient(name = "Saplings (no. stems)", low = "white", high = "#3182bd")+
  geom_point(aes(color = Status, size = DBH))+
  scale_color_manual(name = "Tree Status", labels = c("Dead", "Live"), 
                     values = c("D" = '#808080', "L" = '#31a354'))+
  scale_size_continuous(name = 'Tree DBH (cm)', range = c(.5, 3.5))+
  facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
  labs(x = "Northing (meters)", y = "Easting (meters)")+ 
  theme(axis.text = element_text(size = 10), 
        strip.text = element_text(size = 9), #facet wrap text size
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "#696969", size = 0.01),
        axis.ticks = element_line(color = "#696969", size = 0.5),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5)),
        aspect.ratio = 10/3,
        legend.position = "right")+
  coord_flip()+
  theme_FHM()

treesap_map
#ggsave("./figures/TreeSapMaps.jpg", treesap_map, height = 5, width = 10, dpi = 400)

# Combine 2020/21 trees and saps to match Davis size classes---------------------------
#using Davis size classes, only live trees
saps2 <- saps %>% drop_na() %>% filter(DBH>=2.5)#only taking >2.5cm to match Davis size classes
treesL_df1 <- filter(trees, Status == "L") %>% droplevels()  #only live trees

treesL_df1.5 <- treesL_df1 %>% dplyr::select(Site, Transect, Species, DBH, SampleYear, SampleEventNum)
saps3 <- saps2 %>% dplyr::select(Site, Transect, Species, DBH, SampleYear, SampleEventNum) 
Ltreesap <- rbind(treesL_df1.5, saps3)
Ltreesap$BA_cm2 <- round(pi * ((Ltreesap$DBH/2)^2), 4)

Ltreesap2 <- Ltreesap %>% mutate(size_class = as.factor(case_when(between(DBH, 2.5, 9.9) ~ "d2.5_9.9", 
                                                                  between(DBH, 10, 19.9) ~ "d10_19.9", 
                                                                  between(DBH, 20, 29.9) ~ "d20_29.9", 
                                                                  between(DBH, 30, 39.9) ~ "d30_39.9", 
                                                                  between(DBH, 40, 49.9) ~ "d40_49.9", 
                                                                  between(DBH, 50, 59.9) ~ "d50_59.9", 
                                                                  between(DBH, 60, 69.9) ~ "d60_69.9", 
                                                                  DBH >= 70 ~ "d70", TRUE ~ "unknown")), num_stem = 1)
# Calculating Quadratic Mean Diameter -------------------------------------------------
#Calculating the QMD for each size class using the 2020 diameters, should be better than the mid-point of each size class. 
#Earlier version was by Species and size class
qmd1 <- Ltreesap2 %>% group_by(size_class) %>% summarise(QMD = (sum(DBH^2)),
                                                         avg_dbh = mean(DBH, na.rm = TRUE),
                                                         se_dens = sd(DBH, na.rm = TRUE)/
                                                           sqrt(sum(!is.na(DBH))),
                                                         num_indiv = sum(!is.na(DBH)))
qmd2 <- qmd1
qmd2$QMD <-(qmd2$QMD = sqrt(qmd2$QMD/qmd2$num_indiv)) # had to split QMD equation into 2 parts

#Adding calculated size class QMD to 1959 tree data
list2env(data["trees59"], envir = .GlobalEnv)
trees59b <- trees59 
#combine two smallest size classes to be more consistent with the size of the other size classes
trees59b$size_class <- trees59b$size_class %>% recode(d2.5_4.9 = "d2.5_9.9") %>% 
  recode(d5_9.9 = "d2.5_9.9")

trees59qmd <- left_join(trees59b, qmd2, by = c("size_class"))

midpoints <- data.frame("size_class" = c("d2.5_9.9", "d10_19.9", 
                                         "d20_29.9", "d30_39.9", "d40_49.9", 
                                         "d50_59.9", "d60_69.9", "d70"), 
                        midpoint = c(6.2, 15, 25, 35, 45, 55, 65, 70))

trees59qmd <- left_join(trees59qmd, midpoints, by = "size_class")

trees59qmd2 <- trees59qmd %>% dplyr::select("Site", "Species", "size_class", "num_stem",  
                                     "QMD", "midpoint", "SampleYear", "SampleEventNum")

trees59qmd2$QMD <- coalesce(trees59qmd2$QMD, trees59qmd2$midpoint) #replacing missing QMDs with midpoints for that size class/species
trees59qmd2$BA_cm2 <- round(pi * (((trees59qmd2$QMD/2)^2)*trees59qmd2$num_stem), 4) #basal area per site

# Combine 1959 and 2020/21 >2.4 tree data ---------------------------------
intersect(names(trees59qmd2), names(Ltreesap2))

Ltreesap2.5 <- Ltreesap2 %>% dplyr::select(Site, Species, SampleYear, SampleEventNum, 
                                    size_class, num_stem, BA_cm2, DBH) %>% 
  rename(DBH_QMD = DBH)

trees59qmd2.5 <- trees59qmd2 %>% dplyr::select(Site, Species, SampleYear, SampleEventNum, 
                                        size_class, num_stem, BA_cm2, QMD) %>% 
  rename(DBH_QMD = QMD)

Comb_tree <- rbind(Ltreesap2.5, trees59qmd2.5)

# Combine with event data
list2env(data["events"], envir = .GlobalEnv)
event_tree <- events %>% filter(Module == "Trees")

Comb_tree_event <- left_join(Comb_tree, event_tree, by = c("Site","SampleYear","SampleEventNum"))

sort(unique(Comb_tree_event$size_class))
Comb_tree_event$size_class <- ordered(Comb_tree_event$size_class,
                                      levels = c("d2.5_9.9", "d10_19.9", "d20_29.9", 
                                                 "d30_39.9", "d40_49.9", "d50_59.9", 
                                                 "d60_69.9", "d70_79.9", "d80_89.9", 
                                                 "d90_99.9", "d100p"))
levels(Comb_tree_event$size_class) #Size classes are in the order we want them now

Comb_tree_event$Site <- ordered(Comb_tree_event$Site,
                                levels = c("BH", "OP", "BM", "PM", "BC", "IB", "WP2"))
levels(Comb_tree_event$Site)#Sites are in the ordered by complexity now

Comb_tree_event$SampleEventNum <- as.character(Comb_tree_event$SampleEventNum)
# Calculate basal area and density by site and size class for diameter distributions---------------------------
tree_dist_ha <- Comb_tree_event %>% dplyr::select(-c(NumSubplots, SubplotArea, Module, Species, DBH_QMD)) %>% 
  group_by(Site, size_class, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = round(sum_stems * (10000/TotArea), digits = 0), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  mutate(log10stems_ha = log10(num_stems_ha)) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, SampleYear, 
         size_class, num_stems_ha, BA_m2ha,log10stems_ha) 
  

# Plot diameter distribution comparison -----------------------------------
#Connected lines
Comp_tree_dist_plot <- ggplot(data = tree_dist_ha, aes(color = SampleEventNum, x = size_class, y = num_stems_ha))+
  geom_point()+ 
  geom_line(aes(group = SampleEventNum), linewidth = 2)+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 9), # change axis label size
        strip.text = element_text(size = 10), # change facet text size
        axis.title = element_text(size = 12), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(1,0),
        legend.justification = c(1,0))+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_discrete(labels= c('2.5-10','10-20',
                             '20-30','30-40','40-50',
                             '50-60','60-70','70+'))+ 
  theme_FHM() 

Comp_tree_dist_plot 
#ggsave("./figures/DiameterDist_lines.jpg", Comp_tree_dist_plot, height = 5, width = 10, dpi = 400)

#no saps
dist_plot_no_saps <- tree_dist_ha %>% filter(size_class != "d2.5_9.9") %>% 
  ggplot(aes(color = SampleEventNum, x = size_class, y = num_stems_ha))+
  geom_point()+ 
  geom_line(aes(group = SampleEventNum), linewidth = 2)+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 9), # change axis label size
        strip.text = element_text(size = 10), # change facet text size
        axis.title = element_text(size = 12), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(1,0),
        legend.justification = c(1,0))+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_discrete(labels= c('10-20',
                             '20-30','30-40','40-50',
                             '50-60','60-70','70+'))+ 
  theme_FHM() 

dist_plot_no_saps 


#log10 diameter distribution
Log10_tree_dist_plot <- ggplot(data = tree_dist_ha, aes(color = SampleEventNum, x = size_class, y = log10stems_ha))+
  geom_point()+ 
  geom_line(aes(group = SampleEventNum), linewidth = 2)+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "log10(stems/ha)")+ 
  theme(axis.text = element_text(size = 9), # change axis label size
        strip.text = element_text(size = 10), # change facet text size
        axis.title = element_text(size = 12), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  #legend.position = c(1,0),
  #legend.justification = c(1,0))+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_discrete(labels= c('2.5-10','10-20',
                             '20-30','30-40','40-50',
                             '50-60','60-70','70+'))+ 
  theme_FHM() 

Log10_tree_dist_plot  



#staggered bar plots w/ no saps
#Fill in zeros for missing species detentions (for better plotting)
tree_dist_ha2 <- tree_dist_ha %>% dplyr::select(-c(BA_m2ha, log10stems_ha)) %>% 
                                  pivot_wider(names_from = size_class, values_from = num_stems_ha, values_fill = 0)

tree_dist_ha3 <- tree_dist_ha2 %>% pivot_longer(cols = 6:10, names_to = "size_class", values_to = "num_stems_ha")

#bar graph version
Comp_tree_dist_bar <- tree_dist_ha3 %>% filter(size_class != "d2.5_9.9") %>%  
  ggplot(aes(x = size_class, y = num_stems_ha, fill = SampleEventNum))+
  geom_bar(width = .75, position = position_dodge(), stat = 'identity')+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 9), # change axis label size
        strip.text = element_text(size = 10), # change facet text size
        axis.title = element_text(size = 12), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(0.97,0),
        legend.justification = c(1,0))+
  scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020/21'), values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_discrete(labels= c('10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
  theme_FHM() 

Comp_tree_dist_bar

#ggsave("./figures/DiameterDist_bars.jpg", Comp_tree_dist_bar, height = 5, width = 10, dpi = 400)
# ibuttons ----------------------------------------------------------------
#library(MASS)
list2env(data["ibutDAYall"], envir = .GlobalEnv)
list2env(data["ibutMONall"], envir = .GlobalEnv)

str(ibutDAYall)

p.ibut <- ggplot(data = ibutDAYall,  aes(x = day, y = airTmin))+
                   facet_wrap(~site_NSEW, ncol = 4)+
                   theme_FHM()
                 
p.ibut

#ibut <- stl(ibutDAYall, s.window='periodic')

# Table 2: Calculate density, ba, biomass, carbon mass ---------------------------
#use Comb_tree_event df of trees and saps >2.4 cm created earlier
#DEAD TREES NOT INCLUDED FROM 2020S SAMPLE. NOT COLLECTED IN 1959
SumTable <- Comb_tree_event %>% filter(size_class != 'd2.5_9.9') %>% 
  mutate(carbon = 0) %>% 
  mutate(comp = 0)

SumTable$comp <- mapply(compmass,Species=SumTable$Species, DBH = SumTable$DBH_QMD)
SumTable$carbon <- mapply(compcarb,Species=SumTable$Species, biomass = SumTable$comp)

names(SumTable)
#summing metrics by site and calculating total carbon and biomass for 1959 plots by * by number of stems per site
SumTable_site <- SumTable %>% dplyr::select(-c(NumSubplots, SubplotArea, Module, Species, size_class, DBH_QMD)) %>% 
  group_by(Site, SampleEventNum) %>% 
  summarise(SampleYear = first(SampleYear),
            SiteName = first(SiteName),
            sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            sum_carb = sum(carbon*num_stem),
            sum_comp = sum(comp*num_stem),
            TotArea = first(TotArea)) 

#converting from plot to ha and kg to megagram 
SumTable_ha <- SumTable_site %>%  mutate(density_ha = sum_stems * (10000/TotArea), 
                                         BA_m2ha = sum_BA_cm2/TotArea,
                                         carbonmass_Mgha = (sum_carb * (10000/TotArea))/1000, 
                                         biomass_Mgha = (sum_comp * (10000/TotArea))/1000) %>% 
  mutate(across(where(is.numeric), round, 0)) %>% 
  dplyr::select(Site, SiteName, SampleEventNum, 
         density_ha, BA_m2ha, biomass_Mgha, carbonmass_Mgha) 

SumTable_ha$SampleEventNum <- as.character(SumTable_ha$SampleEventNum)

SumTable_ha2 <- SumTable_ha %>% pivot_wider(names_from = SampleEventNum, 
                                            values_from = c(density_ha, BA_m2ha, 
                                                            biomass_Mgha, carbonmass_Mgha)) 

#Adding sapling stems/ha for both visits to Table 2
sapTable <- Comb_tree_event %>% filter(size_class == 'd2.5_9.9') %>% 
  dplyr::select(-c(NumSubplots, SubplotArea, Module, Species, size_class, DBH_QMD)) %>% 
  group_by(Site, SampleEventNum) %>% 
  summarise(SampleYear = first(SampleYear),
            SiteName = first(SiteName),
            sum_stems = sum(num_stem),
            TotArea = first(TotArea)) 

#converting from plot to ha and kg to megagram 
sapTable_ha <- sapTable %>%  mutate(density_ha = sum_stems * (10000/TotArea)) %>%  
                              mutate(across(where(is.numeric), round, 0)) %>% 
                              dplyr::select(Site, SiteName, SampleEventNum, 
                                                density_ha) 

sapTable_ha$SampleEventNum <- as.character(sapTable_ha$SampleEventNum)

sapTable_ha2 <- sapTable_ha %>% pivot_wider(names_from = SampleEventNum, 
                                            values_from = c(density_ha)) %>% 
                                rename("sap_den_1"= "1") %>% 
                                rename("sap_den_2"= "2")



SumTable_ha3 <- left_join(SumTable_ha2, sapTable_ha2, by = c("Site", "SiteName"))

#Adding 1959 mosses and lichen
list2env(data["BrLi59"], envir = .GlobalEnv)

BrLi59sum <- BrLi59 %>% group_by(Site) %>% summarise(BriLicCov = sum(Cover))

SumTable_ha4 <- left_join(SumTable_ha3, BrLi59sum, by = "Site")

#adding bryophyte/lichen
list2env(data["qd_ch"], envir = .GlobalEnv)
qd_ch <- qd_ch %>% filter(Site != "WP")
CoverClasses <- data.frame("Cover" = 0:9, 
                           cClasses = c("0", "<1", "1-2", "2-5", "5-10", "10-25", "25-50",
                                        "50-75", "75-95", "95-100"),
                           cMidpoints = c("0", "0.5", "1.5", "3.5", "7.5", "17.5", "37.5",
                                          "62.5", "87.5", "97.5"))

qd_ch2 <- left_join(qd_ch, CoverClasses, by = "Cover")
qd_ch2$cMidpoints <- as.numeric(qd_ch2$cMidpoints)
str(qd_ch2)
#add lichen, ns spagnum bryophite, and sphaghum together per quad
BryLic <- qd_ch2 %>% filter(Character == "Lichen" | Character == "Sphagnum"| Character == "NS Bryophyte") %>% 
  group_by(Site, SampleEventNum, Quadrat_Num) %>% 
  summarise(BryLicCov = sum(cMidpoints),
            n = n())

#Summarize by Site
BryLic2 <- BryLic %>% group_by(Site, SampleEventNum) %>% 
  summarise(avgBryLicCov = sum(BryLicCov)/n(),
            num_quads = n(),
            se_cov = sd(BryLicCov)/sqrt(n()),
            sd_cover = sd(BryLicCov))

SumTable_ha5 <- left_join(SumTable_ha4, dplyr::select(BryLic2, avgBryLicCov), by = "Site")

#Adding standing dead tree BA to table
names(events)
names(trees)
dtrees <- trees %>% filter(Status == "D") %>% left_join(event_tree, by = c("Site", "SampleEventNum", "SampleYear"))
dtrees$BA_cm2 <- round(pi * ((dtrees$DBH/2)^2), 4)
names(dtrees)
dtreeTable <- dtrees %>%  dplyr::select(Site,SiteName,BA_cm2,SampleEventNum,SampleYear,TotArea) %>% 
  group_by(Site, SampleEventNum) %>% 
  summarise(SampleYear = first(SampleYear),
            SiteName = first(SiteName),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea)) 

#converting from plot to ha 
dtreeTable_ha <- dtreeTable %>%  mutate(BA_m2ha = sum_BA_cm2/TotArea) %>%  
  mutate(across(where(is.numeric), round, 1)) %>% 
  dplyr::select(Site, SiteName,  BA_m2ha) 

dtreeTable_ha2 <- dtreeTable_ha %>% rename("2020s SDW BA_m2ha"= BA_m2ha)


SumTable_ha6 <- left_join(SumTable_ha5, dtreeTable_ha2, by = c("Site", "SiteName"))

#Adding CWD
list2env(data["cwd"], envir = .GlobalEnv)
cwd <- cwd %>% filter(Site != "WP")

SumTable_ha7 <- left_join(SumTable_ha6, cwd, by = "Site")

SumTable_ha8 <- SumTable_ha7 %>% ungroup() %>% 
  mutate(vol_m3ha = round(vol_m3ha, digits = 1)) %>% 
  mutate(avgBryLicCov = round(avgBryLicCov, digits = 1)) %>% 
  dplyr::select(-c("carbonmass_Mgha_2", "carbonmass_Mgha_1", "num_pieces", "SiteName"))



kable(SumTable_ha8)
#write.csv(Comb_sum_table, './tables/Summary_metrics1959_2020_2.5.csv', row.names = FALSE)


# Table 1 -----------------------------------------------------------------

dSites <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Site_description_table.xlsx")
kable(dSites)

# Table 3: core summary -----------------------------------------------------------------

sCores <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Core_site_summary_table.xlsx")
kable(sCores)


# Table 4: Site Abbreviations ---------------------------------------------
sCodes <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Site_codes_table.xlsx")
kable(sCores)



# Species comparison: Basal area trees and saps, density seeds ----------------------------------------
#Species composition by site: 
#trees
tree_sp_ha <- Comb_tree_event %>% filter(size_class != "d2.5_9.9") %>% 
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

#Fill in zeros for missing species detentions (for better plotting)
tree_sp_ha2 <- tree_sp_ha %>% dplyr::select(-num_stems_ha) %>% pivot_wider(names_from = Species, values_from = BA_m2ha, values_fill = 0)
tree_sp_haC <- tree_sp_ha2 %>% pivot_longer(cols = ABBA:BECO, names_to = "Species", values_to = "BA_m2ha")

#saps
sap_sp_ha <- Comb_tree_event %>% filter(size_class == 'd2.5_9.9') %>% 
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
#Fill in zeros for missing species detentions (for better plotting)
sap_sp_ha2 <- sap_sp_ha %>% dplyr::select(-BA_m2ha) %>% pivot_wider(names_from = Species, values_from = num_stems_ha, values_fill = 0)
sap_sp_haC <- sap_sp_ha2 %>% pivot_longer(cols = ABBA:ACSP, names_to = "Species", values_to = "num_stems_ha")

#seeds
list2env(data["seeds"], envir = .GlobalEnv)
list2env(data["seeds59"], envir = .GlobalEnv)

seeds2 <- seeds %>% filter(Latin_name != "No species") %>%
                    group_by(Site, Latin_name) %>% 
                    summarise(Sum_stem = sum(Count),
                              SampleYear = first(SampleYear),
                              SampleEventNum = first(SampleEventNum)) %>% 
                    mutate(den_m2 = Sum_stem/30) # 30 m2 quadrats, convert to per m2

seeds3 <- seeds2 %>% dplyr::select(Site, Latin_name, den_m2, SampleYear, SampleEventNum) %>%
                      ungroup() 

seeds4 <- rbind(seeds3, seeds59)
seeds4$SampleEventNum <- as.character(seeds4$SampleEventNum)
#Fill in zeros for missing species detentions (for better plotting)

seeds4w <- seeds4 %>% pivot_wider(names_from = Latin_name, values_from = den_m2, values_fill = 0)
seedsC <- seeds4w %>% pivot_longer(cols = 4:13, names_to = "Species", values_to = "den_m2")

table(seedsC$Species)

#Trees
BET_spp <- c("BEAL", "BECO", "BEPA")
LOWCAN <- c("ACPE", "ACSP", "AMELANCHIER", "PRPE", "SODE")
OTHCON <-c("PIST","THOC")
HRDWDS <- c("ACRU", "ACPE", "ACSP", "AMELANCHIER", "PRPE", "SODE","BEAL", "BECO", "BEPA")

tree_sp_haB <- tree_sp_ha %>% filter(Site != "WP") %>% 
                  mutate(sppgrp = case_when(Species %in% HRDWDS ~ paste0("Deciduous spp"),
                                            Species %in% OTHCON ~ paste0("Other conifer spp"),
                                            TRUE ~ paste0(Species)))

tree_sp_haB$sppgrp <- ordered(tree_sp_haB$sppgrp,
                             levels = c("PIRU", "PIGL", "ABBA", "TSCA",
                                        "Other conifer spp", "Deciduous spp"))

pTreeSp <- tree_sp_haB %>% filter(Site != "WP") %>% 
          ggplot(aes(x = SampleEventNum, y = BA_m2ha, fill = sppgrp))+
              geom_col(position = position_stack(reverse = TRUE))+
              facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
  scale_fill_manual(name = "Species", labels = c("PIRU" = 'Picea rubens', "PIGL" = 'Picea glauca',
                                                  'ABBA' = 'Abies balsamea', 'TSCA' = 'Tsuga canadensis'), 
                     values = c("PIRU" = '#35b15a', "PIGL" = '#86ce7e',
                                'ABBA' = '#ceecc5', 'TSCA' = '#d8e6f3', 
                                'Deciduous spp' = "#d9c4ed", 'Other conifer spp' = '#9daee2'))+
  scale_x_discrete(labels= c('1' = '1959', '2' = '2020s'))+ 
  xlab('Year')+
  ylab('Basal area ('~m^2*'/ha)')+
            theme_FHM()

 pTreeSp  
  
#Saps
table(sap_sp_ha$Species)
sap_sp_haB <- sap_sp_ha %>% filter(Site != "WP") %>% 
  mutate(sppgrp = case_when(Species %in% HRDWDS ~ paste0("Deciduous spp"),
                            Species %in% OTHCON ~ paste0("Other conifer spp"),
                            TRUE ~ paste0(Species)))
 
sap_sp_haB$sppgrp <- ordered(sap_sp_haB$sppgrp,
                             levels = c("PIRU", "PIGL", "ABBA", "TSCA",
                                       "Other conifer spp", "Deciduous spp"))

 
 pSapSp <- sap_sp_haB %>% filter(Site != "WP") %>% 
   ggplot(aes(x = SampleEventNum, y = num_stems_ha, fill = sppgrp))+
   geom_col(position = position_stack(reverse = TRUE))+
   facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
   scale_fill_manual(name = "Species", labels = c("PIRU" = 'Picea rubens', "PIGL" = 'Picea glauca',
                                                  'ABBA' = 'Abies balsamea', 'TSCA' = 'Tsuga canadensis'), 
                     values = c("PIRU" = '#35b15a', "PIGL" = '#86ce7e',
                                'ABBA' = '#ceecc5', 'TSCA' = '#d8e6f3', 
                                'Deciduous spp' = "#d9c4ed", 'Other conifer spp' = '#9daee2'))+
   scale_x_discrete(labels= c('1' = '1959', '2' = '2020s'))+ 
   xlab('Year')+
   ylab('Density (stems/ha)')+
   theme_FHM()
 
 pSapSp 
   
 #Seeds
table(seeds4$Latin_name)
sBET_spp <- c("Betula papyrifera")
sLOWCAN <- c("Acer pensylvanicum","Amelanchier","Sorbus decora")
sOTHCON <-c("Pinus strobus")
sHRDWDS <- c("Acer pensylvanicum","Amelanchier","Sorbus decora", "Betula papyrifera", "Acer rubrum")


seed_sp_haB <- seeds4 %>% filter(Site != "WP") %>% 
  mutate(sppgrp = case_when(Latin_name %in% sHRDWDS ~ paste0("Deciduous spp"),
                            Latin_name %in% sOTHCON ~ paste0("Other conifer spp"),
                            TRUE ~ paste0(Latin_name)))

seed_sp_haB$sppgrp <- ordered(seed_sp_haB$sppgrp,
                             levels = c("Picea rubens", "Picea glauca", "Abies balsamea", "Tsuga canadensis",
                                        "Other conifer spp", "Deciduous spp"))

seed_sp_haB$Site <- ordered(seed_sp_haB$Site,
                                levels = c("BH", "OP", "BM", "PM", "BC", "IB", "WP2"))

pSeedSp <- seed_sp_haB %>% filter(Site != "WP") %>% 
  ggplot(aes(x = SampleEventNum, y = den_m2, fill = sppgrp))+
  geom_col(position = position_stack(reverse = TRUE))+
  facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
  scale_fill_manual(name = "Species", values = c("Picea rubens" = '#35b15a', "Picea glauca" = '#86ce7e',
                               'Abies balsamea' = '#ceecc5', 'Tsuga canadensis' = '#d8e6f3', 
                               'Deciduous spp' = "#d9c4ed", 'Other conifer spp' = '#9daee2'))+
  scale_x_discrete(labels= c('1' = '1959', '2' = '2020s'))+ 
  xlab('Year')+
  ylab('Density (stems/'~m^2*')')+
  theme_FHM()

pSeedSp 

#Combine into one plot using patchwork package

pComp <- pTreeSp/pSapSp/pSeedSp +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')
pComp

#Plotting percent relative change in composition: Use __C dataset with added zeros: tree_sp_haC, sap_sp_haC, seedsC
#Calculate total BA/density per site and divide to get % BA by species per visit
#Trees
tree_sp_haC2 <- tree_sp_haC %>% group_by(Site, SampleEventNum) %>% summarise(Tot_BA = sum(BA_m2ha)) 

tree_sp_haC3 <- left_join(tree_sp_haC, tree_sp_haC2, by = c("Site", "SampleEventNum"))

tree_sp_haC4 <- tree_sp_haC3 %>% mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
                                 mutate(stage = "Tree") %>% 
  filter(Site != "WP")
                                 
tree_sp_haC5 <- tree_sp_haC4 %>% dplyr::select(-c(BA_m2ha, Tot_BA, SampleYear)) %>% 
                                 pivot_wider(names_from = SampleEventNum, values_from = Rel_BA) %>% 
                                 dplyr::select("Site", "Species", "stage", "1","2") %>% droplevels()
                              
#Saps                 
sap_sp_haC2 <- sap_sp_haC %>% group_by(Site, SampleEventNum) %>% summarise(Tot_den = sum(num_stems_ha))

sap_sp_haC3 <- left_join(sap_sp_haC, sap_sp_haC2, by = c("Site", "SampleEventNum"))

sap_sp_haC4 <- sap_sp_haC3 %>% mutate(Rel_den = num_stems_ha/Tot_den) %>% 
  mutate(stage = "Sapling")%>% 
  filter(Site != "WP")

sap_sp_haC5 <- sap_sp_haC4 %>% dplyr::select(-c(num_stems_ha, Tot_den, SampleYear)) %>% 
                                pivot_wider(names_from = SampleEventNum, values_from = Rel_den)%>% 
                                dplyr::select("Site", "Species", "stage", "1","2") %>% droplevels()

#Seeds
seedsC2 <- seedsC %>% group_by(Site, SampleEventNum) %>% summarise(Tot_den = sum(den_m2))

seedsC3 <- left_join(seedsC, seedsC2, by = c("Site", "SampleEventNum"))

seedsC4 <- seedsC3 %>% mutate(Rel_den = den_m2/Tot_den) %>% 
  mutate(stage = "Seedling")%>% 
  filter(Site != "WP")

seedsC5 <- seedsC4 %>% dplyr::select(-c(den_m2, Tot_den, SampleYear)) %>% 
  pivot_wider(names_from = SampleEventNum, values_from = Rel_den) %>% 
  dplyr::select("Site", "Species", "stage", "1","2") %>% droplevels()

seedsC5$Species <- recode(seedsC5$Species, 'Picea rubens' = "PIRU")
seedsC5$Species <- recode(seedsC5$Species, "Picea glauca" = "PIGL")
seedsC5$Species <- recode(seedsC5$Species, "Abies balsamea" = "ABBA")
seedsC5$Species <- recode(seedsC5$Species, "Tsuga canadensis" = "TSCA")
seedsC5$Species <- recode(seedsC5$Species, "Amelanchier" = "AMELANCHIER")
seedsC5$Species <- recode(seedsC5$Species, "Acer rubrum" = "ACRU")
seedsC5$Species <- recode(seedsC5$Species, "Betula papyrifera" = "BEPA")
seedsC5$Species <- recode(seedsC5$Species, "Acer pensylvanicum" = "ACPE")
seedsC5$Species <- recode(seedsC5$Species, "Pinus strobus" = "PIST")
seedsC5$Species <- recode(seedsC5$Species, "Sorbus decora" = "SODE")

#Then combine trees/saps/seeds and subtract 1 - 2 to get % change
TSScomb <- rbind(tree_sp_haC5, sap_sp_haC5, seedsC5)
TSScomb2 <- TSScomb %>% rename(Sample1 = "1") %>% rename(Sample2 = "2") %>% 
                        mutate(perRelDiff = Sample2 - Sample1)# %>% 
                        #filter(perRelDiff != 0)

#Fill in zeros for missing species detentions (for better plotting)
TSScomb2w <- TSScomb2 %>% dplyr::select(-c(Sample1, Sample2)) %>%  
                          pivot_wider(names_from = stage, values_from = perRelDiff, values_fill = 0) %>% 
                          mutate(sum = Tree+Sapling+Seedling) %>% 
                          filter(sum != 0)

TSScomb3 <- TSScomb2w %>% dplyr::select(-c(sum)) %>% pivot_longer(cols = 3:5, names_to = "stage", values_to = "perRelDiff")
TSScomb3$stage <- ordered(TSScomb3$stage,levels = c("Tree", "Sapling", "Seedling"))
TSScomb3all <- TSScomb3
TSScomb3gr <- TSScomb3
#plotting
#all species
TSScomb3all$Species <- ordered(TSScomb3all$Species,
                              levels = c("PIRU", "PIGL", "ABBA", "TSCA", "ACRU", "BEAL","BECO", "BEPA",
                                         "PIST","THOC", "ACPE", "ACSP", "AMELANCHIER", "PRPE", "SODE"))

pTSS <- TSScomb3all %>% 
  ggplot(aes(x = Species, y = perRelDiff*100, group = stage, fill = stage, alpha = stage))+
  geom_col(position = position_dodge2())+
  facet_wrap(~Site, ncol = 2, labeller = as_labeller(site_names), scales = "free_x")+
  geom_hline(yintercept = 0)+
  scale_fill_manual(name = "Stage", 
                    values = c("Seedling" = "#e5f5e0", "Sapling" = '#a1d99b', "Tree" = '#31a354'))+
  scale_alpha_manual(values = c(.8,.6,.8), guide = "none")+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        axis.text.x = element_text(face = "italic"))+
  scale_x_discrete(labels = sp_names, guide = guide_axis(n.dodge=2))+
  ylab(bquote('% Relative Change'))+ 
  theme_FHM()

pTSS

#less common species grouped
BET_spp <- c("BEAL", "BECO", "BEPA")
LOWCAN <- c("ACPE", "ACSP", "AMELANCHIER", "PRPE", "SODE")
OTHCON <-c("PIST","THOC")
HRDWDS <- c("ACRU", "ACPE", "ACSP", "AMELANCHIER", "PRPE", "SODE","BEAL", "BECO", "BEPA")

TSScomb3gr2 <- TSScomb3gr %>% mutate(sppgrp = case_when(Species %in% HRDWDS ~ paste0("Deciduous spp"),
                                                        Species %in% OTHCON ~ paste0("Other conifer spp"),
                                                        TRUE ~ paste0(Species)))

TSScomb3gr2$sppgrp <- ordered(TSScomb3gr2$sppgrp,
                              levels = c("PIRU", "PIGL", "ABBA", "TSCA",
                                         "Other conifer spp", "Deciduous spp"))


names(TSScomb3gr2)
TSScomb3gr3 <- TSScomb3gr2 %>% group_by(Site, stage, sppgrp) %>% 
                               summarize(tot_diff = sum(perRelDiff))


pTSSgr <- TSScomb3gr3 %>% filter(Site == "BH"|Site == "OP") %>% 
  ggplot(aes(x = sppgrp, y = tot_diff*100, group = stage, fill = stage, alpha = stage))+
  geom_col(position = position_dodge2())+
  facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names), scales = "free_x")+
  geom_hline(yintercept = 0)+
  scale_fill_manual(name = "Stage", 
                    values = c("Seedling" = "#F0E68C", "Sapling" = '#3182bd', "Tree" = '#31a354'))+
  scale_alpha_manual(values = c(.8,.6,.8), guide = "none")+
  theme(axis.text.x = element_text(face = "italic"))+
  xlab('Species')+
  ylab(bquote('% Relative Change'))+
  scale_x_discrete(labels = sp_names, guide = guide_axis(n.dodge=2))+
  theme_FHM()

pTSSgr

#original colors and horizontal for presentation
TSScomb3gr3$Site <- ordered(TSScomb3gr3$Site, levels = c("BH", "OP", "BM", "PM", "BC", "IB", "WP2"))
table(TSScomb3gr3$Site)
str(TSScomb3gr3)
pTSSgr2 <- TSScomb3gr3 %>%  
  ggplot(aes(x = sppgrp, y = tot_diff*100, group = stage, fill = stage, alpha = stage))+
  geom_col(position = position_dodge2())+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names), scales = "free_x")+
  geom_hline(yintercept = 0)+
  scale_fill_manual(name = "Stage", 
                    values = c("Seedling" = "#e6ecff", "Sapling" = '#809fff', "Tree" = '#3366ff'))+
  scale_alpha_manual(values = c(.8,.6,.8), guide = "none")+
  theme(axis.text.x = element_text(face = "italic"))+
  xlab('Species')+
  ylab(bquote('% Relative Change'))+
  scale_x_discrete(labels = sp_names, guide = guide_axis(n.dodge=2))+
  theme(legend.position = c(0.97,0),
         legend.justification = c(1,0))+
  theme_FHM()

pTSSgr2
ggsave("./figures/SpeciesCompChange.jpg", pTSSgr2, height = 5, width = 11, dpi = 400)

# Table 5 Overstory composition -------------------------------------------
tree_sp_haC6 <- tree_sp_haC4
tree_sp_haC6$Site <- ordered(tree_sp_haC6$Site, levels = c("BH", "OP", "BM", "PM", "BC", "IB", "WP2"))
tree_sp_haC6$Rel_BA <- tree_sp_haC6$Rel_BA*100
tComp <- tree_sp_haC6 %>% ungroup() %>%  unite("SiteYear", 1,3, remove = FALSE) %>% 
                          arrange(Site, SiteYear) %>% 
                          dplyr::select(SiteYear, Species, Rel_BA) %>% 
                          pivot_wider(names_from = SiteYear,
                                      values_from = Rel_BA) 


sp_names2 <- as.data.frame(sp_names)
sp_names2 <- rownames_to_column(sp_names2, "Species")
tComp2 <- left_join(tComp, sp_names2, by = "Species")
str(tComp2)

tComp3 <- tComp2 %>% select(16, 2:15) %>% 
  group_by(sp_names) %>% 
  mutate(Mean = mean(c_across(2:14)))  %>% mutate(across(where(is.numeric), ~round(., 1)))

names(tComp3)
head(tComp3)

tComp3$BC_1 <- paste0(tComp3$BC_1, ".0")

tComp4 <- tComp3 %>% arrange(desc(Mean))
 kable(tComp4)

# Bootstrapping trees -----------------------------------------------------

#remove WP data (first plot, 20 x 100 m) and live trees only to match Davis
names(trees)
trees2 <- trees %>% filter(Site != "WP") %>%  filter(Status == "L") %>% 
  dplyr::select(-c(X,Y)) #may have to drop level, WP still floating around
table(trees2$Site)

#set up dataset for later  calculations
trees3 <- trees2 %>% mutate(carbon = 0) %>% mutate(comp = 0) %>% mutate(num_stem = 1)

trees3$comp <- mapply(compmass,Species=trees3$Species, DBH = trees3$DBH)
trees3$carbon <- mapply(compcarb,Species=trees3$Species, biomass = trees3$comp)
trees3$BA_cm2 <- round(pi * ((trees3$DBH/2)^2), 4)

# Split up tree data into 10x10m subplots --------------------------------
trees4 <- trees3 %>% mutate(SubPlotID = case_when(cY < 10 ~ 1, 
                                                  cY >= 10 & cY < 20 ~ 2,
                                                  cY >= 20 & cY < 30 ~ 3,
                                                  cY >= 30 & cY < 40 ~ 4,
                                                  cY >= 40 & cY < 50 ~ 5,
                                                  cY >= 50 & cY < 60 ~ 6,
                                                  cY >= 60 & cY < 70 ~ 7,
                                                  cY >= 70 & cY < 80 ~ 8,
                                                  cY >= 80 & cY < 90 ~ 9,
                                                  cY >= 90 & cY < 100 ~ 10,
                                                  cY >= 100 & cY < 110 ~ 11,
                                                  cY >= 110 & cY < 120 ~ 12,
                                                  cY >= 120 & cY < 130 ~ 13,
                                                  cY >= 130 & cY < 140 ~ 14,
                                                  cY >= 140 & cY < 150 ~ 15,
                                                  cY >= 150 & cY < 160 ~ 16,
                                                  cY >= 160 & cY < 170 ~ 17,
                                                  cY >= 170 & cY < 180 ~ 18,
                                                  cY >= 180 & cY < 190 ~ 19,
                                                  cY >= 190 & cY <= 200 ~ 20))


#Draw 5 10x10m plots from the 10x200 plot (without replacement) and iterate with purrr
#Use the dataset prepared above and do one site at a time
num_reps = 1000

#boot_mod_BC <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'BC') %>% 
                              #  mutate(boot = .x)) %>%  data.frame()
#boot_mod_BH <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'BH') %>% 
                             #   mutate(boot = .x)) %>%  data.frame()
#boot_mod_OP <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'OP') %>% 
                               # mutate(boot = .x)) %>%  data.frame()
#boot_mod_PM <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'PM') %>% 
                               # mutate(boot = .x)) %>%  data.frame()
#boot_mod_IB <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'IB') %>% 
                              #  mutate(boot = .x)) %>%  data.frame()
#boot_mod_WP2 <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'WP2') %>% 
                              #   mutate(boot = .x)) %>%  data.frame()
#boot_mod_BM <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'BM') %>% 
                               # mutate(boot = .x)) %>%  data.frame()
#boot_mod_all <- rbind(boot_mod_BH, boot_mod_BM, boot_mod_WP2, boot_mod_IB, 
                      #boot_mod_PM, boot_mod_OP, boot_mod_BC)

# keep data in case want to do something else with the distribution
#write.csv(boot_mod_all, './tables/Boot_data_all_sites.csv', row.names = FALSE) 
#using pre-run data set so as not to rerun everytime
boot_mod_all <- read.csv('C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj/tables/Boot_data_all_sites.csv')
#calculate 2.5% and 97.5% confidence intervals

#BC_CI <- bootCI(boot_mod_BC) #function written to calculate confidence intervals on individual site dataframes

boot_mod_all$Site <- ordered(boot_mod_all$Site,
                             levels = c("BH", "OP", "BM", "PM", "BC", "IB", "WP2"))

Box_den <- ggplot(data = boot_mod_all, aes(x = Site, y = density_ha))+
  geom_violin()+
  geom_boxplot(width=0.2, color="grey", alpha=0.2) +
  geom_point(data = SumTable_ha %>% filter(Site != "WP"), 
             aes(x = Site, y = density_ha, shape = SampleEventNum, fill = SampleEventNum))+
  scale_shape_manual(values = c("1" = 24, "2" = 21), 
                     name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'))+
  scale_fill_manual(values = c("1" = '#a1d99b', "2" = '#31a354'), 
                     name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'))+
  labs(y = "stems/ha")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_FHM()+
  ggtitle("Tree density")

Box_den

names(boot_mod_all)

Box_BA<- ggplot(data = boot_mod_all, aes(x = Site, y = BA_m2ha))+
  geom_violin()+
  geom_boxplot(width=0.2, color="grey", alpha=0.2) +
  geom_point(data = SumTable_ha %>% filter(Site != "WP"), 
             aes(x = Site, y = BA_m2ha, shape = SampleEventNum, fill = SampleEventNum))+
  scale_shape_manual(values = c("1" = 24, "2" = 21), 
                     name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'))+
  scale_fill_manual(values = c("1" = '#a1d99b', "2" = '#31a354'), 
                     name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'))+
  xlab('Site')+
  scale_x_discrete(labels = site_names)+
  ylab(bquote('Basal area ('~m^2*'/ha)'))+ 
  ggtitle("Tree basal area")+
  theme_FHM()

Box_BA

Box_bio<- ggplot(data = boot_mod_all, aes(x = Site, y = biomass_Mgha))+
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_point(data = SumTable_ha %>% filter(Site != "WP"), 
             aes(x = Site, y = biomass_Mgha, shape = SampleEventNum))+
  scale_shape_manual(name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = 2, "2" = 1))+
  
  xlab('Site')+
  ylab(bquote('Biomass Mg/ha)'))+ 
  theme_FHM()

Box_bio

#Combine as one figure using the patchwork package
Box_den_ba <- Box_den/Box_BA +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')
Box_den_ba


# NETN Plot Data ----------------------------------------------------------
importData()
#all live trees in the most recent sampling in ACAD
ANPtrees<- joinTreeData(park = "ACAD", from = 2019, QAQC = FALSE, status = "live", output = "verbose")
names(ANPtrees)

#Convert to frequency and basal area by species from the plot level to per hectare. ACAD plot is 225 m2.
ANPtrees_ha <- ANPtrees %>% dplyr::select(PlotCode, SampleYear, ScientificName, num_stems, BA_cm2) %>% 
                            group_by(PlotCode, ScientificName, SampleYear) %>% 
                            summarise(sum_stems = sum(num_stems),
                                      sum_BA_cm2 = sum(BA_cm2),
                                      SampleYear = first(SampleYear),
                                      PlotCode = first(PlotCode)) %>% 
                            mutate(num_stems_ha = sum_stems * (10000/225), 
                                   BA_m2ha = sum_BA_cm2/225) %>% 
                            dplyr::select(PlotCode, SampleYear, 
                                   ScientificName, num_stems_ha, BA_m2ha) 

#Calculating relative basal area, frequency, and importance value for each site by species
ord_ANP <- ANPtrees_ha %>% group_by(PlotCode) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                         Tot_BA = sum(BA_m2ha))
ord_ANP2 <- left_join(ANPtrees_ha, ord_ANP, by = c("PlotCode"))

ord_ANP3 <- ord_ANP2 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
                             mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
                             mutate(Imp_val = (Rel_freq + Rel_BA)/2)


#Selecting NETN spruce plots
#Criteria: >.5 combined importance value for PIRU, PIGL, and Picea sp on plot. PIMA not included.
#Matches Davis criteria for choosing spruce stands: >50% spruce (not black spruce)
table(ord_ANP3$ScientificName)

ANPplots <- ord_ANP3 %>% filter(ScientificName == "Picea rubens" | ScientificName == "Picea glauca" | ScientificName == "Picea") %>% 
                         group_by(PlotCode) %>% summarise(picIV = sum(Imp_val)) %>% 
                         filter(picIV >= .5)

ANPsprucePlots <- ANPplots$PlotCode

#Converting scientific names to 4 letter codes
table(ord_ANP3$ScientificName)

ord_ANP3.5 <- ord_ANP3 %>% add_column(Species = 'blank', .after = 'ScientificName') %>%
  mutate(Species = case_when(ScientificName == 'Abies balsamea' ~ 'ABBA',
                             ScientificName == 'Acer pensylvanicum' ~ 'ACPE',
                             ScientificName == 'Acer rubrum' ~ 'ACRU',
                             ScientificName == 'Acer saccharum' ~ 'ACSA3',
                             ScientificName == 'Amelanchier' ~ 'AMELANCHIER',
                             ScientificName == 'Betula alleghaniensis' ~ 'BEAL',
                             ScientificName == 'Betula cordifolia' ~ 'BECO',
                             ScientificName == 'Betula papyrifera' ~ 'BEPA',
                             ScientificName == 'Betula populifolia' ~ 'BEPO',
                             ScientificName == 'Betula X cearulea' ~ 'BEXC',
                             ScientificName == 'Fagus grandifolia' ~ 'FAGA',
                             ScientificName == 'Fraxinus americana' ~ 'FRAM',
                             ScientificName == 'Fraxinus pennsylvanica' ~ 'FRPE',
                             ScientificName == 'Larix laricina' ~ 'LALA',
                             ScientificName == 'Picea glauca' ~ 'PIGL',
                             ScientificName == 'Picea mariana' ~ 'PIMA',
                             ScientificName == 'Picea rubens' ~ 'PIRU',
                             ScientificName == 'Pinus banksiana' ~ 'PIBA',
                             ScientificName == 'Pinus resinosa' ~ 'PIRE',
                             ScientificName == 'Pinus rigida' ~ 'PIRI',
                             ScientificName == 'Pinus strobus' ~ 'PIST',
                             ScientificName == 'Populus grandidentata' ~ 'POGR',
                             ScientificName == 'Quercus rubra' ~ 'QURU',
                             ScientificName == 'Thuja occidentalis' ~ 'THOC',
                             ScientificName == 'Tsuga canadensis' ~ 'TSCA')) 

#Getting physiographic class for NETN plots
ACAD_events <- joinLocEvent(park = "ACAD", from = 2019, QAQC = FALSE)
ACAD_events2 <- ACAD_events %>% dplyr::select(PlotCode, PhysiographySummary)

ACAD_phys <- left_join(ANPplots, ACAD_events2, by = "PlotCode")
table(ACAD_phys$PhysiographySummary)

# Ordination --------------------------------------------------------------
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
densz_mMDS_plot

#ggsave("./figures/Ord_den_sizeclass.jpg", densz_mMDS_plot, height = 5, width = 9, dpi = 400)

# Species composition ordinations: trees + saps -----------------------------------------
#use tree_sp_ha created earlier: all live stems >2.5cm DBH summed by species, site, and sample event. 
#Make matrix with species on one axis and site+sampling event on the other

ord_tree_sp_ha <- Comb_tree_event %>% #filter(size_class != "d2.5_9.9") %>% 
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
#tried removing saplings and the ordination wouldn't run: too little data

table(Comb_tree_event$size_class)

names(ord_tree_sp_ha)
ord_trees <- ord_tree_sp_ha %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
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

#Calculating relative basal area, frequency, and importance value for each site by species
ord_trees2 <- ord_trees %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                         Tot_BA = sum(BA_m2ha))
ord_trees3 <- left_join(ord_trees, ord_trees2, by = c("Site", "SampleEventNum"))

ord_trees4 <- ord_trees3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2) %>% ungroup() %>% 
  droplevels()

# Importance Value all species------------------------------------------------------
ord_IV_wide <- ord_trees4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

str(ord_trees4)
#ordination
mMDS_IVallz <- metaMDS(ord_IV_wide, distance="bray", k=2)
ordiplot(mMDS_IVallz)
mMDS_IVallz$stress

#Plotting
mMDS_IVallz_gg <- fortify(mMDS_IVallz) #transforms ordination results form ggplot can use

mMDS_IVallz_gg2 <- mMDS_IVallz_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_IVallz_gg2$SampleEventNum <- as.factor(mMDS_IVallz_gg2$SampleEventNum)
mMDS_IVallz_gg2$SiteName <- as.factor(mMDS_IVallz_gg2$SiteName)
str(mMDS_IVallz_gg2)

IVallz_mMDS_plot <- mMDS_IVallz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
IVallz_mMDS_plot

#ggsave("./figures/Ord_sp_comp_all_wsaps.jpg", IVallz_mMDS_plot, height = 5, width = 9, dpi = 400)


# Importance value: least common species dropped ---------------------------
#rank species by summed relative density across all sites+events
IV_sp_rank <- ord_IV_wide %>% map_dfc(mean) %>%
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_IV")
print(arrange(IV_sp_rank, desc(Rel_IV)))
# Only including species with >1% average relative IV: PIRU, ABBA, PIGL, TSCA, BEPA, ACRU
ord_IV_wide_t6 <- ord_IV_wide %>%  dplyr::select (c(PIRU, ABBA, PIGL, TSCA, BEPA, ACRU))

#ordination
mMDS_IVt6z <- metaMDS(ord_IV_wide_t6, distance="bray", k=2)
ordiplot(mMDS_IVt6z)
mMDS_IVt6z$stress

#Plotting
mMDS_IVt6z_gg <- fortify(mMDS_IVt6z) #transforms ordination results form ggplot can use

mMDS_IVt6z_gg2 <- mMDS_IVt6z_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_IVt6z_gg2$SampleEventNum <- as.factor(mMDS_IVt6z_gg2$SampleEventNum)
mMDS_IVt6z_gg2$SiteName <- as.factor(mMDS_IVt6z_gg2$SiteName)

IVt6_mMDS_plot <- mMDS_IVt6z_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
IVt6_mMDS_plot

# ggsave("./figures/Ord_sp_comp_top6_wsaps.jpg", IVt6_mMDS_plot, height = 5, width = 9, dpi = 400)


# Importance value: uncommon species pooled -------------------------------
print(arrange(IV_sp_rank, desc(Rel_IV))) # pool all species <1% as other conifer or other hardwood

ord_trees5 <- ord_trees %>% mutate(Species = case_when(Species == 'THOC' ~ "OCON",
                                                        Species == 'ACPE' ~ "OHWD",
                                                        Species == 'AMELANCHIER' ~ "OHWD",
                                                        Species == 'PIST' ~ "OCON",
                                                        Species == 'BECO' ~ "OHWD",
                                                        Species == 'BEAL' ~ "OHWD",
                                                        Species == 'ACSP' ~ "OHWD",
                                                        Species == 'SODE' ~ "OHWD",
                                                        Species == 'PRPE' ~ "OHWD",
                                                        TRUE ~ Species))
table(ord_trees5$Species)
ord_trees5.5 <- ord_trees5 %>% group_by(Site, SampleEventNum, Species) %>% summarise(num_stems_ha = sum(num_stems_ha),  
                                                                                     BA_m2ha = sum(BA_m2ha),
                                                                                     OrdSite = first(OrdSite))

ord_trees6 <- ord_trees5.5 %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                         Tot_BA = sum(BA_m2ha))

ord_trees7 <- left_join(ord_trees5.5, ord_trees2, by = c("Site", "SampleEventNum"))

ord_trees8 <- ord_trees7 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
                             mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
                             mutate(Imp_val = (Rel_freq + Rel_BA)/2) %>% ungroup() %>% 
                             droplevels()

str(ord_trees8)
ord_IV_hrd_wide <- ord_trees8 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
                                  arrange(OrdSite) %>% 
                                  pivot_wider(names_from = Species,
                                              values_from = Imp_val,
                                              values_fill = 0) %>% 
                                  column_to_rownames(var = "OrdSite")

#ordination
mMDS_IVhrdz <- metaMDS(ord_IV_hrd_wide, distance="bray", k=2)
ordiplot(mMDS_IVhrdz)
mMDS_IVhrdz$stress

#Plotting
mMDS_IVhrdz_gg <- fortify(mMDS_IVhrdz) #transforms ordination results form ggplot can use

mMDS_IVhrdz_gg2 <- mMDS_IVhrdz_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_IVhrdz_gg2$SampleEventNum <- as.factor(mMDS_IVhrdz_gg2$SampleEventNum)
mMDS_IVhrdz_gg2$SiteName <- as.factor(mMDS_IVhrdz_gg2$SiteName)

IVhrd_mMDS_plot <- mMDS_IVhrdz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
IVhrd_mMDS_plot

#ggsave("./figures/Ord_sp_comp_hrdpooled_wsaps.jpg", IVhrd_mMDS_plot, height = 5, width = 9, dpi = 400)


# Basal area all species------------------------------------------------------
ord_BA_wide <- ord_trees4 %>% ungroup() %>% dplyr::select(OrdSite, Species, Rel_BA) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Rel_BA,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

#ordination
mMDS_BAallz <- metaMDS(ord_BA_wide, distance="bray", k=2)
ordiplot(mMDS_BAallz)
mMDS_BAallz$stress

#Plotting
mMDS_BAallz_gg <- fortify(mMDS_BAallz) #transforms ordination results form ggplot can use

mMDS_BAallz_gg2 <- mMDS_BAallz_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_BAallz_gg2$SampleEventNum <- as.factor(mMDS_BAallz_gg2$SampleEventNum)
mMDS_BAallz_gg2$SiteName <- as.factor(mMDS_BAallz_gg2$SiteName)
str(mMDS_BAallz_gg2)

BAallz_mMDS_plot <- mMDS_BAallz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
BAallz_mMDS_plot

ggsave("./figures/Ord_sp_comp_all_wsapsBA.jpg", BAallz_mMDS_plot, height = 5, width = 9, dpi = 400)


# Basal area: least common species dropped ---------------------------
#rank species by summed relative density across all sites+events
BA_sp_rank <- ord_BA_wide %>% map_dfc(mean) %>%
  pivot_longer(cols = everything(), names_to = "Species",
               values_to = "Rel_BA")
print(arrange(BA_sp_rank, desc(Rel_BA)))
# Only including species with >1% average relative BA: PIRU, ABBA, PIGL, TSCA, BEPA, ACRU
ord_BA_wide_t6 <- ord_BA_wide %>%  dplyr::select (c(PIRU, ABBA, PIGL, TSCA, BEPA, ACRU))

#ordination
mMDS_BAt6z <- metaMDS(ord_BA_wide_t6, distance="bray", k=2)
ordiplot(mMDS_BAt6z)
mMDS_BAt6z$stress

#Plotting
mMDS_BAt6z_gg <- fortify(mMDS_BAt6z) #transforms ordination results form ggplot can use

mMDS_BAt6z_gg2 <- mMDS_BAt6z_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_BAt6z_gg2$SampleEventNum <- as.factor(mMDS_BAt6z_gg2$SampleEventNum)
mMDS_BAt6z_gg2$SiteName <- as.factor(mMDS_BAt6z_gg2$SiteName)

BAt6_mMDS_plot <- mMDS_BAt6z_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
BAt6_mMDS_plot

 #ggsave("./figures/Ord_sp_comp_top6_wsapsBA.jpg", BAt6_mMDS_plot, height = 5, width = 9, dpi = 400)


# Basal area: uncommon species pooled -------------------------------
print(arrange(BA_sp_rank, desc(Rel_BA))) # pool all species <1% as other conifer or other hardwood

ord_trees5 <- ord_trees %>% mutate(Species = case_when(Species == 'THOC' ~ "OCON",
                                                       Species == 'ACPE' ~ "OHWD",
                                                       Species == 'AMELANCHIER' ~ "OHWD",
                                                       Species == 'PIST' ~ "OCON",
                                                       Species == 'BECO' ~ "OHWD",
                                                       Species == 'BEAL' ~ "OHWD",
                                                       Species == 'ACSP' ~ "OHWD",
                                                       Species == 'SODE' ~ "OHWD",
                                                       Species == 'PRPE' ~ "OHWD",
                                                       TRUE ~ Species))

table(ord_trees5$Species)
ord_trees5.5 <- ord_trees5 %>% group_by(Site, SampleEventNum, Species) %>% summarise(num_stems_ha = sum(num_stems_ha),  
                                                                                     BA_m2ha = sum(BA_m2ha),
                                                                                     OrdSite = first(OrdSite))

ord_trees6 <- ord_trees5.5 %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                            Tot_BA = sum(BA_m2ha))

ord_trees7 <- left_join(ord_trees5.5, ord_trees2, by = c("Site", "SampleEventNum"))

ord_trees8 <- ord_trees7 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2) %>% ungroup() %>% 
  droplevels()

str(ord_trees8)
ord_BA_hrd_wide <- ord_trees8 %>% ungroup() %>% dplyr::select(OrdSite, Species, Imp_val) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = Species,
              values_from = Imp_val,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

#ordination
mMDS_BAhrdz <- metaMDS(ord_BA_hrd_wide, distance="bray", k=2)
ordiplot(mMDS_BAhrdz)
mMDS_BAhrdz$stress

#Plotting
mMDS_BAhrdz_gg <- fortify(mMDS_BAhrdz) #transforms ordination results form ggplot can use

mMDS_BAhrdz_gg2 <- mMDS_BAhrdz_gg %>% filter(Score == "sites") %>% 
  mutate(SiteName = str_sub(Label, 1,2)) %>% 
  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                    SampleEventNum <50 ~ 2)) 

mMDS_BAhrdz_gg2$SampleEventNum <- as.factor(mMDS_BAhrdz_gg2$SampleEventNum)
mMDS_BAhrdz_gg2$SiteName <- as.factor(mMDS_BAhrdz_gg2$SiteName)

BAhrd_mMDS_plot <- mMDS_BAhrdz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label, color = SampleEventNum))+
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
  theme(legend.position = c(.8,0.97),
        legend.justification = c(0.03,.97))+
  theme_FHM()
BAhrd_mMDS_plot

ggsave("./figures/Ord_sp_comp_hrdpooled_wsapsBA.jpg", BAhrd_mMDS_plot, height = 5, width = 9, dpi = 400)
