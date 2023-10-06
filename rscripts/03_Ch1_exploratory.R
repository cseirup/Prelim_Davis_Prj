
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
nicknames <- c("cwd", "events", "qd_ch", "seeds", "qd_sp", "rwlong", "saps", "seeds59", "soil_d",
               "cores", "trees59", "trees59qmd", "trees","coreIDtrees", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
#list2env(data, envir = .GlobalEnv)

# Common plot labels ------------------------------------------------------
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 
                'PM' = "Pemetic Mtn", 'BH' = 'Bass Harbor Head', 
                'IB' = 'Ironbound Island', 'WP' = 'Western Mtn',
                'WP2' = 'Western Mtn 2')

sp_names <- c('PIRU' = 'Picea rubens', 'PIGL' = 'Picea glauca', 'TSCA' = 'Tsuga canadensis', 
              'ABBA' = "Abies balsamea", 'ACRU' = 'Acer rubrum', 'BEPA' = 'Betula papyrifera',
              'PIST' = 'Pinus strobus', 'BEAL' = 'Betula alleghaniensis', 'THOC' = 'Thuja occidentalis',
              'BECO' = 'Betula cordifolia', 'ACPE' = 'Acer pensylvanicum', 'AMELANCHIER' = 'Amelanchier spp', 
              'ACSP' = 'Acer spicatum')

# Stem map ----------------------------------------------------------------
list2env(data["trees"], envir = .GlobalEnv)
trees$Site <- ordered(trees$Site,
                      levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP", "WP2"))

stem_mapV <- trees %>% filter(Site != 'WP') %>% 
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

stem_mapWP <- trees %>% filter(Site == 'WP') %>% 
  ggplot(aes(x = cY, y = cX))+
  geom_point(aes(color = Status, size = DBH))+
  #coord_equal()+# makes x and y the same scale
  scale_color_manual(name = "Status", labels = c("Dead", "Live"), 
                     values = c("D" = '#808080', "L" = '#31a354'))+
  scale_size_continuous(range = c(1, 5))+
  facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names))+
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
        aspect.ratio = 3/10,
        legend.position="bottom")+
  geom_hline(yintercept = 5, linetype = 2)+
  geom_hline(yintercept = 10, linetype = 2)+
  geom_hline(yintercept = 15, linetype = 2)+
  #geom_text(aes(label=Tag), size = 3)+ #add tag number labels
  #coord_flip()+
  theme_FHM()
stem_mapWP
# Sapling Heat Maps -------------------------------------------------------
list2env(data["saps"], envir = .GlobalEnv)
saps$Site <- ordered(saps$Site,
                     levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP", "WP2"))

sapsb <- saps %>% mutate(X = case_when(Transect == 1|Transect == 3 ~ 2.5, 
                                       Transect == 2|Transect == 4 ~ 7.5)) %>%  
  mutate(Y = case_when(Transect == 3|Transect == 4 ~ Quadrat+100,
                       Transect == 1|Transect == 2 ~ Quadrat+0)) %>%
  drop_na() %>% 
  mutate(num_stems = 1) #drop_na deletes quadrats with no saplings recorded

#summarize by site and remove WP (oddly shaped plot, would have to be handled seperately)
sapsc <- sapsb %>% group_by(Site, Transect, Quadrat) %>% summarise(sum_stems = sum(num_stems),
                                                                   X = first(X),
                                                                   Y = first(Y)) %>% 
  filter(Site != 'WP')

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
sapsd <- sapsc %>% rename(cX = X, cY = Y) %>% select(-Quadrat) #renameing columns to match trees
treesB <- trees %>% select(Site, Transect, Status, DBH, cX, cY) %>%  filter(Site!='WP')

treesap_map <- ggplot(treesB, aes(x = cY, y = cX))+
  geom_raster(data = sapsd, aes(fill = sum_stems), hjust = 1, vjust = 0.5)+
  geom_hline(yintercept = 5, linetype = 2)+
  geom_vline(xintercept = 100, linetype = 2)+ 
  scale_fill_gradient(name = "Sapling stems", low = "white", high = "#3182bd")+
  geom_point(aes(color = Status, size = DBH))+
  scale_color_manual(name = "Tree Status", labels = c("Dead", "Live"), 
                     values = c("D" = '#808080', "L" = '#31a354'))+
  scale_size_continuous(name = 'Tree DBH', range = c(.5, 3.5))+
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

treesap_map

# Combine 2020/21 trees and saps to match Davis size classes---------------------------
#using Davis size classes, only live trees
list2env(data["saps"], envir = .GlobalEnv)
saps2 <- saps %>% drop_na() %>% filter(DBH>=2.5)#only taking >2.5cm to match Davis size classes
treesL_df1 <- filter(trees, Status == "L") %>% droplevels()  #only live trees

treesL_df1.5 <- treesL_df1 %>% select(Site, Transect, Species, DBH, SampleYear, SampleEventNum)
saps3 <- saps2 %>% select(Site, Transect, Species, DBH, SampleYear, SampleEventNum) 
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

trees59qmd2 <- trees59qmd %>% select("Site", "Species", "size_class", "num_stem",  
                                     "QMD", "midpoint", "SampleYear", "SampleEventNum")

trees59qmd2$QMD <- coalesce(trees59qmd2$QMD, trees59qmd2$midpoint) #replacing missing QMDs with midpoints for that size class/species
trees59qmd2$BA_cm2 <- round(pi * (((trees59qmd2$QMD/2)^2)*trees59qmd2$num_stem), 4) #basal area per site

# Combine 1959 and 2020/21 >2.4 tree data ---------------------------------
intersect(names(trees59qmd2), names(Ltreesap2))

Ltreesap2.5 <- Ltreesap2 %>% select(Site, Species, SampleYear, SampleEventNum, 
                                    size_class, num_stem, BA_cm2, DBH) %>% 
  rename(DBH_QMD = DBH)

trees59qmd2.5 <- trees59qmd2 %>% select(Site, Species, SampleYear, SampleEventNum, 
                                        size_class, num_stem, BA_cm2, QMD) %>% 
  rename(DBH_QMD = QMD)

Comb_tree <- rbind(Ltreesap2.5, trees59qmd2.5)

#Duplicate 1959 WP data so it will plot against 2022 WP2 data
Comb_tree2 <- Comb_tree %>% filter(SampleYear == 1959 & Site == "WP")
Comb_tree2$Site <- recode(Comb_tree2$Site, WP = "WP2")
#add back to the full tree data
Comb_tree3 <- rbind(Comb_tree, Comb_tree2)

# Combine with event data
list2env(data["events"], envir = .GlobalEnv)
event_tree <- events %>% filter(Module == "Trees")

intersect(names(Comb_tree3), names(event_tree))

Comb_tree_event <- left_join(Comb_tree3, event_tree, by = c("Site","SampleYear","SampleEventNum"))

sort(unique(Comb_tree_event$size_class))
Comb_tree_event$size_class <- ordered(Comb_tree_event$size_class,
                                      levels = c("d2.5_9.9", "d10_19.9", "d20_29.9", 
                                                 "d30_39.9", "d40_49.9", "d50_59.9", 
                                                 "d60_69.9", "d70_79.9", "d80_89.9", 
                                                 "d90_99.9", "d100p"))
levels(Comb_tree_event$size_class) #Size classes are in the order we want them now

Comb_tree_event$Site <- ordered(Comb_tree_event$Site,
                                levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP", "WP2"))
levels(Comb_tree_event$Site)#Sites are in the ordered by complexity now

Comb_tree_event$SampleEventNum <- as.character(Comb_tree_event$SampleEventNum)
# Calculate basal area and density by site and size class for diameter distributions---------------------------
tree_dist_ha <- Comb_tree_event %>% select(-c(NumSubplots, SubplotArea, Module, Species, DBH_QMD)) %>% 
  group_by(Site, size_class, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = round(sum_stems * (10000/TotArea), digits = 0), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  mutate(log10stems_ha = log10(num_stems_ha)) %>% 
  select(Site, SiteName, SampleEventNum, SampleYear, 
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
        legend.title = element_text(size = 12))+
        #legend.position = c(1,0),
        #legend.justification = c(1,0))+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_discrete(labels= c('2.5-10','10-20',
                             '20-30','30-40','40-50',
                             '50-60','60-70','70+'))+ 
  theme_FHM() 

Comp_tree_dist_plot  

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

#staggered bar plots
Comp_tree_dist_bar <- ggplot(data = tree_dist_ha, aes(x = size_class, y = num_stems_ha, fill = SampleEventNum))+
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
        legend.position = c(0.97,0.3),
        legend.justification = c(1,0))+
  scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020/21'), values = c("1" = '#a1d99b', "2" = '#31a354'))+
  scale_x_discrete(labels= c('2.5 \U2013 10','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
  theme_FHM() 

Comp_tree_dist_bar

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
SumTable_site <- SumTable %>% select(-c(NumSubplots, SubplotArea, Module, Species, size_class, DBH_QMD)) %>% 
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
  select(Site, SiteName, SampleEventNum, 
         density_ha, BA_m2ha, biomass_Mgha, carbonmass_Mgha) 

SumTable_ha$SampleEventNum <- as.character(SumTable_ha$SampleEventNum)

SumTable_ha2 <- SumTable_ha %>% pivot_wider(names_from = SampleEventNum, 
                                            values_from = c(density_ha, BA_m2ha, 
                                                            biomass_Mgha, carbonmass_Mgha)) 

#Adding sapling stems/ha for both visits to Table 2
sapTable <- Comb_tree_event %>% filter(size_class == 'd2.5_9.9') %>% 
  select(-c(NumSubplots, SubplotArea, Module, Species, size_class, DBH_QMD)) %>% 
  group_by(Site, SampleEventNum) %>% 
  summarise(SampleYear = first(SampleYear),
            SiteName = first(SiteName),
            sum_stems = sum(num_stem),
            TotArea = first(TotArea)) 

#converting from plot to ha and kg to megagram 
sapTable_ha <- sapTable %>%  mutate(density_ha = sum_stems * (10000/TotArea)) %>%  
                              mutate(across(where(is.numeric), round, 0)) %>% 
                              select(Site, SiteName, SampleEventNum, 
                                                density_ha) 

sapTable_ha$SampleEventNum <- as.character(sapTable_ha$SampleEventNum)

sapTable_ha2 <- sapTable_ha %>% pivot_wider(names_from = SampleEventNum, 
                                            values_from = c(density_ha)) %>% 
                                rename("sap_den_1"= "1") %>% 
                                rename("sap_den_2"= "2")



SumTable_ha3 <- left_join(SumTable_ha2, sapTable_ha2, by = c("Site", "SiteName"))

#Adding standing dead tree BA to table
names(events)
names(trees)
dtrees <- trees %>% filter(Status == "D") %>% left_join(event_tree, by = c("Site", "SampleEventNum", "SampleYear"))
dtrees$BA_cm2 <- round(pi * ((dtrees$DBH/2)^2), 4)
names(dtrees)
dtreeTable <- dtrees %>%  select(Site,SiteName,BA_cm2,SampleEventNum,SampleYear,TotArea) %>% 
  group_by(Site, SampleEventNum) %>% 
  summarise(SampleYear = first(SampleYear),
            SiteName = first(SiteName),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea)) 

#converting from plot to ha 
dtreeTable_ha <- dtreeTable %>%  mutate(BA_m2ha = sum_BA_cm2/TotArea) %>%  
  mutate(across(where(is.numeric), round, 0)) %>% 
  select(Site, SiteName,  BA_m2ha) 

dtreeTable_ha2 <- dtreeTable_ha %>% rename("2020s SDW BA_m2ha"= BA_m2ha)


SumTable_ha4 <- left_join(SumTable_ha3, dtreeTable_ha2, by = c("Site", "SiteName"))

#Adding CWD
list2env(data["cwd"], envir = .GlobalEnv)

SumTable_ha5 <- left_join(SumTable_ha4, cwd, by = "Site")

SumTable_ha6 <- SumTable_ha5 %>% ungroup() %>% 
  mutate(vol_m3ha = round(vol_m3ha)) %>% 
  select(-c("carbonmass_Mgha_2", "carbonmass_Mgha_1", "num_pieces","Site"))

kable(SumTable_ha6)
#write.csv(Comb_sum_table, './tables/Summary_metrics1959_2020_2.5.csv', row.names = FALSE)


# Table 1 -----------------------------------------------------------------

dSites <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Site_description_table.xlsx")
kable(dSites)

# Table 3 -----------------------------------------------------------------

sCores <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Core_site_summary_table.xlsx")
kable(sCores)


# Species comparison: Basal area trees and saps, density seeds ----------------------------------------
#Species composition by site: 
#trees
tree_sp_ha <- Comb_tree_event %>% filter(size_class != "d2.5_9.9") %>% 
                                    select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
                                    group_by(Site, Species, SampleEventNum) %>% 
                                    summarise(sum_stems = sum(num_stem),
                                              sum_BA_cm2 = sum(BA_cm2),
                                              TotArea = first(TotArea),
                                              SampleYear = first(SampleYear),
                                              SiteName = first(SiteName)) %>% 
                                    mutate(num_stems_ha = sum_stems * (10000/TotArea), 
                                           BA_m2ha = sum_BA_cm2/TotArea) %>% 
                                    select(Site, SiteName, SampleEventNum, SampleYear, 
                                           Species, num_stems_ha, BA_m2ha) 

#Fill in zeros for missing species detections (for better plotting)
tree_sp_ha2 <- tree_sp_ha %>% select(-num_stems_ha) %>% pivot_wider(names_from = Species, values_from = BA_m2ha, values_fill = 0)
tree_sp_haC <- tree_sp_ha2 %>% pivot_longer(cols = ABBA:BECO, names_to = "Species", values_to = "BA_m2ha")

#saps
sap_sp_ha <- Comb_tree_event %>% filter(size_class == 'd2.5_9.9') %>% 
                                  select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
                                  group_by(Site, Species, SampleEventNum) %>% 
                                  summarise(sum_stems = sum(num_stem),
                                            sum_BA_cm2 = sum(BA_cm2),
                                            TotArea = first(TotArea),
                                            SampleYear = first(SampleYear),
                                            SiteName = first(SiteName)) %>% 
                                  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
                                         BA_m2ha = sum_BA_cm2/TotArea) %>% 
                                  select(Site, SiteName, SampleEventNum, SampleYear, 
                                         Species, num_stems_ha, BA_m2ha) 
#Fill in zeros for missing species detentions (for better plotting)
sap_sp_ha2 <- sap_sp_ha %>% select(-num_stems_ha) %>% pivot_wider(names_from = Species, values_from = BA_m2ha, values_fill = 0)
sap_sp_haC <- sap_sp_ha2 %>% pivot_longer(cols = ABBA:ACSP, names_to = "Species", values_to = "BA_m2ha")

#seeds
list2env(data["seeds"], envir = .GlobalEnv)
list2env(data["seeds59"], envir = .GlobalEnv)

seeds2 <- seeds %>% filter(Latin_name != "No species") %>%
                    group_by(Site, Latin_name) %>% 
                    summarise(Sum_stem = sum(Count),
                              SampleYear = first(SampleYear),
                              SampleEventNum = first(SampleEventNum)) %>% 
                    mutate(den_m2 = Sum_stem/30) # 30 m2 quadrats, convert to per m2

seeds3 <- seeds2 %>% select(Site, Latin_name, den_m2, SampleYear, SampleEventNum) %>%
                      ungroup() 

seeds4 <- rbind(seeds3, seeds59)
seeds4$SampleEventNum <- as.character(seeds4$SampleEventNum)
#Fill in zeros for missing species detentions (for better plotting)

seeds4w <- seeds4 %>% pivot_wider(names_from = Latin_name, values_from = den_m2, values_fill = 0)
seedsC <- seeds4w %>% pivot_longer(cols = 4:13, names_to = "Species", values_to = "den_m2")

table(seedsC$Species)

#setting species order for all plots for consitency
tree_sp_haC$Species <- ordered(tree_sp_haC$Species,
                               levels = c("PIRU", "PIGL", "ABBA", "TSCA", "AMELANCHIER", 
                                          "BEPA", "ACRU", "ACPE", "THOC", "PIST", "BECO", "BEAL"))
sap_sp_haC$Species <- ordered(sap_sp_haC$Species,
                               levels = c("PIRU", "PIGL", "ABBA", "TSCA", "AMELANCHIER", 
                                          "BEPA", "ACRU", "ACPE", "THOC", "BECO", "ACSP","PRPE","SODE"))

seedsC$Species <- ordered(seedsC$Species,
                              levels = c("Picea rubens", "Picea glauca", "Abies balsamea", "Tsuga canadensis", "Amelanchier", 
                                         "Acer rubrum", "Acer pensylvanicum ", "Pinus strobus", "Sorbus decora"))
###Beech Mtn####: doing each site separately so species selection is customized
#trees
rankBA(tree_sp_haC, "BM")
BM_tr <- c("PIRU", "BEPA", "TSCA", "ACRU", "THOC") # top 4 for each visit, only kicking out ABBA
#saps
rankBA(sap_sp_haC, "BM")
BM_sp <- c("PIRU", "THOC","ABBA","TSCA","ACPE") #top 4 for each visit, doesn't kick anything out
#seeds
rankDEN(seedsC, "BM")
sp_names
BM_sd <- c("Picea rubens", "Abies balsamea","Tsuga canadensis") #top 4 for each visit, doesn't kick anything out

#BM tree plot
pBM_tree <- Comptreeplot(tree_sp_haC, "BM", BM_tr)
pBM_tree

#BM sap plot
pBM_sap <- Compsapplot(sap_sp_haC, "BM", BM_sp)
pBM_sap

#BM seed plot
pBM_seed <- Compseedplot(seedsC, "BM", BM_sd)
pBM_seed

####Pemetic Mtn###
#trees
rankBA(tree_sp_haC, "PM")
PM_tr <- c("PIRU","TSCA", "BEPA", "ACRU", "PIST") # top 4 for each visit, only kicking out ABBA
#saps
rankBA(sap_sp_haC, "PM")
PM_sp <- c("PIRU", "ACRU","TSCA","ACPE") #top 4 for each visit, doesn't kick anything out
#seeds
rankDEN(seedsC, "PM")
sp_names
PM_sd <- c("Picea rubens", "Betula papyrifera","Tsuga canadensis","Acer pensylvanicum") #top 4 for each visit, doesn't kick anything out

#PM tree plot
pPM_tree <- Comptreeplot(tree_sp_haC, "PM", PM_tr)
pPM_tree

#PM sap plot
pPM_sap <- Compsapplot(sap_sp_haC, "PM", PM_sp)
pPM_sap

#PM seed plot
pPM_seed <- Compseedplot(seedsC, "PM", PM_sd)
pPM_seed

####Pemetic Mtn###
#trees
rankBA(tree_sp_haC, "PM")
PM_tr <- c("PIRU","TSCA", "BEPA", "ACRU", "PIST") # top 4 for each visit, didn't kick any out
#saps
rankBA(sap_sp_haC, "PM")
PM_sp <- c("PIRU", "ACRU","TSCA","ACPE") #top 4 for each visit, doesn't kick anything out
#seeds
rankDEN(seedsC, "PM")
sp_names
PM_sd <- c("Picea rubens", "Betula papyrifera","Tsuga canadensis","Acer pensylvanicum") #top 4 for each visit, doesn't kick anything out

#PM tree plot
pPM_tree <- Comptreeplot(tree_sp_haC, "PM", PM_tr)
pPM_tree

#PM sap plot
pPM_sap <- Compsapplot(sap_sp_haC, "PM", PM_sp)
pPM_sap

#PM seed plot
pPM_seed <- Compseedplot(seedsC, "PM", PM_sd)
pPM_seed

####Blackwoods###
#trees
rankBA(tree_sp_haC, "BC")
BC_tr <- c("PIRU","BEPA", "THOC", "ACRU", "PIST") # top 4 for each visit, didn't kick any out
#saps
rankBA(sap_sp_haC, "BC")
BC_sp <- c("PIRU", "ABBA", "AMELANCHIER", "BEPA", "ACPE") #top 4 for each visit, doesn't kick anything out
#seeds
rankDEN(seedsC, "BC")
sp_names
BC_sd <- c("Picea rubens",  "Abies balsamea","Pinus strobus") #top 4 for each visit, doesn't kick anything out

#BC tree plot
pBC_tree <- Comptreeplot(tree_sp_haC, "BC", BC_tr)
pBC_tree

#BC sap plot
pBC_sap <- Compsapplot(sap_sp_haC, "BC", BC_sp)
pBC_sap

#BC seed plot
pBC_seed <- Compseedplot(seedsC, "BC", BC_sd)
pBC_seed

####Otter Point###
#trees
rankBA(tree_sp_haC, "OP")
OP_tr <- c("PIRU","BEPA", "PIGL", "ABBA") # top 4 for each visit, didn't kick any out
#saps
rankBA(sap_sp_haC, "OP")
OP_sp <- c("PIRU", "ABBA", "AMELANCHIER", "BEPA") #top 4 for each visit, doesn't kick anything out
#seeds
rankDEN(seedsC, "OP")
sp_names
OP_sd <- c("Picea rubens",  "Abies balsamea","Amelanchier") #top 4 for each visit, doesn't kick anything out

#OP tree plot
pOP_tree <- Comptreeplot(tree_sp_haC, "OP", OP_tr)
pOP_tree

#OP sap plot
pOP_sap <- Compsapplot(sap_sp_haC, "OP", OP_sp)
pOP_sap

#OP seed plot
pOP_seed <- Compseedplot(seedsC, "OP", OP_sd)
pOP_seed

####Bass Harbor###
#trees
rankBA(tree_sp_haC, "BH")
BH_tr <- c("PIRU","BEPA", "PIGL", "ABBA") # top 4 for each visit, Kicked out Amelanchier
#saps
rankBA(sap_sp_haC, "BH")
BH_sp <- c("PIRU", "PIGL", "ABBA", "BECO") #top 4 for each visit, kicked out PRPR, BEPA, SODE
#seeds
rankDEN(seedsC, "BH")
sp_names
BH_sd <- c("Picea rubens", "Picea glauca", "Abies balsamea", "Sorbus decora") #top 4 for each visit, doesn't kick anything out

#BH tree plot
pBH_tree <- Comptreeplot(tree_sp_haC, "BH", BH_tr)
pBH_tree

#BH sap plot
pBH_sap <- Compsapplot(sap_sp_haC, "BH", BH_sp)
pBH_sap

#BH seed plot
pBH_seed <- Compseedplot(seedsC, "BH", BH_sd)
pBH_seed

####Ironbound###
#trees
rankBA(tree_sp_haC, "IB")
IB_tr <- c("PIRU","BEAL", "ACRU", "ABBA") # top 4 for each visit, none kicked out
#saps
rankBA(sap_sp_haC, "IB")
IB_sp <- c("PIRU","AMELANCHIER", "ACRU", "ABBA", "ACPE") #top 4 for each visit, kicked out ACSP, THOC
#seeds
rankDEN(seedsC, "IB")
sp_names
IB_sd <- c("Picea rubens", "Acer rubrum", "Abies balsamea", "Pinus strobus") #top 4 for each visit, doesn't kick anything out

#IB tree plot
pIB_tree <- Comptreeplot(tree_sp_haC, "IB", IB_tr)
pIB_tree

#IB sap plot
pIB_sap <- Compsapplot(sap_sp_haC, "IB", IB_sp)
pIB_sap

#IB seed plot
pIB_seed <- Compseedplot(seedsC, "IB", IB_sd)
pIB_seed

####Western Mtn 2###
#trees
rankBA(tree_sp_haC, "WP2")
WP2_tr <- c("PIRU","BEAL", "ACRU", "ABBA", "BECO") # top 4 for each visit, none kicked out
#saps
rankBA(sap_sp_haC, "WP2")
WP2_sp <- c("PIRU","BECO", "ACRU", "ABBA") #top 4 for each visit, kicked out ACSP, THOC
#seeds
rankDEN(seedsC, "WP2")
sp_names
WP2_sd <- c("Picea rubens", "Abies balsamea") #top 4 for each visit, doesn't kick anything out

#WP2 tree plot
pWP2_tree <- Comptreeplot(tree_sp_haC, "WP2", WP2_tr)
pWP2_tree

#WP2 sap plot
pWP2_sap <- Compsapplot(sap_sp_haC, "WP2", WP2_sp)
pWP2_sap

#WP2 seed plot
pWP2_seed <- Compseedplot(seedsC, "WP2", WP2_sd)
pWP2_seed

#Combining all plots with patchwork

spComp <- (pBH_tree+pBH_sap+pBH_seed)/(pOP_tree+pOP_sap+pOP_seed)/(pPM_tree+pPM_sap+pPM_seed)/
          (pBC_tree+pBC_sap+pBC_seed)/(pBM_tree+pBM_sap+pBM_seed)/(pIB_tree+pIB_sap+pIB_seed)/
          (pWP2_tree+pWP2_sap+pWP2_seed)+
          plot_layout(guides = 'collect')

TR <- pBH_tree/pOP_tree/pPM_tree/pBC_tree/pBM_tree/pIB_tree/pWP2_tree
TRy <- wrap_elements(panel = TR) +
       labs(tag = bquote('Basal area ('~m^2*'/ha)')) +
       theme(plot.tag = element_text(size = rel(1.5), angle = 90),
        plot.tag.position = "left")+
  plot_annotation(title = 'Trees', theme = theme(plot.title = element_text(hjust = .5)))
TRy

SP <- pBH_sap/pOP_sap/pPM_sap/pBC_sap/pBM_sap/pIB_sap/pWP2_sap
SPy <- wrap_elements(panel = SP) +
        labs(tag = bquote('Basal area ('~m^2*'/ha)')) +
         theme(plot.tag = element_text(size = rel(1.5), angle = 90),
               plot.tag.position = "left")+
  plot_annotation(title = 'Saplings', theme = theme(plot.title = element_text(hjust = .5)))

SD <- pBH_seed/pOP_seed/pPM_seed/pBC_seed/pBM_seed/pIB_seed/pWP2_seed +
  plot_layout(guides = 'collect')
SDy <- wrap_elements(panel = SD) +
  labs(tag = bquote('Density (stems/'~m^2*')')) +
  theme(plot.tag = element_text(size = rel(1.5), angle = 90),
        plot.tag.position = "left")+
  plot_annotation(title = 'Seedlings', theme = theme(plot.title = element_text(hjust = .5)))

spComp2 <- wrap_elements(TRy)+wrap_elements(SPy)+wrap_elements(SDy)

spComp3 <- wrap_elements(panel = spComp2) +
  labs(tag = "Species") +
theme(plot.tag = element_text(size = rel(1.5)),
      plot.tag.position = "bottom")
spComp3


# Bootstrapping trees -----------------------------------------------------

#remove WP data (first plot, 20 x 100 m) and live trees only to match Davis
names(trees)
trees2 <- trees %>% filter(Site != "WP") %>%  filter(Status == "L") %>% 
  select(-c(X,Y)) #may have to drop level, WP still floating around
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
                             levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP2"))

Box_den <- ggplot(data = boot_mod_all, aes(x = Site, y = density_ha))+
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_point(data = SumTable_ha %>% filter(Site != "WP"), 
             aes(x = Site, y = density_ha, shape = SampleEventNum))+
  scale_shape_manual(name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = 2, "2" = 1))+
  labs(x = "Site", y = "stems/ha")+ 
  theme_FHM()

Box_den

names(boot_mod_all)

Box_BA<- ggplot(data = boot_mod_all, aes(x = Site, y = BA_m2ha))+
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_point(data = SumTable_ha %>% filter(Site != "WP"), 
             aes(x = Site, y = BA_m2ha, shape = SampleEventNum))+
  scale_shape_manual(name = "Non-bootstrapped metric", labels = c("1" = '1959', "2" = '2020-2022'), 
                     values = c("1" = 2, "2" = 1))+
  xlab('Site')+
  ylab(bquote('Basal area ('~m^2*'/ha)'))+ 
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


# NETN Plot Data ----------------------------------------------------------
importData()
#all live trees in the most recent sampling in ACAD
ANPtrees<- joinTreeData(park = "ACAD", from = 2019, QAQC = FALSE, status = "live", output = "verbose")
names(ANPtrees)

#Convert to frequency and basal area by species from the plot level to per hectare. ACAD plot is 225 m2.
ANPtrees_ha <- ANPtrees %>% select(PlotCode, SampleYear, ScientificName, num_stems, BA_cm2) %>% 
                            group_by(PlotCode, ScientificName, SampleYear) %>% 
                            summarise(sum_stems = sum(num_stems),
                                      sum_BA_cm2 = sum(BA_cm2),
                                      SampleYear = first(SampleYear),
                                      PlotCode = first(PlotCode)) %>% 
                            mutate(num_stems_ha = sum_stems * (10000/225), 
                                   BA_m2ha = sum_BA_cm2/225) %>% 
                            select(PlotCode, SampleYear, 
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
ACAD_events2 <- ACAD_events %>% select(PlotCode, PhysiographySummary)

ACAD_phys <- left_join(ANPplots, ACAD_events2, by = "PlotCode")
table(ACAD_phys$PhysiographySummary)

# Ordination --------------------------------------------------------------

# Species composition ordinations -----------------------------------------
#use tree_sp_ha created earlier: all live stems >2.5cm DBH summed by species, site, and sample event. 
#Includes duplicated 1959 WP data named as WP2
#Make matrix with species on one axis and site+sampling event on the other
ord_tree_sp_ha <- Comb_tree_event %>% select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
  group_by(Site, Species, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  select(Site, SiteName, SampleEventNum, SampleYear, 
         Species, num_stems_ha, BA_m2ha) 

names(ord_tree_sp_ha)
ord_trees <- ord_tree_sp_ha %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
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
ord_trees2 <- ord_trees %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                         Tot_BA = sum(BA_m2ha))
ord_trees3 <- left_join(ord_trees, ord_trees2, by = c("Site", "SampleEventNum"))

ord_trees4 <- ord_trees3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2)


# Relative Basal Area -----------------------------------------------------

ord_BA_wide <- ord_trees4 %>% ungroup() %>% select(OrdSite, Species, Rel_BA) %>% 
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
ord_BA_wide2 <- ord_BA_wide %>% select(-c(SODE, PRPE, BECO, AMELANCHIER, ACPE, BEAL, ACSP)) 

sort(names(ord_BA_wide2))

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA

#All species included
BA_pcaALL <- princomp(ord_BA_wide, cor=T)
summary(BA_pcaALL)#gives proportion of variance explained by each PC; 2 gives you 41%, 3 gives you 56%
plot(BA_pcaALL)
biplot(BA_pcaALL)
loadings(BA_pcaALL)

#using ggfortify package + ggplot to get a decent plot
BA_pca_ALL_df <- fortify(BA_pcaALL) %>% select(Comp.1, Comp.2) %>%  rownames_to_column(var = "OrdSite")

#ggrepel package for optimizing labels
options(ggrepel.max.overlaps = Inf)#Prints all labels even if overlap

BA_ord_plotALL <- BA_pca_ALL_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
  theme_FHM()

BA_ord_plotALL

#Top 8 species
BA_pca <- princomp(ord_BA_wide2, cor=T)
summary(BA_pca)#gives proportion of variance explained by each PC; 2 gives you 59%, 3 gives you 74%
plot(BA_pca)
biplot(BA_pca)
loadings(BA_pca)

#using ggfortify package + ggplot to get a decent plot
BA_pca_df <- fortify(BA_pca) %>% select(Comp.1, Comp.2) %>%  rownames_to_column( var = "OrdSite")

#ggrepel package for optimizing labels

BA_ord_plot <- BA_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
  theme_FHM()

BA_ord_plot

# Relative Frequency ------------------------------------------------------
ord_F_wide <- ord_trees4 %>% ungroup() %>% select(OrdSite, Species, Rel_freq) %>% 
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
#ord_BA_wide2 <- ord_BA_wide %>% select(-c(SODE, PRPE, BECO, AMELANCHIER, ACPE, BEAL, ACSP)) 

#PCA: first option in AB hw #7
#Based on slide 41 of wk10 Multivariate, should use PCA
F_pca <- princomp(ord_F_wide, cor=T) # all species included
summary(F_pca)#gives proportion of variance explained by each PC; 2 gives you 47%, 3 gives you 62%
plot(F_pca)
biplot(F_pca)
loadings(F_pca)

#using ggfortify package + ggplot to get a decent plot
F_pca_df <- fortify(F_pca) %>% select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")

F_ord_plot <- F_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
  theme_FHM()

F_ord_plot

# Importance Value ------------------------------------------------------
ord_IV_wide <- ord_trees4 %>% ungroup() %>% select(OrdSite, Species, Imp_val) %>% 
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
IV_pca <- princomp(ord_IV_wide, cor=T) # all species included
summary(IV_pca)#gives proportion of variance explained by each PC; 2 gives you 45%, 3 gives you 59%
plot(IV_pca)
biplot(IV_pca)
loadings(IV_pca)

#using ggfortify package + ggplot to get a decent plot
IV_pca_df <- fortify(IV_pca) %>% select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")

IV_ord_plot <- IV_pca_df %>% ggplot(aes(x=Comp.1, y=Comp.2, label = OrdSite))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(-Inf, NA),
                  ylim = c(-Inf, Inf))+ #labels can overlap plot borders
  theme_FHM()

IV_ord_plot

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
                                          select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
                                          group_by(Site, Species, SampleEventNum) %>% 
                                          summarise(sum_stems = sum(num_stem),
                                                    sum_BA_cm2 = sum(BA_cm2),
                                                    TotArea = first(TotArea),
                                                    SampleYear = first(SampleYear),
                                                    SiteName = first(SiteName)) %>% 
                                          mutate(num_stems_ha = sum_stems * (10000/TotArea), 
                                                 BA_m2ha = sum_BA_cm2/TotArea) %>% 
                                          select(Site, SiteName, SampleEventNum, SampleYear, 
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

ord_IV_sap_wide <- ord_sap_BET4 %>% ungroup() %>% select(OrdSite, Species, Imp_val) %>% 
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
ord_IV_sap_wide2 <- ord_IV_sap_wide %>% select(-c(AMELANCHIER, ACSP, SODE, PRPE)) 

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
sp_ha_BET <- Comb_tree_event_BET %>% select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
  group_by(Site, Species, SampleEventNum) %>% 
  summarise(sum_stems = sum(num_stem),
            sum_BA_cm2 = sum(BA_cm2),
            TotArea = first(TotArea),
            SampleYear = first(SampleYear),
            SiteName = first(SiteName)) %>% 
  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
         BA_m2ha = sum_BA_cm2/TotArea) %>% 
  select(Site, SiteName, SampleEventNum, SampleYear, 
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

ord_IV_all_wide <- ord_all_BET4 %>% ungroup() %>% select(OrdSite, Species, Imp_val) %>% 
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
ord_IV_all_wide2 <- ord_IV_all_wide %>% select(-c(AMELANCHIER, PIST, BEAL, ACSP, SODE, PRPE)) 

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

tree10_sp_ha <- Comb_tree_event %>% filter(size_class != 'd2.5_9.9') %>% select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
                                    group_by(Site, Species, SampleEventNum) %>% 
                                    summarise(sum_stems = sum(num_stem),
                                              sum_BA_cm2 = sum(BA_cm2),
                                              TotArea = first(TotArea),
                                              SampleYear = first(SampleYear),
                                              SiteName = first(SiteName)) %>% 
                                    mutate(num_stems_ha = sum_stems * (10000/TotArea), 
                                           BA_m2ha = sum_BA_cm2/TotArea) %>% 
                                    select(Site, SiteName, SampleEventNum, SampleYear, 
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
ord_10IV_wide <- ord_10trees4 %>% ungroup() %>% select(OrdSite, Species, Imp_val) %>% 
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
IV10_pca_df <- fortify(IV10_pca) %>% select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")

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

#ord_ANP5 <- ord_ANP4 %>% rename(OrdSite =  PlotCode) %>% ungroup() %>% select(-ScientificName)

#ord_10trees5 <- ord_10trees4 %>% ungroup() %>% select(-c(Site, SiteName, SampleEventNum))

#names(ord_10trees5)
#names(ord_ANP5)

#ordNETN <- rbind(ord_10trees5, ord_ANP5)

#ordNETN_wide <- ordNETN %>% ungroup() %>% select(OrdSite, Species, Imp_val) %>% 
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
#NETN_pca_df <- fortify(NETN_pca) %>% select(Comp.1, Comp.2) %>% rownames_to_column(var = "OrdSite")

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
  select(OrdSite, Species, Rel_BA) %>% 
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
ord_10BA_wide2 <- ord_10BA_wide %>% select(-c(AMELANCHIER, ACPE, BECO, BEAL, PIST)) 
sort(names(ord_10BA_wide2))

#metaMDS: Shawn says commonly used for vegetation 
mMDS_BA10 <- metaMDS(ord_10BA_wide2, distance="bray", k=2)# error stress nearly zero
ordiplot(mMDS_BA10)
text(mMDS_BA10)


# Relative Basal Area >2.5cm DBH, splitting by event-----------------------------------------------------
#Basal Area MNDS 1959 sites, top 8 species
ord_BA_wideOLD <- ord_trees4 %>% ungroup() %>% filter(SampleEventNum == 1) %>% select(OrdSite, Species, Rel_BA) %>% 
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
ord_BA_wideOLD2 <- ord_BA_wideOLD %>% select(-c(BEAL, THOC, ACSP)) 

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
#will pool BEPA/BECO

sz_ha_BET <- Comb_tree_event_BET %>% select(-c(NumSubplots, SubplotArea, Module, DBH_QMD)) %>% 
                                         group_by(Site, size_class, SampleEventNum) %>% 
                                         summarise(sum_stems = sum(num_stem),
                                                   sum_BA_cm2 = sum(BA_cm2),
                                                   TotArea = first(TotArea),
                                                   SampleYear = first(SampleYear),
                                                   SiteName = first(SiteName)) %>% 
                                         mutate(num_stems_ha = sum_stems * (10000/TotArea), 
                                                BA_m2ha = sum_BA_cm2/TotArea) %>% 
                                         select(Site, SiteName, SampleEventNum, SampleYear, 
                                                size_class, num_stems_ha, BA_m2ha) 

ord_sz_BET <- sz_ha_BET %>% add_column(OrdSite = 'blank', .after = 'Site') %>%
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
ord_sz_BET2 <- ord_sz_BET %>% group_by(Site, SampleEventNum) %>% summarise(Tot_stems = sum(num_stems_ha),
                                                                             Tot_BA = sum(BA_m2ha))
ord_sz_BET3 <- left_join(ord_sz_BET, ord_sz_BET2, by = c("Site", "SampleEventNum"))

ord_sz_BET4 <- ord_sz_BET3 %>% mutate(Rel_freq = num_stems_ha/Tot_stems) %>% 
  mutate(Rel_BA = BA_m2ha/Tot_BA) %>% 
  mutate(Imp_val = (Rel_freq + Rel_BA)/2)


# Diameter dist. Relative Basal Area NMDS-----------------------------------------------------

ord_BA_sz_wide <- ord_sz_BET4 %>% ungroup() %>% select(OrdSite, size_class, Rel_BA) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = size_class,
              values_from = Rel_BA,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

mMDS_BAsz <- metaMDS(ord_BA_sz_wide, distance="bray", k=2)
ordiplot(mMDS_BAsz)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_BAsz_gg <- fortify(mMDS_BAsz) #transforms ordination results form ggplot can use

mMDS_BAsz_gg2 <- mMDS_BAsz_gg %>% filter(Score == "sites") %>% 
                                  mutate(SiteName = str_sub(Label, 1,2)) %>% 
                                  mutate(SampleEventNum = str_sub(Label, -2)) %>% 
                                  mutate(SampleEventNum = case_when(SampleEventNum > 50 ~ 1,
                                                                    SampleEventNum <50 ~ 2)) 

mMDS_BAsz_gg2$SampleEventNum <- as.factor(mMDS_BAsz_gg2$SampleEventNum)
mMDS_BAsz_gg2$SiteName <- as.factor(mMDS_BAsz_gg2$SiteName)
str(mMDS_BAsz_gg2)
  

BAsz_mMDS_plot <- mMDS_BAsz_gg2 %>% arrange(desc(Label)) %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point(aes(color = SampleEventNum, size = 2), show.legend = FALSE)+
  coord_cartesian(clip = "off") +
  geom_path(aes(group = SiteName), color="black",
            arrow = arrow(type = "closed",
                          length=unit(0.075, "inches")))+
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
BAsz_mMDS_plot

# Diameter dist. Relative Density Area NMDS-----------------------------------------------------

ord_den_sz_wide <- ord_sz_BET4 %>% ungroup() %>% select(OrdSite, size_class, Rel_freq) %>% 
  arrange(OrdSite) %>% 
  pivot_wider(names_from = size_class,
              values_from = Rel_freq,
              values_fill = 0) %>% 
  column_to_rownames(var = "OrdSite")

mMDS_densz <- metaMDS(ord_den_sz_wide, distance="bray", k=2)
ordiplot(mMDS_densz)

#Plotting
library(ggvegan)#this needs to be loaded here in order for the fortify() to work on an NMDS object
mMDS_densz_gg <- fortify(mMDS_densz) #transforms ordination results form ggplot can use

mMDS_densz_gg2 <- mMDS_densz_gg %>% filter(Score == "sites")

densz_mMDS_plot <- mMDS_densz_gg2 %>% ggplot(aes(x=NMDS1, y=NMDS2, label=Label))+
  geom_point()+
  coord_cartesian(clip = "off") +
  geom_text_repel(xlim = c(NA, NA),
                  ylim = c(NA, NA))+ 
  scale_x_continuous(breaks = 1:2, # add buffer to make space for labels
                     expand = expansion(mult = 0.5)) + #labels can overlap plot borders
  theme_FHM()
densz_mMDS_plot

# Diameter dist. importance value  NMDS-----------------------------------------------------

ord_IV_sz_wide <- ord_sz_BET4 %>% ungroup() %>% select(OrdSite, size_class, Imp_val) %>% 
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

ord_IV_sz_wide2 <- ord_sz_BET4 %>% filter(SampleEventNum == 2) %>%  ungroup() %>% select(OrdSite, size_class, Imp_val) %>% 
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

ord_IV_sz_wide1 <- ord_sz_BET4 %>% filter(SampleEventNum == 1) %>%  ungroup() %>% select(OrdSite, size_class, Imp_val) %>% 
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
