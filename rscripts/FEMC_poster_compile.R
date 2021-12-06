#Source code for figures and tables for FEMC poster 12/2021

#install.packages("LaTeX")
#devtools::install_github("katemmiller/forestMIDN")
library(ggalt)
library(tidyverse)
library(forestMIDN)
library(knitr)
library(kableExtra)


# Load data ---------------------------------------------------------------

QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("cwd", "events", "qd_ch", "seeds", "qd_sp", "saps", "soil_d", "trees59", "trees59qmd", "trees", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
list2env(data, envir = .GlobalEnv)

# Common plot labels ------------------------------------------------------
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 
                'PM' = "Pemetic Mtn", 'BH' = 'Bass Harbor Head', 
                'IB' = 'Ironbound Island', 'WP' = 'Western Mtn')

sp_names <- c('PIRU' = 'Picea rubens', 'PIGL' = 'Picea glauca', 'TSCA' = 'Tsuga canadensis', 
              'ABBA' = "Abies balsamea", 'ACRU' = 'Acer rubrum', 'BEPA' = 'Betula papyrifera',
              'PIST' = 'Pinus strobus', 'BEAL' = 'Betula alleghaniensis', 'THOC' = 'Thuja occidentalis',
              'BECO' = 'Betula cordifolia', 'ACPE' = 'Acer pensylvanicum', 'AMELANCHIER' = 'Amelanchier spp')

# Stem map ----------------------------------------------------------------
#Only one site
stem_mapV <- trees %>% filter(Site == 'BM') %>% 
  ggplot(aes(x = cY, y = cX))+
  geom_point(aes(color = Status, size = DBH))+
  #coord_equal()+# makes x and y the same scale
  scale_color_manual(name = "Status", labels = c("Dead", "Live"), 
                       values = c("D" = '#808080', "L" = '#228B22'))+
  scale_size_continuous(range = c(1, 5))+
  facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names))+
  labs(x = "Easting (meters)", y = "Northing (meters)")+ 
  theme(axis.text = element_text(size = 12), 
        strip.text = element_text(size = 16), #facet wrap text size
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "#696969", size = 0.01),
        axis.ticks = element_line(color = "#696969", size = 0.5),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5)),
        aspect.ratio = 10/3)+
  geom_hline(yintercept = 5, linetype = 2)+
  geom_vline(xintercept = 100, linetype = 2)+
  #geom_text(aes(label=Tag), size = 3)+ #add tag # labels
  coord_flip()+
  theme_FVM()
stem_mapV

stem_mapH <- trees %>% filter(Site == 'BM') %>% 
  ggplot(aes(x = cY, y = cX))+
  geom_point(aes(color = Status, size = DBH))+
  #coord_equal()+# makes x and y the same scale
  scale_color_manual(name = "Status", labels = c("Dead", "Live"), 
                     values = c("D" = '#808080', "L" = '#228B22'))+
  scale_size_continuous(range = c(1, 5))+
  facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names))+
  labs(x = "Easting (meters)", y = "Northing (meters)")+ 
  theme(axis.text = element_text(size = 12), 
        strip.text = element_text(size = 16), #facet wrap text size
        axis.title = element_text(size = 16),
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
  geom_vline(xintercept = 100, linetype = 2)+
  #geom_text(aes(label=Tag), size = 3)+ #add tag number labels
  #coord_flip()+
  theme_FVM()
stem_mapH

# Combine 2020/21 trees and saps to match Davis size classes---------------------------
#using Davis size classes, only live trees
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
#Calculating the QMD for each size class and species using the 2020 diameters, should be better than the mid-point of each size class
#
qmd1 <- Ltreesap2 %>% group_by(Species, size_class) %>% summarise(QMD = (sum(DBH^2)),
                                                                  avg_dbh = mean(DBH, na.rm = TRUE),
                                                                  se_dens = sd(DBH, na.rm = TRUE)/
                                                                    sqrt(sum(!is.na(DBH))),
                                                                  num_indiv = sum(!is.na(DBH)))
qmd2 <- qmd1
qmd2$QMD <-(qmd2$QMD = sqrt(qmd2$QMD/qmd2$num_indiv)) # had to split QMD equation into 2 parts

#Adding calculated size class QMD to 1959 tree data
trees59b <- trees59 #combine two smallest size classes to be more consistent with the size of the other size classes
trees59b$size_class <- trees59b$size_class %>% recode(d2.5_4.9 = "d2.5_9.9") %>% 
                                               recode(d5_9.9 = "d2.5_9.9")

trees59qmd <- left_join(trees59b, qmd2, by = c("Species", "size_class"))

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


# Combine with event data
event_tree <- events %>% filter(Module == "Trees")

intersect(names(Comb_tree), names(event_tree))

Comb_tree_event <- left_join(Comb_tree, event_tree, by = c("Site","SampleYear","SampleEventNum"))

sort(unique(Comb_tree_event$size_class))
Comb_tree_event$size_class <- ordered(Comb_tree_event$size_class,
                                levels = c("d2.5_9.9", "d10_19.9", "d20_29.9", 
                                           "d30_39.9", "d40_49.9", "d50_59.9", 
                                           "d60_69.9", "d70_79.9", "d80_89.9", 
                                           "d90_99.9", "d100p"))
levels(Comb_tree_event$size_class) #Size classes are in the order we want them now

Comb_tree_event$Site <- ordered(Comb_tree_event$Site,
                                levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP"))
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
                                    select(Site, SiteName, SampleEventNum, SampleYear, 
                                            size_class, num_stems_ha, BA_m2ha) 

# Plot diameter distribution comparison -----------------------------------
#Connected lines
Comp_tree_dist_plot <- ggplot(data = tree_dist_ha, aes(color = SampleEventNum, x = size_class, y = num_stems_ha))+
  geom_point()+ 
  geom_line(aes(group = SampleEventNum), size = 2)+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 10), # change axis label size
        strip.text = element_text(size = 14), # change facet text size
        axis.title = element_text(size = 16), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020/21'), 
                     values = c("1" = '#6DA346', "2" = '#228B22'))+
  scale_x_discrete(labels= c('2.5 \U2013 10','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
  theme_FVM() 
    
Comp_tree_dist_plot  

#staggered bar plots
Comp_tree_dist_bar <- ggplot(data = tree_dist_ha, aes(x = size_class, y = num_stems_ha, fill = SampleEventNum))+
  geom_bar(width = .75, position = position_dodge(), stat = 'identity')+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 10), # change axis label size
        strip.text = element_text(size = 14), # change facet text size
        axis.title = element_text(size = 16), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))+
  scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020/21'), values = c("1" = '#6DA346', "2" = '#228B22'))+
  scale_x_discrete(labels= c('2.5 \U2013 10','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
  theme_FVM() 

Comp_tree_dist_bar

# Allometric equations ----------------------------------------------------

#Complete tree biomass: Young's complete tree, Aaron Teets code; DBH=diameter at breast height, includes bark, to biomass in kilograms, comp = complete
compmass=function(Species, DBH){
  if(Species=='PIRU'){comp.b0=1.10651;comp.b1=2.298388}
  else if(Species=='ABBA'){comp.b0=0.8161535;comp.b1=2.413978}
  else if(Species=='PIGL'){comp.b0=1.10651;comp.b1=2.298388}#using PIRU values
  else if(Species=='PIST'){comp.b0=0.5724865;comp.b1=2.467798}
  else if(Species=='TSCA'){comp.b0=0.8645225;comp.b1=2.38591}
  else if(Species=='THOC'){comp.b0=1.32942;comp.b1=1.919051}
  else if(Species=='BEAL'){comp.b0=1.345053;comp.b1=2.335473}
  else if(Species=='BEPA'){comp.b0=0.74343;comp.b1=2.639837}
  else if(Species=='BECO'){comp.b0=0.74343;comp.b1=2.639837}#using BEPA values
  else if(Species=='ACRU'){comp.b0=1.187596;comp.b1=2.37025}
  else if(Species=='AMELANCHIER'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='ACPE'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='ACSP'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='ACER'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='SODE'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='PRPE'){comp.b0=1.187596;comp.b1=2.37025} #using ACRU
  else if(Species=='UNKCON'){comp.b0=1.10651;comp.b1=2.298388} #using PIRU
  else {comp.b0=0;comp.b1=0}
  comp=exp(comp.b0+comp.b1*(log(DBH/2.54)))*0.4536
  return(comp)}

#calculating carbon mass from carbon mass using % from Lamlom and Savidge 2003
compcarb=function(Species, biomass){
  if(Species=='PIRU'){per.carb = 0.5039} #using picea glauca 0.5039
  else if(Species=='PIGL'){per.carb = 0.5039}
  else if(Species=='ABBA'){per.carb = 0.5008}
  else if(Species=='PIST'){per.carb = 0.4974}
  else if(Species=='TSCA'){per.carb = 0.5033}
  else if(Species=='THOC'){per.carb = 0.5172}
  else if(Species=='BEAL'){per.carb = 0.4627}
  else if(Species=='BEPA'){per.carb = 0.4837}
  else if(Species=='BECO'){per.carb = 0.4837}#using BEPA
  else if(Species=='ACRU'){per.carb = 0.4864}
  else if(Species=='ACER'){per.carb = 0.4864} #using ACRU
  else if(Species=='AMELANCHIER'){per.carb = 0.4864} #using ACRU
  else if(Species=='ACPE'){per.carb = 0.4864} #using ACRU
  else if(Species=='ACSP'){per.carb = 0.4864} #using ACRU
  else if(Species=='SODE'){per.carb = 0.4864} #using ACRU
  else if(Species=='PRPE'){per.carb = 0.4953} #using Prunus serotina
  else if(Species=='UNKCON'){per.carb = 0.5039} #using Picea glauca
  else {per.carb = 0.0}
  carbon=per.carb*biomass
  return(carbon)}

# Table 1: Calculate density, ba, biomass, carbon mass ---------------------------
#use Comb_tree_event df of trees and saps >2.4 cm created earlier

SumTable <- Comb_tree_event %>% filter(size_class != 'd2.5_9.9') %>% 
                                 mutate(carbon = 0) %>% 
                                 mutate(comp = 0)

SumTable$comp <- mapply(compmass,Species=SumTable$Species, DBH = SumTable$DBH_QMD)
SumTable$carbon <- mapply(compcarb,Species=SumTable$Species, biomass = SumTable$comp)

names(SumTable)
SumTable_site <- SumTable %>% select(-c(NumSubplots, SubplotArea, Module, Species, size_class, DBH_QMD)) %>% 
                                 group_by(Site, SampleEventNum) %>% 
                                 summarise(SampleYear = first(SampleYear),
                                           SiteName = first(SiteName),
                                           sum_stems = sum(num_stem),
                                           sum_BA_cm2 = sum(BA_cm2),
                                           sum_carb = sum(carbon),
                                           sum_comp = sum(comp),
                                           TotArea = first(TotArea)) 

SumTable_ha <- SumTable_site %>%  mutate(num_stems_ha = sum_stems * (10000/TotArea), 
                                         BA_m2ha = sum_BA_cm2/TotArea,
                                         carbonmass_Mgha = (sum_carb * (10000/TotArea))/1000, 
                                         biomass_Mgha = (sum_comp * (10000/TotArea))/1000) %>% 
                                  mutate(across(where(is.numeric), round, 0)) %>% 
                                  select(Site, SiteName, SampleEventNum, 
                                         num_stems_ha, BA_m2ha, biomass_Mgha, carbonmass_Mgha) #converting from plot to ha and kg to megagram 

SumTable_ha$SampleEventNum <- as.character(SumTable_ha$SampleEventNum)

SumTable_ha2 <- SumTable_ha %>% pivot_wider(names_from = SampleEventNum, 
                                           values_from = c(num_stems_ha, BA_m2ha, 
                                                           biomass_Mgha, carbonmass_Mgha))  
                                        
SumTable_ha3 <- left_join(SumTable_ha2, cwd, by = "Site")
names(SumTable_ha3)
SumTable_ha4 <- SumTable_ha3 %>% ungroup() %>% 
                                 mutate(vol_m3ha = round(vol_m3ha)) %>% 
                                 select(-c("carbonmass_Mgha_2", "carbonmass_Mgha_1", "num_pieces","Site"))
                                  
                                      
kable(SumTable_ha4)

#write.csv(Comb_sum_table, './tables/Summary_metrics1959_2020_2.5.csv', row.names = FALSE)

# Basal area by species comparison ----------------------------------------
tree_sp_ha <- Comb_tree_event %>% select(-c(NumSubplots, SubplotArea, Module, size_class, DBH_QMD)) %>% 
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

#removing very infrequent species
sp_rank <- tree_sp_ha %>% group_by(Species, SampleEventNum) %>% 
                        summarize(BA_m2ha = sum(BA_m2ha))
too_few <- c("SODE", "PRPE", "ACSP", "ACER")
top_sev <- c("PIRU", "PIGL", "ABBA", "TSCA", "ACRU", "BEPA", "THOC")
top_four <- c("PIRU", "PIGL", "ABBA", "TSCA")

Comb_sp4 <- tree_sp_ha %>% ungroup() %>% 
                           add_row(Site = 'OP', Species = 'PIGL', BA_m2ha = 0, SampleEventNum = '2') %>% 
                           add_row(Site = 'OP', Species = 'ACRU', BA_m2ha = 0, SampleEventNum = '1') %>%
                           add_row(Site = 'BC', Species = 'ACRU', BA_m2ha = 0, SampleEventNum = '1') %>% 
                           add_row(Site = 'BC', Species = 'THOC', BA_m2ha = 0, SampleEventNum = '1') %>% 
                           add_row(Site = 'BH', Species = 'BEPA', BA_m2ha = 0, SampleEventNum = '1') %>% 
                           add_row(Site = 'BM', Species = 'ACRU', BA_m2ha = 0, SampleEventNum = '1') %>% 
                           add_row(Site = 'BM', Species = 'BEPA', BA_m2ha = 0, SampleEventNum = '2') %>% 
                           add_row(Site = 'IB', Species = 'THOC', BA_m2ha = 0, SampleEventNum = '2') %>% 
                           add_row(Site = 'PM', Species = 'THOC', BA_m2ha = 0, SampleEventNum = '1') %>%
                           add_row(Site = 'WP', Species = 'ACRU', BA_m2ha = 0, SampleEventNum = '1') %>% 
                           filter(Species %in% top_four)# only top four species

unique(Comb_sp4$Species)

Comb_sp4$Site <- ordered(Comb_sp4$Site,
                                levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP"))
levels(Comb_sp4$Site)#Sites are in the ordered by complexity now

Comb_sp_plot <- ggplot(data = Comb_sp4, aes(color = Species, x = SampleEventNum, y = BA_m2ha))+
  geom_point()+ 
  geom_line(aes(group = Species), size = .5)+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  xlab('Sample Year')+
  ylab(bquote('Basal area ('~m^2*'/ha)'))+ 
  theme(axis.text = element_text(size = 10), # change axis label size
        strip.text = element_text(size = 14), # change facet text size
        axis.title = element_text(size = 16), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))+
  scale_x_discrete(labels= c("1" = '1959', "2" = '2020/21'))+ 
  theme_FVM() 

Comb_sp_plot

#Alternate visualizations
#Dumbell plot

#added 0 record for OP PIGL so the change is more obvious
dumbbell <- Comb_sp4 %>% 
  ggplot(aes(x= BA_m2ha, y= reorder(Species, BA_m2ha))) +
  geom_line(aes(group = Species))+
  geom_point(aes(color=SampleEventNum), size=4) +
  theme(legend.position="bottom")+
  facet_wrap(~Site, ncol = 2, labeller = as_labeller(site_names))+
  xlab(bquote('Basal area ('~m^2*'/ha)'))+
  ylab('Species')+ 
  scale_y_discrete(labels = sp_names)+
  scale_color_manual(name = "Year", labels = c("1" = '1959', "2" = '2020/21'), 
                     values = c("1" = '#6DA346', "2" = '#228B22'))+
  theme_FVM()
dumbbell

#staggered bar plot

st_bar <- ggplot(Comb_sp4, aes(fill=SampleEventNum, x=reorder(Species, BA_m2ha), y=BA_m2ha))+
  geom_bar(position = "dodge", stat = "identity")+
  theme(legend.position="bottom")+
  theme(axis.text = element_text(size = 10), # change axis label size
        strip.text = element_text(size = 14), # change facet text size
        axis.title = element_text(size = 16), # change axis title size
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  xlab('Species')+
  ylab(bquote('Basal area ('~m^2*'/ha)'))+ 
  scale_x_discrete(labels = sp_names)+
  scale_fill_manual(name = "Year", labels = c("1" = '1959', "2" = '2020/21'), 
                     values = c("1" = '#6DA346', "2" = '#228B22'))+
  theme_FVM()
st_bar
