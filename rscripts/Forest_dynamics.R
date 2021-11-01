# Code for calculating various forest dynamics metrics
#katemmiller/forestNETNarch
library(tidyverse)
library(forestMIDNarch)
library(scales)

# Load data ---------------------------------------------------------------

saps <- read.csv("../data/Davis_saps_20200825.csv")
trees <- read.csv("../data/Davis_trees_20200825.csv")
CWD <- readxl::read_xlsx("../data/CWD_data.xlsx")
trees59 <- readxl::read_xlsx("../data/1959_trees_20210409.xlsx")

# 2020 diameter distributions NETN comparion)--------------------------------------------------
  #treesL <- filter(trees, Status == "L") %>% droplevels() #should do live and dead separately later
  #treesD <- filter(trees, Status == "D") %>% droplevels()

trees_df1 <- trees #all trees
trees_df1$BA_cm2 <- round(pi * ((trees_df1$DBH/2)^2), 4)
trees_df2 <- trees_df1 %>% mutate(size_class = as.factor(case_when(between(DBH, 
                                    10, 19.9) ~ "d10_19.9", between(DBH, 20, 29.9) ~ "d20_29.9", 
                                    between(DBH, 30, 39.9) ~ "d30_39.9", between(DBH, 40, 
                                    49.9) ~ "d40_49.9", between(DBH, 50, 59.9) ~ "d50_59.9", 
                                    between(DBH, 60, 69.9) ~ "d60_69.9", between(DBH, 70, 
                                    79.9) ~ "d70_79.9", between(DBH, 80, 89.9) ~ "d80_89.9", 
                                    between(DBH, 90, 99.9) ~ "d90_99.9", DBH >= 100 ~ "d100p", 
                                    TRUE ~ "unknown")), stem = 1)

trees_dist <- trees_df2 %>% group_by(Site, size_class) %>% 
                              summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000) #plot area is 2000 m2

trees_dist_wide_dens <- trees_dist %>% select(Site, size_class, num_stems_ha) %>% 
                                 spread(size_class, num_stems_ha, fill = 0)

trees_dist_wide_ba <- trees_dist %>% select(Site, size_class, BA_m2ha) %>% spread(size_class, BA_m2ha, fill = 0)
  

trees_dist_long <- trees_dist_wide_dens %>% select(Site, d10_19.9:d50_59.9) %>% 
  pivot_longer(cols = c(d10_19.9:d50_59.9), names_to = "size_class",
               values_to = "density")

trees_dist_sum <- trees_dist_long %>% group_by(Site, size_class) %>% #not actually doing anything right now but keeping in for names
  summarize(avg_dens = mean(density, na.rm = TRUE),
            se_dens = sd(density, na.rm = TRUE)/
             sqrt(sum(!is.na(density))))
sort(unique(trees_dist_sum$size_class))

trees_dist_sum$size_class <- ordered(trees_dist_sum$size_class,
                                         levels = c("d10_19.9", "d20_29.9", "d30_39.9",
                                                    "d40_49.9", "d50_59.9", "d60_69.9",
                                                    "d70_79.9", "d80_89.9", "d90_99.9",
                                                    "d100p"))
levels(trees_dist_sum$size_class) #Size classes are in the order we want them now

# Now we need to arrange the data using the new ordered factor levels
trees_dist_sum <- trees_dist_sum %>% arrange(Site, size_class)
head(trees_dist_sum)
# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 'PM' = "Pemetic Mtn")

# Make ggplot graph
trees_dist_plot <- ggplot(data = trees_dist_sum, aes(x = size_class, y = avg_dens))+
  geom_bar(stat = 'identity', fill = '#82B07A')+ #revised 7/30/20
  #geom_errorbar(aes(ymin = avg_dens - se_dens, 
                    #ymax = avg_dens + se_dens, x = size_class),
                #color = "#696969", 
                #width = 0.5,
                #size = 0.3,
                #position = position_dodge(0.9))+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ #revised 7-30-20
  theme(axis.text = element_text(size = 6), # change axis label size
        strip.text = element_text(size = 8), # change facet text size
        axis.title = element_text(size = 9), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+   #added 8/3/20                                       
  scale_x_discrete(labels= c('10 \U2013 20', '20 \U2013 30', '30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70 \U2013 80',
                             '80 \U2013 90','90 \U2013 100','100+'))+ #revised 7-30-20
  theme_FVM()
print(trees_dist_plot)
ggsave("./figures/All_trees_diam_dist.jpg", trees_dist_plot, dpi = 300, 
       width = 5.7, height = 4, units = 'in')

#only live trees

treesL_df1 <- filter(trees, Status == "L") %>% droplevels()  #all trees
treesL_df1$BA_cm2 <- round(pi * ((treesL_df1$DBH/2)^2), 4)
treesL_df2 <- treesL_df1 %>% mutate(size_class = as.factor(case_when(between(DBH, 10, 19.9) ~ "d10_19.9", 
                                                                       between(DBH, 20, 29.9) ~ "d20_29.9", 
                                                                       between(DBH, 30, 39.9) ~ "d30_39.9", 
                                                                       between(DBH, 40, 49.9) ~ "d40_49.9", 
                                                                       between(DBH, 50, 59.9) ~ "d50_59.9", 
                                                                       between(DBH, 60, 69.9) ~ "d60_69.9", 
                                                                       between(DBH, 70, 79.9) ~ "d70_79.9", 
                                                                       between(DBH, 80, 89.9) ~ "d80_89.9", 
                                                                       between(DBH, 90, 99.9) ~ "d90_99.9", 
                                                                       DBH >= 100 ~ "d100p", TRUE ~ "unknown")), stem = 1)

treesL_dist <- treesL_df2 %>% group_by(Site, size_class) %>% 
  summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000)#plot area is 2000 m2

treesL_dist_wide_dens <- treesL_dist %>% select(Site, size_class, num_stems_ha) %>% 
  spread(size_class, num_stems_ha, fill = 0)

treesL_dist_wide_ba <- treesL_dist %>% select(Site, size_class, BA_m2ha) %>% spread(size_class, BA_m2ha, fill = 0)


treesL_dist_long <- treesL_dist_wide_dens %>% select(Site, d10_19.9:d50_59.9) %>% 
  pivot_longer(cols = c(d10_19.9:d50_59.9), names_to = "size_class",
               values_to = "density")

treesL_dist_sum <- treesL_dist_long %>% group_by(Site, size_class) %>% # not actually doing anything right now
                   summarize(avg_dens = mean(density, na.rm = TRUE),
                     se_dens = sd(density, na.rm = TRUE)/
                       sqrt(sum(!is.na(density))))
sort(unique(treesL_dist_sum$size_class))

treesL_dist_sum$size_class <- ordered(treesL_dist_sum$size_class,
                                     levels = c("d10_19.9", "d20_29.9", "d30_39.9",
                                                "d40_49.9", "d50_59.9", "d60_69.9",
                                                "d70_79.9", "d80_89.9", "d90_99.9",
                                                "d100p"))
levels(treesL_dist_sum$size_class) #Size classes are in the order we want them now

# Now we need to arrange the data using the new ordered factor levels
treesL_dist_sum <- treesL_dist_sum %>% arrange(Site, size_class)
head(treesL_dist_sum)
# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 'PM' = "Pemetic Mtn")

# Make ggplot graph
treesL_dist_plot <- ggplot(data = treesL_dist_sum, aes(x = size_class, y = avg_dens))+
  geom_bar(stat = 'identity', fill = '#82B07A')+ #revised 7/30/20
  #geom_errorbar(aes(ymin = avg_dens - se_dens, 
  #ymax = avg_dens + se_dens, x = size_class),
  #color = "#696969", 
  #width = 0.5,
  #size = 0.3,
  #position = position_dodge(0.9))+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ #revised 7-30-20
  theme(axis.text = element_text(size = 6), # change axis label size
        strip.text = element_text(size = 8), # change facet text size
        axis.title = element_text(size = 9), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  scale_x_discrete(labels= c('10 \U2013 20', '20 \U2013 30', '30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70 \U2013 80',
                             '80 \U2013 90','90 \U2013 100','100+'))+ #revised 7-30-20
  ggtitle("Only Live Trees")+
  theme_FVM()     
print(treesL_dist_plot)
ggsave("./figures/Live_treesL_diam_dist.jpg", treesL_dist_plot, dpi = 300, 
       width = 5.7, height = 4, units = 'in')


# 2020 Diameter Distribution (Davis comparison) --------------------------------------
#using Davis size classes, only live trees
saps2 <- saps %>% drop_na() %>% add_column(Tag = NA)
treesL_df1 <- filter(trees, Status == "L") %>% droplevels()  #all trees
names(treesL_df1)
treesL_df1.5 <- treesL_df1 %>% add_column(Quadrat = NA) %>% select(Site, Transect, Quadrat, Tag, Species, DBH)
saps3 <- saps2 %>% select(Site, Transect, Quadrat, Tag, Species, DBH) 
Ltreesap <- rbind(treesL_df1.5, saps3)
Ltreesap$BA_cm2 <- round(pi * ((Ltreesap$DBH/2)^2), 4)

Ltreesap2 <- Ltreesap %>% mutate(size_class = as.factor(case_when(between(DBH, 1.0, 2.4) ~ "d1_2.4", 
                                                                     between(DBH, 2.5, 4.9) ~ "d2.5_4.9", 
                                                                     between(DBH, 5, 9.9) ~ "d5_9.9", 
                                                                     between(DBH, 10, 19.9) ~ "d10_19.9", 
                                                                     between(DBH, 20, 29.9) ~ "d20_29.9", 
                                                                     between(DBH, 30, 39.9) ~ "d30_39.9", 
                                                                     between(DBH, 40, 49.9) ~ "d40_49.9", 
                                                                     between(DBH, 50, 59.9) ~ "d50_59.9", 
                                                                     between(DBH, 60, 69.9) ~ "d60_69.9", 
                                                                     DBH >= 70 ~ "d70", TRUE ~ "unknown")), stem = 1)

treesapL_dist <- Ltreesap2 %>% group_by(Site, size_class) %>% 
  summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000)#plot area is 2000 m2; not currently visualizing the ba

treesapL_dist_wide_dens <- treesapL_dist %>% select(Site, size_class, num_stems_ha) %>% #
  spread(size_class, num_stems_ha, fill = 0)

treessapL_dist_wide_ba <- treesapL_dist %>% select(Site, size_class, BA_m2ha) %>% spread(size_class, BA_m2ha, fill = 0)

treesapL_dist_long <- treesapL_dist_wide_dens %>% select(Site, d1_2.4:d50_59.9) %>% 
  pivot_longer(cols = c(d1_2.4:d50_59.9), names_to = "size_class",
               values_to = "density")

treesapL_dist_sum <- treesapL_dist_long %>% group_by(Site, size_class) %>% # not actually doing anything right now
  summarize(avg_dens = mean(density, na.rm = TRUE),
            se_dens = sd(density, na.rm = TRUE)/
              sqrt(sum(!is.na(density))))

sort(unique(treesapL_dist_sum$size_class))

treesapL_dist_sum$size_class <- ordered(treesapL_dist_sum$size_class,
                                      levels = c("d1_2.4", "d2.5_4.9", "d5_9.9", "d10_19.9", "d20_29.9", "d30_39.9",
                                                 "d40_49.9", "d50_59.9", "d60_69.9",
                                                 "d70_79.9", "d80_89.9", "d90_99.9",
                                                 "d100p"))
levels(treesapL_dist_sum$size_class) #Size classes are in the order we want them now

# Now we need to arrange the data using the new ordered factor levels
treesapL_dist_sum <- treesapL_dist_sum %>% arrange(Site, size_class)
head(treesapL_dist_sum)
# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 'PM' = "Pemetic Mtn")

# Make ggplot graph
treesapL_dist_plot <- ggplot(data = treesapL_dist_sum, aes(x = size_class, y = avg_dens))+
  geom_bar(stat = 'identity', fill = '#82B07A')+ #revised 7/30/20
  #geom_errorbar(aes(ymin = avg_dens - se_dens, 
  #ymax = avg_dens + se_dens, x = size_class),
  #color = "#696969", 
  #width = 0.5,
  #size = 0.3,
  #position = position_dodge(0.9))+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ #revised 7-30-20
  theme(axis.text = element_text(size = 6), # change axis label size
        strip.text = element_text(size = 8), # change facet text size
        axis.title = element_text(size = 9), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  scale_x_discrete(labels= c('1 \U2013 2.4', '2.5 \U2013 4.9', '5 \U2013 9.9','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ #revised 7-30-20
  ggtitle("Only Live Trees - Davis size classes")+
  theme_FVM()     
print(treesapL_dist_plot)
ggsave("./figures/Live_treesap_davis_size_classes.jpg", treesapL_dist_plot, dpi = 300, 
       width = 5.7, height = 4, units = 'in')


# 1959 basal area/density -------------------------------------------------
trees59qmd2$BA_cm2 <- round(pi * (((trees59qmd2$QMD/2)^2)*trees59qmd2$num_stem), 4) #basal area per site
tree59_dist_allo <- trees59qmd2 %>% group_by(Site, size_class, Species) %>% 
  summarise(num_stems_ha = sum(num_stem) * 20, BA_m2ha = sum(BA_cm2)/500)#Davis plots with 5 10x10m, for allometric equations

tree59_dist <- trees59qmd2 %>% group_by(Site, size_class) %>% 
  summarise(num_stems_ha = sum(num_stem) * 20, BA_m2ha = sum(BA_cm2)/500)# for diameter distribution, not split by species
# 1959 diameter distribution ----------------------------------------------
tree59_distL_dist_wide_dens <- tree59_dist %>% select(Site, size_class, num_stems_ha) %>% 
  spread(size_class, num_stems_ha, fill = 0)

tree59_dist_wide_ba <- tree59_dist %>% select(Site, size_class, BA_m2ha) %>% spread(size_class, BA_m2ha, fill = 0)# not visualizing this right now

tree59_dist_long <- tree59_distL_dist_wide_dens %>% select(Site, d10_19.9:d50_59.9) %>% 
  pivot_longer(cols = c(d10_19.9:d50_59.9), names_to = "size_class",
               values_to = "density")

tree59_dist_sum <- tree59_dist_long %>% group_by(Site, size_class) %>% # not actually doing anything right now
  summarize(avg_dens = mean(density, na.rm = TRUE),
            se_dens = sd(density, na.rm = TRUE)/
              sqrt(sum(!is.na(density))))

sort(unique(tree59_dist_sum$size_class))

#Attempt to combine 2020 and 1959 data together for visulization
tree59_dist_sum2 <- tree59_dist_sum %>% mutate(year = 1959)
treesapL_dist_sum2 <- treesapL_dist_sum %>% filter(size_class != "d1_2.4") %>%  mutate(year = 2020)
Comp_tree_dist <- rbind(tree59_dist_sum2, treesapL_dist_sum2)
Comp_tree_dist$year <- as.character(Comp_tree_dist$year)
sort(unique(Comp_tree_dist$size_class))

Comp_tree_dist$size_class <- ordered(Comp_tree_dist$size_class,
                                      levels = c("d1_2.4", "d2.5_4.9", "d5_9.9", "d10_19.9", "d20_29.9", "d30_39.9",
                                      "d40_49.9", "d50_59.9", "d60_69.9",
                                     "d70_79.9", "d80_89.9", "d90_99.9",
                                      "d100p"))
levels(Comp_tree_dist$size_class) #Size classes are in the order we want them now


# Now we need to arrange the data using the new ordered factor levels
Comp_tree_dist <- Comp_tree_dist %>% arrange(Site, size_class)
head(Comp_tree_dist)
# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 'PM' = "Pemetic Mtn", 'WM' = "Western Mtn", 'IB' = 'Ironbound', 'BH' = "Bass Harbor Head")

# Make ggplot graph
Comp_tree_dist_plot <- ggplot(data = Comp_tree_dist, aes(fill = year, x = size_class, y = avg_dens))+
  geom_bar(width = .75, position = position_dodge(.75), stat = 'identity')+ #revised 7/30/20
  facet_wrap(~Site, ncol = 7, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ #revised 7-30-20
  theme(axis.text = element_text(size = 8), # change axis label size
        strip.text = element_text(size = 10), # change facet text size
        axis.title = element_text(size = 11), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  scale_x_discrete(labels= c('1 \U2013 2.4', '2.5 \U2013 4.9', '5 \U2013 9.9','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ #revised 7-30-20
  ggtitle("1959 vs 2020 - live trees only")+
  theme_FVM()     
print(Comp_tree_dist_plot)
ggsave("./figures/Comp_tree_dist_1959_2020.jpg", Comp_tree_dist_plot, dpi = 300, 
       width = 5.7, height = 4, units = 'in')

#something going wrong with this figure
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
  else if(Species=='ACRU'){comp.b0=1.187596;comp.b1=2.37025}
  else if(Species=='AMELANCHIER'){comp.b0=1.187596;comp.b1=2.37025} #same as ACRU
  else if(Species=='ACPE'){comp.b0=1.187596;comp.b1=2.37025} #same as ACRU
  else if(Species=='PRPE2'){comp.b0=1.187596;comp.b1=2.37025} #same as ACRU
  else if(Species=='UNKCON'){comp.b0=1.10651;comp.b1=2.298388} #same as PIRU
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
  else if(Species=='ACRU'){per.carb = 0.4864}
  else if(Species=='AMELANCHIER'){per.carb = 0.4864} #using ACRU
  else if(Species=='ACPE'){per.carb = 0.4864} #using ACRU
  else if(Species=='PRPE2'){per.carb = 0.4953} #using Prunus serotina
  else if(Species=='UNKCON'){per.carb = 0.5039} #using Picea glauca
  else {per.carb = 0.0}
  carbon=per.carb*biomass
  return(carbon)}

# 2020 trees: density, ba, biomass, carbon mass ---------------------------
#combining 2020 trees and sapling (>2.4cm) 
#Trying without saps
saps2 <- saps %>% drop_na() %>% mutate(stem=1) %>% select(Site, Transect, Species, DBH, stem) 
treesL <- filter(trees, Status == "L") %>% droplevels() %>% mutate(stem=1) %>% select(Site, Transect, Species, DBH, stem)
names(saps2);names(treesL)
Ltreesap <- rbind(treesL, saps2)
Ltreesap2 <- Ltreesap %>% filter(DBH > 2.4) # Davis included trees/saps over 1in DBH
#Ltreesap2$BA_cm2 <- round(pi * ((Ltreesap2$DBH/2)^2), 4)# add column for basal area
treesL$BA_cm2 <- round(pi * ((treesL$DBH/2)^2), 4)# add column for basal area

#Calculate biomass and carbon mass for 2020 data
#Mass2020 <- Ltreesap2 # live trees, including saps over 2.5 cm
Mass2020 <- treesL # live trees >10cm
Mass2020$carbon <- 0
Mass2020$comp <- mapply(compmass,Species=Mass2020$Species, DBH = Mass2020$DBH)
Mass2020$carbon <- mapply(compcarb,Species=Mass2020$Species, biomass = Mass2020$comp)
head(Mass2020)

Mass2020_ha <- Mass2020 %>% group_by(Site, Species) %>% 
                  summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000, 
                            biomass_mgha = (sum(comp)*5)/1000, carbonmass_mgha = (sum(carbon)*5)/1000) 
#2020 plot area is 2000 m2; coversion facto to hectare is 5; converting from kg to megagram

Table2020 <- Mass2020_ha %>% pivot_wider(names_from = Site, values_from = c( num_stems_ha, BA_m2ha, biomass_mgha, carbonmass_mgha), values_fill = 0) %>% arrange(desc(num_stems_ha_BC))# dispay as table
                                                            
write.csv(Table2020, './tables/Summary_metrics2020_greater10cm.csv', row.names = FALSE)


# Quadratic Mean Diameter -------------------------------------------------
#Calculating the QMD for each Davis size class using the 2020 diameters, should be better than the mid-point of each size class
Ltreesap3 <- Ltreesap2 %>% mutate(size_class = as.factor(case_when(between(DBH, 1.0, 2.4) ~ "d1_2.4", 
                                                                  between(DBH, 2.5, 4.9) ~ "d2.5_4.9", 
                                                                  between(DBH, 5, 9.9) ~ "d5_9.9", 
                                                                  between(DBH, 10, 19.9) ~ "d10_19.9", 
                                                                  between(DBH, 20, 29.9) ~ "d20_29.9", 
                                                                  between(DBH, 30, 39.9) ~ "d30_39.9", 
                                                                  between(DBH, 40, 49.9) ~ "d40_49.9", 
                                                                  between(DBH, 50, 59.9) ~ "d50_59.9", 
                                                                  between(DBH, 60, 69.9) ~ "d60_69.9", 
                                                                  DBH >= 70 ~ "d70", TRUE ~ "unknown")), stem = 1)


qmd1 <- Ltreesap3 %>% group_by(Species, size_class) %>% summarise(QMD = (sum(DBH^2)),
                                                                  avg_dbh = mean(DBH, na.rm = TRUE),
                                                                  se_dens = sd(DBH, na.rm = TRUE)/
                                                                    sqrt(sum(!is.na(DBH))),
                                                                  num_indiv = sum(!is.na(DBH)))
qmd2 <- qmd1
qmd2$QMD <-(qmd2$QMD = sqrt(qmd2$QMD/qmd2$num_indiv)) # had to split QMD equation into 2 parts

#Adding calculated size class QMD to 1959 tree data
trees59qmd <- left_join(trees59, qmd2, by = c("Species", "size_class"))

midpoints <- data.frame("size_class" = c("d1_2.4", "d2.5_4.9", "d5_9.9", "d10_19.9", "d20_29.9", "d30_39.9", "d40_49.9", "d50_59.9", "d60_69.9", "d70"), midpoint = c(1.7, 3.7, 7.45, 15, 25, 35, 45, 55, 65, 70))

trees59qmd <- left_join(trees59qmd, midpoints, by = "size_class")

trees59qmd2 <- trees59qmd %>% select("Site", "Species", "size_class", "num_stem", "QMD", "midpoint")

trees59qmd2$QMD <- coalesce(trees59qmd2$QMD, trees59qmd2$midpoint) #replacing missing QMDs with midpoints for that size class/species

# 1959 trees: density, ba, biomass, carbon mass ---------------------------
Mass1959 <- trees59qmd2 # live trees, including saps over 2.5 cm
Mass1959$BA_cm2 <- round(pi * ((Mass1959$QMD/2)^2), 4)*Mass1959$num_stem# add column for basal area
Mass1959$comp <- 0
Mass1959$comp <- mapply(compmass, Species=Mass1959$Species, DBH = Mass1959$QMD)
Mass1959$comp2 <- Mass1959$comp*Mass1959$num_stem
Mass1959$carbon <- 0
Mass1959$carbon <- mapply(compcarb,Species=Mass1959$Species, biomass = Mass1959$comp2)
head(Mass1959)
sort(unique(Mass1959$size_class))
g10 <- c( "d10_19.9", "d20_29.9", "d30_39.9", "d40_49.9", "d50_59.9")
Mass1959 <-  filter(Mass1959, size_class %in% g10)# selecting only stems >10cm
Mass1959_ha <- Mass1959 %>% group_by(Site, Species) %>% 
                            summarise(num_stems_ha = sum(num_stem) * 20, BA_m2ha = sum(BA_cm2)/500, 
                                      biomass_mgha = (sum(comp2)*20/1000), carbonmass_mgha = (sum(carbon)*20)/1000)
                           
#1959 plot area is 500 m2; conversion factor to hectare is 20; converting to megagrams

Table1959 <- Mass1959_ha %>% pivot_wider(names_from = Site, values_from = c( num_stems_ha, BA_m2ha, biomass_mgha, carbonmass_mgha), values_fill = 0) %>% arrange(desc(num_stems_ha_BC))# display as table

write.csv(Table1959, './tables/Summary_metrics1959b_10cm.csv', row.names = FALSE)

# Graphing summary metrics for comparison -----------------------------------------
step1 <- Mass1959_ha %>% mutate(year = 1959) 
step2 <- Mass2020_ha %>%  mutate(year = 2020)
Combo_sum <- rbind(step1, step2)
Combo_sum$year <- as.character(Combo_sum$year)
sampled_sites <- c("BC", "BM", "OP", "PM")
Combo_sum2 <- Combo_sum %>% filter(Site %in% sampled_sites)
table(Combo_sum2$Site)

#Sort species in a useful order
sort(unique(Combo_sum2$Species))
Combo_sum2$Species <- ordered(Combo_sum2$Species, levels = c("PIRU", "ABBA", "TSCA", "BEPA", "THOC", "ACRU", "ACPE", "PIGL", "AMELANCHIER", "PIST"))
levels(Combo_sum2$Species)

# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 'PM' = "Pemetic Mtn", 'WM' = "Western Mtn", 'IB' = 'Ironbound', 'BH' = "Bass Harbor Head")

names(Combo_sum2)

# Make ggplot graphs
Combo_sum_plot <- ggplot(data = Combo_sum2, aes(x = Species, y = num_stems_ha, fill = year))+
  geom_bar(width = .75, position = position_dodge(.75), stat = 'identity')+ 
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Species", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 6), 
        strip.text = element_text(size = 8), 
        axis.title = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
        ggtitle("1959 vs 2020 - Tree density (live trees >10cm)")+
        theme_FVM()     
print(Combo_sum_plot)
ggsave("./figures/Combo_sum_plot_density2.jpg", Combo_sum_plot, dpi = 300, 
       width = 7, height = 4, units = 'in')

Combo_sum_plot <- ggplot(data = Combo_sum2, aes(x = Species, y = BA_m2ha, fill = year))+
  geom_bar(width = .75, position = position_dodge(.75), stat = 'identity')+ 
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Species", y = "BA_m2ha")+ 
  theme(axis.text = element_text(size = 4), 
        strip.text = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  ggtitle("1959 vs 2020 - Basal Area per ha (live trees >10cm)")+
  theme_FVM()     
print(Combo_sum_plot)
ggsave("./figures/Combo_sum_plot_BA_m2ha2.jpg", Combo_sum_plot, dpi = 300, 
       width = 7, height = 4, units = 'in')

Combo_sum_plot <- ggplot(data = Combo_sum2, aes(x = Species, y = biomass_mgha, fill = year))+
  geom_bar(width = .75, position = position_dodge(.75), stat = 'identity')+ 
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Species", y = "biomass_mgha")+ 
  scale_y_continuous(labels = comma)+
  theme(axis.text = element_text(size = 6), 
        strip.text = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10)))+
  ggtitle("1959 vs 2020 - Biomass per ha (live trees >10cm)")+
  theme_FVM()     
print(Combo_sum_plot)
ggsave("./figures/Combo_sum_plot_biomass_mgha2.jpg", Combo_sum_plot, dpi = 300, 
       width = 7, height = 4, units = 'in')

Combo_sum_plot <- ggplot(data = Combo_sum2, aes(x = Species, y = carbonmass_mgha, fill = year))+
  geom_bar(width = .75, position = position_dodge(.75), stat = 'identity')+ 
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Species", y = "carbonmass_mgha")+ 
  scale_y_continuous(labels = comma)+
  theme(axis.text = element_text(size = 6), 
        strip.text = element_text(size = 6), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  ggtitle("1959 vs 2020 - Carbon mass per ha (live trees >10cm)")+
  theme_FVM()     
print(Combo_sum_plot)
ggsave("./figures/Combo_sum_plot_carbonmass_mgha2.jpg", Combo_sum_plot, dpi = 300, 
       width = 7, height = 4, units = 'in')
