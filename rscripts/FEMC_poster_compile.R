#Source code for figures and tables for FEMC poster 12/2021

#install.packages(c("scales"))
#devtools::install_github("katemmiller/forestMIDN")

library(tidyverse)
library(forestMIDN)
library(knitr)
library(kableExtra)


# Load data ---------------------------------------------------------------

QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("qd_ch", "seeds", "qd_sp", "saps", "soil_d", "trees59", "trees", "mapped_ibuttons")#preferred names for each df
data <- files_list %>% map(read.csv) %>% set_names(nm = nicknames) #load datafiles into a list and rename with nicknames
list2env(data, envir = .GlobalEnv)

# 2020/21 basal area and density for Davis comparison ---------------------------
#using Davis size classes, only live trees
saps2 <- saps %>% drop_na() %>% add_column(Tag = NA) %>% filter(DBH>=2.5)#only taking >2.5cm to match Davis size classes
treesL_df1 <- filter(trees, Status == "L") %>% droplevels()  #only live trees
names(treesL_df1)
treesL_df1.5 <- treesL_df1 %>% select(Site, Transect, Tag, Species, DBH)
saps3 <- saps2 %>% select(Site, Transect, Tag, Species, DBH) 
Ltreesap <- rbind(treesL_df1.5, saps3)
Ltreesap$BA_cm2 <- round(pi * ((Ltreesap$DBH/2)^2), 4)

Ltreesap2 <- Ltreesap %>% mutate(size_class = as.factor(case_when(between(DBH, 2.5, 9.9) ~ "d2.5_9.9", 
                                                                  between(DBH, 10, 19.9) ~ "d10_19.9", 
                                                                  between(DBH, 20, 29.9) ~ "d20_29.9", 
                                                                  between(DBH, 30, 39.9) ~ "d30_39.9", 
                                                                  between(DBH, 40, 49.9) ~ "d40_49.9", 
                                                                  between(DBH, 50, 59.9) ~ "d50_59.9", 
                                                                  between(DBH, 60, 69.9) ~ "d60_69.9", 
                                                                  DBH >= 70 ~ "d70", TRUE ~ "unknown")), stem = 1)

Ltreesap3 <- Ltreesap2 %>% group_by(Site, size_class) %>% 
  summarise(num_stems_ha = sum(stem) * 5, BA_m2ha = sum(BA_cm2)/2000)#plot area is 2000 m2

Ltreesap_den <- Ltreesap3 %>% select(Site, size_class, num_stems_ha) %>% # frequency data for diameter dist
                                             arrange(Site, size_class) %>% 
                                             pivot_wider(names_from = size_class, 
                                                         values_from = num_stems_ha, 
                                                         values_fill = NA)

Ltreesap_ba <- Ltreesap3 %>% select(Site, size_class, BA_m2ha) %>%  #basal area data for table + mass calc
                                 arrange(Site, size_class) %>% 
                                 pivot_wider(names_from = size_class, 
                                 values_from = BA_m2ha, 
                                 values_fill = 0) 



# 2020/2021 diameter distribution for Davis comparison --------------------
Ltreesap_den_long <- Ltreesap_den %>% select(everything()) %>% 
                                      pivot_longer(cols = c(-Site), 
                                                   names_to = "size_class",
                                                   values_to = "density")

sort(unique(Ltreesap_den_long$size_class))

Ltreesap_den_long$size_class <- ordered(Ltreesap_den_long$size_class,
                                        levels = c("d2.5_9.9", 
                                                   "d10_19.9", "d20_29.9", "d30_39.9",
                                                   "d40_49.9", "d50_59.9", "d60_69.9",
                                                   "d70_79.9", "d80_89.9", "d90_99.9",
                                                   "d100p"))
levels(Ltreesap_den_long$size_class) #Size classes are in the order we want them now

# Now we need to arrange the data using the new ordered factor levels
Ltreesap_den_long <- Ltreesap_den_long %>% arrange(Site, size_class)
head(Ltreesap_den_long)
# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 
                'PM' = "Pemetic Mtn", 'BH' = 'Bass Harbor Head', 
                'IB' = 'Ironbound Island', 'WP' = 'Western Mtn')


# Make ggplot graph
Ltreesap_DBH_plot <- ggplot(data = Ltreesap_den_long, aes(x = size_class, y = density))+
  geom_bar(stat = 'identity', fill = '#4DB6D0')+ 
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ #
  theme(axis.text = element_text(size = 6), # change axis label size
        strip.text = element_text(size = 8), # change facet text size
        axis.title = element_text(size = 9), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  scale_x_discrete(labels= c('2.5 \U2013 10','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
   theme_FVM() 
   #ggtitle("Only Live Trees - Davis size classes")+
     
print(Ltreesap_DBH_plot)

#ggsave("./figures/Live_treesap_davis_size_classes.jpg", Ltreesap_DBH_plot, dpi = 300, 
      # width = 5.7, height = 4, units = 'in')


# Calculating Quadratic Mean Diameter -------------------------------------------------
#Calculating the QMD for each Davis size class using the 2020 diameters, should be better than the mid-point of each size class
Ltreesap3.5 <- Ltreesap2 %>% mutate(size_class = as.factor(case_when(between(DBH, 2.5, 9.9) ~ "d2.5_9.9",
                                                                   between(DBH, 10, 19.9) ~ "d10_19.9", 
                                                                   between(DBH, 20, 29.9) ~ "d20_29.9", 
                                                                   between(DBH, 30, 39.9) ~ "d30_39.9", 
                                                                   between(DBH, 40, 49.9) ~ "d40_49.9", 
                                                                   between(DBH, 50, 59.9) ~ "d50_59.9", 
                                                                   between(DBH, 60, 69.9) ~ "d60_69.9", 
                                                                   DBH >= 70 ~ "d70", TRUE ~ "unknown")), stem = 1)


qmd1 <- Ltreesap3.5 %>% group_by(Species, size_class) %>% summarise(QMD = (sum(DBH^2)),
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

trees59qmd2 <- trees59qmd %>% select("Site", "Species", "size_class", "num_stem", "QMD", "midpoint")

trees59qmd2$QMD <- coalesce(trees59qmd2$QMD, trees59qmd2$midpoint) #replacing missing QMDs with midpoints for that size class/species

# 1959 basal area and density -------------------------------------------------
trees59qmd2$BA_cm2 <- round(pi * (((trees59qmd2$QMD/2)^2)*trees59qmd2$num_stem), 4) #basal area per site
tree59_dist_allo <- trees59qmd2 %>% group_by(Site, size_class, Species) %>% 
  summarise(num_stems_ha = sum(num_stem) * 20, BA_m2ha = sum(BA_cm2)/500)#Davis plots with 5 10x10m, split by species for allometric equations

tree59_dist <- trees59qmd2 %>% group_by(Site, size_class) %>% # not by species for diameter dist
  summarise(num_stems_ha = sum(num_stem) * 20, BA_m2ha = sum(BA_cm2)/500)

tree59_den <- tree59_dist %>% select(Site, size_class, num_stems_ha) %>% # frequency data for diameter dist
                              arrange(Site, size_class) %>% 
                              pivot_wider(names_from = size_class, 
                              values_from = num_stems_ha, 
                              values_fill = NA)

tree59_ba <- tree59_dist %>% select(Site, size_class, BA_m2ha) %>%  #basal area data for table + mass calc
                             arrange(Site, size_class) %>% 
                             pivot_wider(names_from = size_class, 
                             values_from = BA_m2ha, 
                             values_fill = 0) 

# 1959 diameter distribution ----------------------------------------------
tree59_den_long <- tree59_den %>% select(everything()) %>% 
                                  pivot_longer(c(-Site), 
                                  names_to = "size_class",
                                  values_to = "density")

sort(unique(tree59_den_long$size_class))

tree59_den_long$size_class <- ordered(tree59_den_long$size_class,
                                        levels = c("d2.5_9.9", 
                                                   "d10_19.9", "d20_29.9", "d30_39.9",
                                                   "d40_49.9", "d50_59.9", "d60_69.9",
                                                   "d70_79.9", "d80_89.9", "d90_99.9",
                                                   "d100p"))
levels(tree59_den_long$size_class) #Size classes are in the order we want them now

# Now we need to arrange the data using the new ordered factor levels
tree59_den_long <- tree59_den_long %>% arrange(Site, size_class)
head(tree59_den_long)
# Set up labels for facet wrap
site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 
                'PM' = "Pemetic Mtn", 'BH' = 'Bass Harbor Head', 
                'IB' = 'Ironbound Island', 'WP' = 'Western Mtn')

# Make ggplot graph
tree59_DBH_plot <- ggplot(data = tree59_den_long, aes(x = size_class, y = density))+
  geom_bar(stat = 'identity', fill = '#D9717D')+ 
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ #
  theme(axis.text = element_text(size = 6), # change axis label size
        strip.text = element_text(size = 8), # change facet text size
        axis.title = element_text(size = 9), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  scale_x_discrete(labels= c('2.5 \U2013 10','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
   theme_FVM()
#ggtitle("1959: Only Live Trees - Davis size classes")+
     
print(tree59_DBH_plot)

#ggsave("./figures/1959_Live_treesap_davis_size_classes.jpg", tree59_DBH_plot, dpi = 300, 
# width = 5.7, height = 4, units = 'in')


# Diameter distributions with 1959 + 2020/21 ------------------------------

comp1 <- tree59_den_long %>% add_column(year = 1959)
comp1$year <- as.character(comp1$year)

comp2 <- Ltreesap_den_long %>% mutate(year = '2020/21') 
comp2$year <- as.character(comp2$year)

Comp_tree_dist <- rbind(comp1, comp2)
sort(unique(Comp_tree_dist$size_class))

#Comp_tree_dist <- Comp_tree_dist %>% arrange(Site, size_class)
head(Comp_tree_dist)

# Make ggplot graph
Comp_tree_dist_plot <- ggplot(data = Comp_tree_dist, aes(color = year, x = size_class, y = density))+
  geom_point()+ 
  geom_line(aes(group = year))+
  facet_wrap(~Site, ncol = 4, labeller = as_labeller(site_names))+
  labs(x = "Tree Diameter Class (cm)", y = "stems/ha")+ 
  theme(axis.text = element_text(size = 6), # change axis label size
        strip.text = element_text(size = 8), # change facet text size
        axis.title = element_text(size = 9), # change axis title size
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 5)))+
  scale_x_discrete(labels= c('2.5 \U2013 10','10 \U2013 20',
                             '20 \U2013 30','30 \U2013 40','40 \U2013 50',
                             '50 \U2013 60','60 \U2013 70','70+'))+ 
  theme_FVM() 
#ggtitle("1959 vs 2020/21 - live trees only")+
      
print(Comp_tree_dist_plot)
#ggsave("./figures/Comp_tree_dist_1959_2020.jpg", Comp_tree_dist_plot, dpi = 300, 
      # width = 5.7, height = 4, units = 'in')
