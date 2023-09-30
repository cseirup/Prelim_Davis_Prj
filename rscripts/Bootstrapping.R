
#Code to subsample the 10x200m plot as 5 10x10m plots (same as Davis) in order to get a measure of
#variablity across the stand for stem density, basal area, and biomass

#once complete, moved into Ch1_exploratory.R



#devtools::install_github("katemmiller/forestTrends")
library(tidyverse)
library(forestNETN)
library(forestTrends)

source("C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj/rscripts/davis_functions.R")

QC_data_folder = "C:/01_NETN/Forest_Health/R_Dev/Davis_data" #location of QC'd datafiles
files_list <- list.files(path = QC_data_folder, pattern="*.csv", full.names=TRUE) #read names of data files in folder
files_list # review list in order to choose nicknames
nicknames <- c("cwd", "events", "qd_ch", "seeds", "qd_sp", "rwlong", "saps", "soil_d",
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

# Load 2020s tree data ----------------------------------------------------------------
list2env(data["trees"], envir = .GlobalEnv)

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
boot_mod_BC <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'BC') %>% 
                             mutate(boot = .x)) %>%  data.frame()
boot_mod_BH <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'BH') %>% 
                                mutate(boot = .x)) %>%  data.frame()
boot_mod_OP <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'OP') %>% 
                                mutate(boot = .x)) %>%  data.frame()
boot_mod_PM <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'PM') %>% 
                                mutate(boot = .x)) %>%  data.frame()
boot_mod_IB <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'IB') %>% 
                                mutate(boot = .x)) %>%  data.frame()
boot_mod_WP2 <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'WP2') %>% 
                                mutate(boot = .x)) %>%  data.frame()
boot_mod_BM <- purrr::map_dfr(seq_len(num_reps), ~sample_fun_trees(trees4, 'BM') %>% 
                                mutate(boot = .x)) %>%  data.frame()
boot_mod_all <- rbind(boot_mod_BH, boot_mod_BM, boot_mod_WP2, boot_mod_IB, 
                      boot_mod_PM, boot_mod_OP, boot_mod_BC)

# keep data in case want to do something else with the distribution
write.csv(boot_mod_all, './tables/Boot_data_all_sites.csv', row.names = FALSE) 

#calculate 2.5% and 97.5% confidence intervals

BC_CI <- bootCI(boot_mod_BC) #function written to calculate confidence intervals on individual site dataframes

boot_mod_all$Site <- ordered(boot_mod_all$Site,
                      levels = c("BM", "PM", "BC",  "OP", "BH","IB", "WP2"))

Hist_den <- ggplot(data = boot_mod_all, aes(x = Site, y = density_ha))+
              geom_boxplot()

  
Hist_den



