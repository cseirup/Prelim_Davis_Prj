#Quality control (QC) of Davis site data collecting 2020/2021

library(tidyverse)
library(forestMIDN)
library(readxl)
library(lubridate)

setwd('C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj')

# Load data ---------------------------------------------------------------
saps <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_saps_20221206.csv")
trees <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_trees_20221206.csv")
seeds <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_seedlings_long_20221209.csv")
qd_sp <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_species_long_20221209.csv")
qd_ch <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_char_long_20221209.csv")
soil <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_soil_depth_data_20221011.xlsx")
#need 2022 ibutton data
ibutton_d <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/ibutton_dat_daily_2021-06-29.csv")
ibutton_m <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/ibutton_dat_monthly_2021-06-29.csv")
trees_59<- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/1959_trees_20210409.xlsx")
#still to check: CWD, ibutton, digitized davis data (trees done)

# QC checks for tree data -------------------------------------------------
str(trees)
#Site
names(trees)
unique(trees$Site) # 8 site codes

#Transect
unique(trees$Transect) # 4 transects

#Tag
unique(trees$Tag)# One stump <1.3m recorded -> delete
stump <- trees %>% filter(trees$Tag == "Stump")
trees1 <- anti_join(trees, stump) #removed stump from PM trees
ibuttons <- trees1 %>% filter(is.na(trees1$Tag)) # 10 ibutton records
trees2 <- anti_join(trees1, ibuttons)
unique(trees2$Tag)

write.csv(ibuttons, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Mapped_ibuttons.csv", row.names = FALSE)

#Status
unique(trees2$Status)
table(trees2$Status) # 2 live trees recorded as l instead of L
trees2$Status <- recode(trees2$Status, l = "L")

#Species
unique(trees2$Species)
table(trees2$Species) #fix spelling errors
trees2$Species <- recode(trees2$Species, SWSPP = "UNKCON")
trees2$Species <- recode(trees2$Species, AMEL = "AMELANCHIER")
trees2$Species <- recode(trees2$Species, AMALANCHIER = "AMELANCHIER")
trees2$Species <- recode(trees2$Species, ACAC = "ACRU") #Checked datasheet, ACRU makes the most sense

#DBH
unique(trees2$DBH)
not_tree <- trees2 %>% filter(DBH < 10) # none < 10 cm DBH
DBH_plot <- ggplot(trees2, aes(x = DBH))+
              geom_histogram()+
              facet_wrap(~Site)+
              theme_FVM()
DBH_plot

DBH_plot2 <- ggplot(trees2, aes(x = Site, y = DBH))+
  geom_boxplot()
DBH_plot2

# X and Y
unique(trees2$X)#typo 25 should be 2.5, fixed on datasheet
unique(trees2$Y)

XY_plot <- ggplot(trees2, aes(x=X, y=Y))+
             geom_point()+
             theme_FVM()
XY_plot #looks good, all btwn 0-100 and 0-5

#Crown
trees2L <- filter(trees2, Status == "L")
trees2D <- filter(trees2, Status == "D")
table(trees2L$Crown)
sum(is.na(trees2L$Crown))#one NA crown class for a live tree
trees2L_na <- filter(trees2L, is.na(Crown)) # Tree 47 at IB is missing crown class (checked datasheet). Leaving in for now

unique(trees2D$Crown)# no crown classes for dead trees - good

treesWP2 <- trees2 %>% filter(Site == "WP2")
#cX and cY
unique(treesWP2$cX)
unique(treesWP2$cY)

cXcY_plot2 <- ggplot(treesWP2, aes(x=cX, y=cY))+
  geom_point()+
  theme_FVM()
cXcY_plot2 #looks good

#cX and cY
unique(trees2$cX)
unique(trees2$cY)

cXcY_plot <- ggplot(trees2, aes(x=cX, y=cY))+
  geom_point()+
  theme_FVM()
cXcY_plot #looks good

#adding year, sample event, number of subplots, and subplot area to data
s2020 <- c("OP", "BM", "BC", "PM") # 4 sites sampled in 2020
s2021 <- c("IB", "WP", "BH") # 3 sites sampled in 2021
s2022 <- c("WP2") # 1 site sampled in 2022

trees3 <- trees2 %>% mutate(Date = case_when(Site == 'WP2' ~ '2022-10-01',
                                             Site == 'IB' ~ '2021-07-07',
                                             Site == 'WP' ~ '2021-07-21',
                                             Site == 'OP' ~ '2020-08-06',
                                             Site == 'PM' ~ '2020-08-20',
                                             Site == 'BH' ~ '2021-06-18',
                                             Site == 'BC' ~ '2020-08-17',
                                             Site == 'BM' ~ '2020-08-13')) %>% 
                     mutate(SampleEventNum = 2) %>% 
                     mutate(SampleYear = year(Date))
                     

write.csv(trees3, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_trees_QC.csv", row.names = FALSE)

#Adding Core_ID for Jay so can join to core data and calculate BAI
trees4 <- trees3
trees4$Tag<-str_pad(trees4$Tag,width=3,side="left",pad=0) #Pad plot number so retains 3-digits
trees4$Core_ID<-paste(trees4$Site, trees4$Tag, sep="_")

write.csv(trees4, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_trees_QC_Core_ID.csv", row.names = FALSE)
# QC checks for sapling data ----------------------------------------------
str(saps)
#Site
names(saps)
unique(saps$Site) # 8 site codes

#Transect
unique(saps$Transect) # 4 transects

#Quadrat
unique(saps$Quadrat)# BM has a 3.9 quadrat? Data entry error, should be Quadrat 60 (checked datasheet). Fixed

#Species
unique(saps$Species)
table(saps$Species) #fix spelling errors
saps$Species <- recode(saps$Species, BETCOR = "BECO")
saps$Species <- recode(saps$Species, PRUPEN = "PRPE")

#DBH
unique(saps$DBH)#NAs are for quadrats with no saplings
not_sap <- saps %>% filter(DBH >= 10) # none >= 10 cm DBH
DBHs_plot <- ggplot(saps, aes(x = DBH))+
  geom_histogram(binwidth=1)+
  facet_wrap(~Site)+
  theme_FVM()
DBHs_plot

DBHs_plot2 <- ggplot(saps, aes(x = Site, y = DBH))+
  geom_boxplot()
DBHs_plot2

saps2 <- saps %>% mutate(Date = case_when(Site == 'WP2' ~ '2022-10-01',
                                          Site == 'IB' ~ '2021-07-07',
                                          Site == 'WP' ~ '2021-07-21',
                                          Site == 'OP' ~ '2020-08-06',
                                          Site == 'PM' ~ '2020-08-20',
                                          Site == 'BH' ~ '2021-06-18',
                                          Site == 'BC' ~ '2020-08-17',
                                          Site == 'BM' ~ '2020-08-13')) %>% 
                 mutate(SampleEventNum = 2) %>% 
                 mutate(SampleYear = year(Date)) 

write.csv(saps2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_saps_QC.csv", row.names = FALSE)

# QC checks for seedlings -------------------------------------------------
str(seeds)
#Site
names(seeds)
unique(seeds$Site) # 8 site codes

#Date
unique(seeds$Date)

#Quad
unique(seeds$Quadrat_Num)
table(seeds$Quadrat_Num)

#Latin_Name
unique(seeds$Latin_name)

#DS_quad
table(seeds$DS_quad)

#Locations_m
table(seeds$Location_m)

#Size class
table(seeds$Size_class)

#Count
table(seeds$Count)

table(complete.cases(seeds))# no na's

seeds2 <- seeds %>% select(-"DS_quad", -"Location_m") %>% mutate(SampleYear = year(Date)) %>% 
                    mutate(SampleEventNum = 2) 

write.csv(seeds2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_quad_seedings_QC.csv", row.names = FALSE)

# QC checks for quad species data -----------------------------------------
str(qd_sp)
#Site
names(qd_sp)
unique(qd_sp$Site) # 7 site codes

#Date
unique(qd_sp$Date)
qd_sp$Date <- mdy(qd_sp$Date)

#Latin_Name
unique(qd_sp$Latin_name)

#Cover data 
qd_sp %>% map(table) #iterate table() over all columns, prints output for each in the console

table(complete.cases(qd_sp))# no na's
qd_sp2 <- qd_sp %>% select(-"DS_quad", -"Location_m") %>% mutate(SampleYear = year(Date)) %>% 
  mutate(SampleEventNum = 2)

write.csv(qd_sp2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_quad_species_QC.csv", row.names = FALSE)
# QC checks for quad character data ---------------------------------------
str(qd_ch)
qd_ch %>% map(table) #iterate table() over all columns, prints output for each in the console

table(complete.cases(qd_ch))# no na's
qd_ch2 <- qd_ch %>% select(-"DS_quad", -"Location_m") %>% mutate(SampleYear = year(Date)) %>% 
  mutate(SampleEventNum = 2) 

write.csv(qd_ch2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_quad_character_QC.csv", row.names = FALSE)
# QC checks for Soil Depth --------------------------------------------------------------
str(soil)
names(soil)
soil1 <- soil %>% rename(Depth = `Soil Depth (cm)`) %>% select(-Notes)

soil1 %>% map(table)# IB has a depth of 100+, changing to 100 cm so numeric
soil1$Depth <- recode(soil1$Depth, '100+' = "100")

str(soil1)
soil1$Depth <- as.numeric(soil1$Depth)

depth_plot <- ggplot(soil1, aes(x = Site, y = Depth))+
                  geom_boxplot()
            
depth_plot

soil2 <- soil1 %>% select(-"DS point*") %>% mutate(SampleYear = year(Date)) 

table(complete.cases(soil2))# no na's
write.csv(soil2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_soil_depth_QC.csv", row.names = FALSE)

# QC checks for 1959 tree data --------------------------------------------
#Cover data 
str(trees_59)
trees_59 %>% map(table) #iterate table() over all columns, prints output for each in the console

table(complete.cases(trees_59))# no na's

trees_592 <- trees_59 %>% mutate(SampleYear = 1959) %>% 
  mutate(SampleEventNum = 1) 

write.csv(trees_592, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_trees_1959_QC.csv", row.names = FALSE)


