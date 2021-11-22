#Quality control (QC) of Davis site data collecting 2020/2021

library(tidyverse)
library(forestMIDN)
library(readxl)

setwd('C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj')

# Load data ---------------------------------------------------------------
saps <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_saps_20211116.csv")
trees <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_trees_20211116.csv")
seeds <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_seeding_data_20210814.xlsx")
qd_sp <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_species_data_20210814.xlsx")
qd_ch <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_character_data_20210814.xlsx")
soil <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_soil_depth_data_20210808.xlsx")
ibutton_d <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/ibutton_dat_daily_2021-06-29.csv")
ibutton_m <- read.csv("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/ibutton_dat_monthly_2021-06-29.csv")
trees_59<- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/1959_trees_20210409.xlsx")
#still to check: CWD, ibutton, digitized davis data (trees done)

# QC checks for tree data -------------------------------------------------
#Site
names(trees)
unique(trees$Site) # 7 site codes

#Transect
unique(trees$Transect) # 4 transects

#Tag
unique(trees$Tag)# One stump <1.3m recorded -> delete
stump <- trees %>% filter(trees$Tag == "Stump")
trees1 <- anti_join(trees, stump) #removed stump from PM trees
ibuttons <- trees1 %>% filter(is.na(trees1$Tag))
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
trees2$Species <- recode(trees2$Species, ACAC = "ACER")

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

#cX and cY
unique(trees2$cX)
unique(trees2$cY)

cXcY_plot <- ggplot(trees2, aes(x=cX, y=cY))+
  geom_point()+
  theme_FVM()
cXcY_plot #looks good

write.csv(trees2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_trees_QC.csv", row.names = FALSE)

# QC checks for sapling data ----------------------------------------------
#Site
names(saps)
unique(saps$Site) # 7 site codes

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
unique(saps$DBH)
not_sap <- saps %>% filter(DBH >= 10) # none >= 10 cm DBH
DBHs_plot <- ggplot(saps, aes(x = DBH))+
  geom_histogram()+
  facet_wrap(~Site)+
  theme_FVM()
DBHs_plot

DBHs_plot2 <- ggplot(saps, aes(x = Site, y = DBH))+
  geom_boxplot()
DBHs_plot2

write.csv(saps, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_saps_QC.csv", row.names = FALSE)

# QC checks for seedlings -------------------------------------------------
#Site
names(seeds)
unique(seeds$Site) # 7 site codes

#Date
unique(seeds$Date)

#Quad
unique(seeds$Quad)
table(seeds$Quad)

#Latin_Name
unique(seeds$Latin_name)

#15-30cm
table(seeds$`15-30cm`)# confirmed outliers (25 and 27) are correct

#30-100cm
table(seeds$`30-100cm`)# confirmed outlier (15) is correct

#1-1.5m
table(seeds$`1-1.5m`)

#>1.5m
table(seeds$`>1.5m`)

#Saps >1cm, ≤2.5 cm
table(seeds$'Saps >1cm, ≤2.5 cm')

table(complete.cases(seeds))# no na's
write.csv(seeds, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_quad_seedings_QC.csv", row.names = FALSE)

# QC checks for quad species data -----------------------------------------
#Site
names(qd_sp)
unique(qd_sp$Site) # 7 site codes

#Date
unique(qd_sp$Date)

#Latin_Name
unique(qd_sp$Latin_name)

#Cover data 
qd_sp %>% map(table) #iterate table() over all columns, prints output for each in the console

table(complete.cases(qd_sp))# no na's
write.csv(qd_sp, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_quad_species_QC.csv", row.names = FALSE)
# QC checks for quad character data ---------------------------------------
qd_ch %>% map(table) #iterate table() over all columns, prints output for each in the console

table(complete.cases(qd_ch))# no na's
write.csv(qd_ch, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_quad_character_QC.csv", row.names = FALSE)
# QC checks for Soil Depth --------------------------------------------------------------
names(soil)
soil1 <- soil %>% rename(Depth = `Soil Depth (cm)`) %>% select(-Notes)

soil1 %>% map(table)# IB has a depth of 100+, changing to 100 cm so numeric
soil1$Depth <- recode(soil1$Depth, '100+' = "100")

str(soil1)
soil1$Depth <- as.numeric(soil1$Depth)

depth_plot <- ggplot(soil1, aes(x = Site, y = Depth))+
                  geom_boxplot()
            
depth_plot

table(complete.cases(soil1))# no na's
write.csv(soil1, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_soil_depth_QC.csv", row.names = FALSE)


# QC checks for 1959 tree data --------------------------------------------
#Cover data 
trees_59 %>% map(table) #iterate table() over all columns, prints output for each in the console

table(complete.cases(trees_59))# no na's
write.csv(trees_59, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_trees_1959_QC.csv", row.names = FALSE)
# QC checks for quad character data ---------------------------------------

