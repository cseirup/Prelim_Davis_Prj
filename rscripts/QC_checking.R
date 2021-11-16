#Quality control (QC) of Davis site data collecting 2020/2021

library(tidyverse)
library(forestMIDN)

setwd('C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj')

# Load data ---------------------------------------------------------------
saps <- read.csv("../Davis_data/Davis_saps_20211116.csv")
trees <- read.csv("../Davis_data/Davis_trees_20211116.csv")
#still to check: CWD, Soil depth data, Quad Data, Quadrat Species, Seedlings, ibutton, digitized davis data


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


