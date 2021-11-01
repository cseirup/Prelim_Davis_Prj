# Code to combine datasheet data in order to generate analysis dataset

# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)

# Set up folder structure ------------------------------------------------

subfolders<-c("shapefiles", "figures", "tables", "rscripts", "chrono_graphs")

invisible(lapply(subfolders, function(x){
if(!dir.exists(x)){dir.create(paste(getwd(), x, sep = '/'))}
}))

if(!dir.exists("../data")){dir.create('../data')} # creating a data folder for datasets when created

# Import datasheet data ---------------------------------------------------

BC_trees<- read_excel('../datasheets/Davis_BC_trees_20200817.xlsx')
BM_trees<- read_excel('../datasheets/Davis_BM_trees_20200813.xlsx')
OP_trees<- read_excel('../datasheets/Davis_OP_trees_20200806.xlsx')
PM_trees<- read_excel('../datasheets/Davis_PM_trees_20200808.xlsx')

BC_saps<- read_excel('../datasheets/Davis_BC_saps_20200817.xlsx')
BM_saps<- read_excel('../datasheets/Davis_BM_saps_20200813.xlsx')
OP_saps<- read_excel('../datasheets/Davis_OP_saps_20200806.xlsx')
PM_saps<- read_excel('../datasheets/Davis_PM_saps_20200808.xlsx')

# Combine datasets --------------------------------------------------------
# Adding site column 
BC_trees1 <- BC_trees %>% add_column(Site = 'BC', .before = 'Tr')
BM_trees1 <- BM_trees %>% add_column(Site = 'BM', .before = 'Tr')
OP_trees1 <- OP_trees %>% add_column(Site = 'OP', .before = 'Tr')
PM_trees1 <- PM_trees %>% add_column(Site = 'PM', .before = 'Tr')
  #adjust PM Y coordinates over 205m instead of 200m  
  PM_trees1$Y <- PM_trees1$Y*.97087

BC_saps1 <- BC_saps %>% add_column(Site = 'BC', .before = 'Trans')
BM_saps1 <- BM_saps %>% add_column(Site = 'BM', .before = 'Trans')
OP_saps1 <- OP_saps %>% add_column(Site = 'OP', .before = 'Trans')
PM_saps1 <- PM_saps %>% add_column(Site = 'PM', .before = 'Trans')

# Combine
Trees <- rbind(BC_trees1, BM_trees1, OP_trees1, PM_trees1)
Saps <- rbind(BC_saps1, BM_saps1, OP_saps1, PM_saps1)

# Standardize column names
names(Trees)
names(Saps)
Trees1 <- Trees %>% rename(Transect = Tr, Tag = TAG, Species = SPP) 
Trees1$Tag <- str_pad(Trees1$Tag, 3, side = "left", pad = "0")
Saps1 <- Saps %>% rename(Transect = Trans, Quadrat = Quad, Species = SPP)

#Calculate x and y values for whole plot based on transect number

Trees2 <- Trees1 %>% mutate(cX = case_when(Transect == 1 ~ X+0, Transect == 2 ~ X+5, Transect == 3 ~ X+0, Transect == 4 ~ X+5 ))                       
                   %>% mutate(cY = case_when(Transect == 1 ~ Y+0, Transect == 2 ~ Y+0, Transect == 3 ~ Y+100, Transect == 4 ~ Y+100))

# Remove ibuttons + stump (9 rows)
table(complete.cases(Trees2$Species))
Trees3 <- drop_na(Trees2, Species)
table(complete.cases(Trees3$Species))
table(complete.cases(Saps1))

# Export combined datasets as .csv
write.csv(Trees3, "../data/Davis_trees_20200825.csv")
write.csv(Saps1, "../data/Davis_saps_20200825.csv")


