# Code to combine datasheet data in order to generate analysis dataset

# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)

# Set up folder structure ------------------------------------------------

subfolders<-c("shapefiles", "figures", "tables", "rscripts", "chrono_graphs")

invisible(lapply(subfolders, function(x){
if(!dir.exists(x)){dir.create(paste(getwd(), x, sep = '/'))}
}))

#if(!dir.exists("../data")){dir.create('../data')} # creating a data folder for datasets when created

# Import datasheet data ---------------------------------------------------

BC_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_BC_trees_20200817.xlsx')
BM_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_BM_trees_20200813.xlsx')
OP_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_OP_trees_20210714.xlsx')
PM_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_PM_trees_20200808.xlsx')
WP_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_WP_trees_20210808.xlsx')
IB_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_IB_trees_20210714.xlsx')
BH_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_BH_trees_20210618.xlsx')

BC_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_BC_saps_20200817.xlsx')
BM_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_BM_saps_20200813.xlsx')
OP_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_OP_saps_20200806.xlsx')
PM_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_PM_saps_20200808.xlsx')
WP_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_WP_saps_20210808.xlsx')
IB_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_IB_saps_20210711.xlsx')
BH_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/datasheets/Davis_BH_saps_20210618.xlsx')

# Combine datasets --------------------------------------------------------
# Adding site column 
BC_trees1 <- BC_trees %>% add_column(Site = 'BC', .before = 'Tr')
BM_trees1 <- BM_trees %>% add_column(Site = 'BM', .before = 'Tr')
OP_trees1 <- OP_trees %>% add_column(Site = 'OP', .before = 'Tr')
PM_trees1 <- PM_trees %>% add_column(Site = 'PM', .before = 'Tr')
  #adjust PM Y coordinates over 205m instead of 200m  
  PM_trees1$Y <- round(PM_trees1$Y*.97087, digits = 2)
WP_trees1 <- WP_trees %>% add_column(Site = 'WP', .before = 'Tr')  
 #adjust WP Y coordinates over 104m instead of 100m  
  WP_trees1$Y <- round(WP_trees1$Y*.96154, digits = 2)
IB_trees1 <- IB_trees %>% select(Tr, TAG, Status, SPP, DBH, X, Y, Crown, Comments) %>%  add_column(Site = 'IB', .before = 'Tr') # remove DS trans column
BH_trees1 <- BH_trees %>% add_column(Site = 'BH', .before = 'Tr')

BC_saps1 <- BC_saps %>% add_column(Site = 'BC', .before = 'Trans')
BM_saps1 <- BM_saps %>% add_column(Site = 'BM', .before = 'Trans')
OP_saps1 <- OP_saps %>% add_column(Site = 'OP', .before = 'Trans')
PM_saps1 <- PM_saps %>% add_column(Site = 'PM', .before = 'Trans')
WP_saps1 <- WP_saps %>% add_column(Site = 'WP', .before = 'Trans')  
IB_saps1 <- IB_saps %>% select(Trans, Quad, SPP, DBH) %>% add_column(Site = 'IB', .before = 'Trans') #remove DS trans column
BH_saps1 <- BH_saps %>% add_column(Site = 'BH', .before = 'Trans')

# Combine all except WP
Trees <- rbind(BC_trees1, BM_trees1, OP_trees1, PM_trees1, IB_trees1, BH_trees1)
Saps <- rbind(BC_saps1, BM_saps1, OP_saps1, PM_saps1, IB_saps1, BH_saps1, WP_saps1)# WP included 

#Calculate x and y values for whole plot based on transect number
Trees1 <- Trees %>% mutate(cX = case_when(Tr == 1 ~ X+0, Tr == 2 ~ X+5, Tr == 3 ~ X+0, Tr == 4 ~ X+5 )) %>%                      
                     mutate(cY = case_when(Tr == 1 ~ Y+0, Tr == 2 ~ Y+0, Tr == 3 ~ Y+100, Tr == 4 ~ Y+100))

#Calculate X values for WP - 20 x 100 plot sampled instead of 10 x 200 plot
WP_trees2 <- WP_trees1 %>% mutate(cX = case_when(Tr == 1 ~ X+0, Tr == 2 ~ X+5, Tr == 3 ~ X+10, Tr == 4 ~ X+15 )) %>%  # 1 to 4 arranged from S to N, W to E
                            mutate(cY = case_when(Tr == 1 ~ Y+0, Tr == 2 ~ Y+0, Tr == 3 ~ Y+0, Tr == 4 ~ Y+0))

#combine WP with other plots
Trees2 <- rbind(Trees1, WP_trees2)

# Standardize column names
names(Trees2)
names(Saps)
Trees3 <- Trees2 %>% rename(Transect = Tr, Tag = TAG, Species = SPP) 
Trees3$Tag <- str_pad(Trees3$Tag, 3, side = "left", pad = "0")
Saps1 <- Saps %>% rename(Transect = Trans, Quadrat = Quad, Species = SPP)

# Remove ibuttons + stump (11 rows)
#table(complete.cases(Trees3$Species))
#Trees4 <- drop_na(Trees3, Species)
#table(complete.cases(Trees4$Species)) # moved this to QC checks
table(complete.cases(Saps1$DBH))# Quadrats with no saps are entered as Species = 'NA' and DBH = 'NA'. Could change to NO SPP and leave DBH null?

# Export combined datasets as .csv
write.csv(Trees3, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_trees_20211116.csv", row.names = FALSE)
write.csv(Saps1, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_saps_20211116.csv", row.names = FALSE) # Saplings not transformed to interpret spatial info included in transect and quadrat number. Remember that WP plot layout is diff


