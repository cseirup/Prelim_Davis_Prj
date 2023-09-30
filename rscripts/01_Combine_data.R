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

BC_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_BC_trees_20200817.xlsx')
BM_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_BM_trees_20200813.xlsx')
OP_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_OP_trees_20210714.xlsx')
PM_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_PM_trees_20200808.xlsx')
WP_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP_trees_20210808.xlsx')
IB_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_IB_trees_20210714.xlsx')
BH_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_BH_trees_20210618.xlsx')
WP2_trees<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP2_trees_20221001.xlsx')

BC_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_BC_saps_20200817.xlsx')
BM_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_BM_saps_20200813.xlsx')
OP_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_OP_saps_20200806.xlsx')
PM_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_PM_saps_20200808.xlsx')
WP_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP_saps_20210808.xlsx')
IB_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_IB_saps_20210711.xlsx')
BH_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_BH_saps_20210618.xlsx')
WP2_saps<- read_excel('C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP2_saps_20221026.xlsx')

seeds <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_quad_seeding_data_20221011.xlsx")
qd_sp <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_quad_species_data_20221011.xlsx")
qd_ch <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_quad_character_data_20221011.xlsx")

WP2_seeds <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP2_quad_seeding_20221209.xlsx")
WP2_qd_sp <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP2_quad_species_20221209.xlsx")
WP2_qd_ch <- read_excel("C:/Users/cseirup/Documents/Personal/grad school/Davis project/data - raw/Davis_WP2_quad_character_20221209.xlsx")


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
IB_trees1 <- IB_trees %>% select(Tr, TAG, Status, SPP, DBH, X, Y, Crown, Comments) %>%  
  add_column(Site = 'IB', .before = 'Tr') # remove DS trans column
BH_trees1 <- BH_trees %>% add_column(Site = 'BH', .before = 'Tr')
WP2_trees1 <- WP2_trees %>% add_column(Site = 'WP2', .before = 'Tr')

BC_saps1 <- BC_saps %>% add_column(Site = 'BC', .before = 'Trans')
BM_saps1 <- BM_saps %>% add_column(Site = 'BM', .before = 'Trans')
OP_saps1 <- OP_saps %>% add_column(Site = 'OP', .before = 'Trans')
PM_saps1 <- PM_saps %>% add_column(Site = 'PM', .before = 'Trans')
WP_saps1 <- WP_saps %>% add_column(Site = 'WP', .before = 'Trans')  
IB_saps1 <- IB_saps %>% select(Trans, Quad, SPP, DBH) %>% add_column(Site = 'IB', .before = 'Trans') #remove DS trans column
BH_saps1 <- BH_saps %>% add_column(Site = 'BH', .before = 'Trans')
WP_saps1 <- WP_saps %>% add_column(Site = 'WP', .before = 'Trans')
WP2_saps1 <- WP2_saps %>% add_column(Site = 'WP2', .before = 'Trans')

# Combine all except WP, WP2
Trees <- rbind(BC_trees1, BM_trees1, OP_trees1, PM_trees1, IB_trees1, BH_trees1) 
Saps <- rbind(BC_saps1, BM_saps1, OP_saps1, PM_saps1, IB_saps1, BH_saps1, WP_saps1)# WP2 not included 

#Calculate x and y values for whole plot based on transect number
Trees1 <- Trees %>% mutate(cX = case_when(Tr == 1 ~ X+0, Tr == 2 ~ X+5, Tr == 3 ~ X+0, Tr == 4 ~ X+5 )) %>%                      
                     mutate(cY = case_when(Tr == 1 ~ Y+0, Tr == 2 ~ Y+0, Tr == 3 ~ Y+100, Tr == 4 ~ Y+100))

#Calculate X values for WP - 20 x 100 plot sampled instead of 10 x 200 plot
WP_trees2 <- WP_trees1 %>% mutate(cX = case_when(Tr == 1 ~ X+0, Tr == 2 ~ X+5, Tr == 3 ~ X+10, Tr == 4 ~ X+15 )) %>%  # 1 to 4 arranged from S to N, W to E
                            mutate(cY = case_when(Tr == 1 ~ Y+0, Tr == 2 ~ Y+0, Tr == 3 ~ Y+0, Tr == 4 ~ Y+0))

#Calculate WP2 x and y values for whole plot based on transect number
#Plot was sampled from 0 to 170m than 0 to -30m because the plot because to steep to safely sample at 170m
#Converting original x, y to the corresponding values had the plot be sampled in the traditional way (0m to 200m)
WP2_trees2 <- WP2_trees1 %>% mutate(DS_Tr = Tr) %>% 
                             mutate(cX = case_when(Tr == 1 ~ X+0, Tr == 2 ~ X+5, 
                                                   Tr == 3 ~ X+0, Tr == 4 ~ X+5,
                                                   Tr == 5 ~ X+0, Tr == 6 ~ X+5)) %>%  
                             mutate(cY = case_when(Tr == 5 ~ 30-Y, Tr == 6 ~ 30-Y,
                                                   Tr == 1 ~ Y+30, Tr == 2 ~ Y+30, 
                                                   Tr == 3 ~ Y+130, Tr == 4 ~ Y+130))
#Converting original transect #'s to the corresponding values had the plot be sampled in the traditional way (0m to 200m)
WP2_trees3 <- WP2_trees2 %>% mutate(Tr = case_when(cY <= 100 & cX <= 5 ~ 1, 
                                                   cY > 100 & cX <= 5 ~ 3, 
                                                   cY <= 100 & cX > 5 ~ 2,
                                                   cY > 100 & cX > 5 ~ 4))

WP2_trees4 <- WP2_trees3 %>%  select(-DS_Tr)#removed original datasheet transect #'s.


#Converting sapling Quad # + Transect #to the corresponding #'s if had it been sampled in the traditional way
WP2_saps2 <- WP2_saps1 %>% mutate(DS_Quad = Quad) %>%
                          mutate(DS_Trans = Trans) %>%
  mutate(Quad = case_when(Trans == 1 ~ Quad+30, Trans == 2 ~ Quad+30, 
                          Trans == 3 ~ Quad+30, Trans == 4 ~ Quad+30,
                          Trans == 5 ~ 20-Quad, Trans == 6 ~ 20-Quad))

WP2_saps3 <- WP2_saps2 %>% mutate(Trans = case_when(Quad >= 100 & DS_Trans == 1 ~ 3, 
                                                   Quad >= 100 & DS_Trans == 2 ~ 4, 
                                                   DS_Trans == 5 ~ 1, 
                                                   DS_Trans == 6 ~ 2, 
                                                   TRUE ~ Trans))

WP2_saps4 <- WP2_saps3 %>% mutate(Quad = case_when(Quad >= 100 & DS_Trans == 1 ~ Quad-100, 
                                                   Quad >= 100 & DS_Trans == 2 ~ Quad-100,
                                                   TRUE ~ Quad))

plot(WP2_saps4$Trans, WP2_saps4$Quad)
plot(WP2_saps4$DS_Trans, WP2_saps4$DS_Quad)

WP2_saps5 <- WP2_saps4 %>%  select(-DS_Trans, -DS_Quad)

#combine WP, WP2 with other plots
Trees2 <- rbind(Trees1, WP_trees2, WP2_trees4)
Saps1 <- rbind(Saps, WP2_saps5)


# Standardize column names
names(Trees2)
names(Saps)
Trees3 <- Trees2 %>% rename(Transect = Tr, Tag = TAG, Species = SPP) 
Trees3$Tag <- str_pad(Trees3$Tag, 3, side = "left", pad = "0")
Saps2 <- Saps1 %>% rename(Transect = Trans, Quadrat = Quad, Species = SPP)

table(complete.cases(Saps2$DBH))# Quadrats with no saps are entered as Species = 'NA' and DBH = 'NA'. Could change to NO SPP and leave DBH null?

# Export combined datasets as .csv
write.csv(Trees3, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_trees_20221206.csv", row.names = FALSE)
write.csv(Saps2, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_saps_20221206.csv", row.names = FALSE) 
# Saplings not transformed to interpret spatial info included in transect and quadrat number. Remember that WP plot layout is diff


# Combine Quadrat Data ----------------------------------------------------
#move quadfiles to raw data folder
#transform to long
#add column for quadrat location
#Rename quadrats 1-30
#combine datasets

#WP2 plot was sampled out of order: 30 to 170m, than 30 to 0m. Quadrats were not placed at the same locations in the final plot as other
#Davis plots, but in an unbiased way. This code will take the original quadrat names recorded on datasheets and convert them to the
#real location in meters in the final plot. A number 1-30 will be assigned to each quad in ascending order


#Quadrat Character Data
  #WP2 plot first
WP2_qd_ch1 <- WP2_qd_ch %>% pivot_longer(4:33, names_to = "DS_quad", values_to = "Cover")
WP2_qd_ch2 <- WP2_qd_ch1 %>% mutate(Location_m = DS_quad)
WP2_qd_ch3 <- WP2_qd_ch2 %>% mutate(Location_m = case_when(DS_quad == 'neg5m' ~ 25, 
                                                        DS_quad == 'neg15m' ~ 15,
                                                        DS_quad == 'neg16m' ~ 14,
                                                        DS_quad == 'neg25m' ~ 5,
                                                        DS_quad == '4m' ~ 34,
                                                        DS_quad == '5m' ~ 35,
                                                        DS_quad == '14m' ~ 44,
                                                        DS_quad == '24m' ~ 54,
                                                        DS_quad == '25m' ~ 55,
                                                        DS_quad == '34m' ~ 64,
                                                        DS_quad == '44m' ~ 74,
                                                        DS_quad == '45m' ~ 75,
                                                        DS_quad == '54m' ~ 84,
                                                        DS_quad == '64m' ~ 94,
                                                        DS_quad == '65m' ~ 95,
                                                        DS_quad == '74m' ~ 104,
                                                        DS_quad == '84m' ~ 114,
                                                        DS_quad == '85m' ~ 115,
                                                        DS_quad == '94m' ~ 124,
                                                        DS_quad == '104m' ~ 134,
                                                        DS_quad == '105m' ~ 135,
                                                        DS_quad == '114m' ~ 144,
                                                        DS_quad == '124m' ~ 154,
                                                        DS_quad == '125m' ~ 155,
                                                        DS_quad == '134m' ~ 164,
                                                        DS_quad == '144m' ~ 174,
                                                        DS_quad == '145m' ~ 175,
                                                        DS_quad == '154m' ~ 184,
                                                        DS_quad == '164m' ~ 194,
                                                        DS_quad == '165m' ~ 195))

WP2_qd_ch4 <- WP2_qd_ch3 %>% arrange(Character, Location_m) %>% 
                             add_column(Quadrat_Num = rep(1:30, times = 7))
#All other sites
qd_ch1 <- qd_ch %>% pivot_longer(4:33, names_to = "DS_quad", values_to = "Cover")
qd_ch2 <- qd_ch1 %>% mutate(Location_m = DS_quad)
qd_ch2$Location_m <- str_replace(qd_ch2$Location_m, "m", "")

str(qd_ch2)
qd_ch2$Location_m <- as.numeric(qd_ch2$Location_m)
qd_ch3 <- qd_ch2 %>% arrange(Site, Character, Location_m) %>% 
  add_column(Quadrat_Num = rep(1:30, times = 49))

names(qd_ch3)
names(WP2_qd_ch4)
qd_ch4 <-rbind(WP2_qd_ch4, qd_ch3)

write.csv(qd_ch4, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_char_long_20221209.csv", row.names = FALSE) 


#Quadrat Species Data
#WP2 plot first
WP2_qd_sp1 <- WP2_qd_sp %>% pivot_longer(4:33, names_to = "DS_quad", values_to = "Cover")
WP2_qd_sp2 <- WP2_qd_sp1 %>% mutate(Location_m = DS_quad)
WP2_qd_sp3 <- WP2_qd_sp2 %>% mutate(Location_m = case_when(DS_quad == 'neg5m' ~ 25, 
                                                           DS_quad == 'neg15m' ~ 15,
                                                           DS_quad == 'neg16m' ~ 14,
                                                           DS_quad == 'neg25m' ~ 5,
                                                           DS_quad == '4m' ~ 34,
                                                           DS_quad == '5m' ~ 35,
                                                           DS_quad == '14m' ~ 44,
                                                           DS_quad == '24m' ~ 54,
                                                           DS_quad == '25m' ~ 55,
                                                           DS_quad == '34m' ~ 64,
                                                           DS_quad == '44m' ~ 74,
                                                           DS_quad == '45m' ~ 75,
                                                           DS_quad == '54m' ~ 84,
                                                           DS_quad == '64m' ~ 94,
                                                           DS_quad == '65m' ~ 95,
                                                           DS_quad == '74m' ~ 104,
                                                           DS_quad == '84m' ~ 114,
                                                           DS_quad == '85m' ~ 115,
                                                           DS_quad == '94m' ~ 124,
                                                           DS_quad == '104m' ~ 134,
                                                           DS_quad == '105m' ~ 135,
                                                           DS_quad == '114m' ~ 144,
                                                           DS_quad == '124m' ~ 154,
                                                           DS_quad == '125m' ~ 155,
                                                           DS_quad == '134m' ~ 164,
                                                           DS_quad == '144m' ~ 174,
                                                           DS_quad == '145m' ~ 175,
                                                           DS_quad == '154m' ~ 184,
                                                           DS_quad == '164m' ~ 194,
                                                           DS_quad == '165m' ~ 195))

WP2_qd_sp4 <- WP2_qd_sp3 %>% arrange(Latin_name, Location_m) %>% 
  add_column(Quadrat_Num = rep(1:30, times = 13))
#All other sites
qd_sp1 <- qd_sp %>% pivot_longer(4:33, names_to = "DS_quad", values_to = "Cover")
qd_sp2 <- qd_sp1 %>% mutate(Location_m = DS_quad)
qd_sp2$Location_m <- str_replace(qd_sp2$Location_m, "m", "")

str(qd_sp2)
qd_sp2$Location_m <- as.numeric(qd_sp2$Location_m)
qd_sp3 <- qd_sp2 %>% arrange(Site, Latin_name, Location_m) %>% 
  add_column(Quadrat_Num = rep(1:30, times = 85))

names(qd_sp3)
names(WP2_qd_sp4)
qd_sp4 <-rbind(WP2_qd_sp4, qd_sp3)

write.csv(qd_sp4, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_species_long_20221209.csv", row.names = FALSE) 


#Quadrat Seedling Data
#WP2 plot first
WP2_seeds1 <- WP2_seeds %>% pivot_longer(5:9, names_to = "Size_class", values_to = "Count")

WP2_seeds2 <- WP2_seeds1 %>% rename(DS_quad = Quad) %>%  mutate(Location_m = DS_quad)
WP2_seeds3 <- WP2_seeds2 %>% mutate(Location_m = case_when(DS_quad == '-5' ~ 25, 
                                                           DS_quad == '-15' ~ 15,
                                                           DS_quad == '-16' ~ 14,
                                                           DS_quad == '-25' ~ 5,
                                                           DS_quad == '4m' ~ 34,
                                                           DS_quad == '5m' ~ 35,
                                                           DS_quad == '14m' ~ 44,
                                                           DS_quad == '24m' ~ 54,
                                                           DS_quad == '25m' ~ 55,
                                                           DS_quad == '34m' ~ 64,
                                                           DS_quad == '44m' ~ 74,
                                                           DS_quad == '45m' ~ 75,
                                                           DS_quad == '54m' ~ 84,
                                                           DS_quad == '64m' ~ 94,
                                                           DS_quad == '65m' ~ 95,
                                                           DS_quad == '74m' ~ 104,
                                                           DS_quad == '84m' ~ 114,
                                                           DS_quad == '85m' ~ 115,
                                                           DS_quad == '94m' ~ 124,
                                                           DS_quad == '104m' ~ 134,
                                                           DS_quad == '105m' ~ 135,
                                                           DS_quad == '114m' ~ 144,
                                                           DS_quad == '124m' ~ 154,
                                                           DS_quad == '125m' ~ 155,
                                                           DS_quad == '134m' ~ 164,
                                                           DS_quad == '144m' ~ 174,
                                                           DS_quad == '145m' ~ 175,
                                                           DS_quad == '154m' ~ 184,
                                                           DS_quad == '164m' ~ 194,
                                                           DS_quad == '165m' ~ 195))

str(WP2_seeds3)
WP2_seeds4 <- WP2_seeds3 %>% arrange(Location_m, Latin_name) # no change, already in the right order
qnum <- WP2_seeds4 %>% select(Site, Location_m) %>% 
                          unique() %>%  add_column(Quadrat_Num = seq(1:30))
WP2_seeds5 <- left_join(WP2_seeds4, qnum, by = c("Site","Location_m"))

#All other sites
seeds1 <- seeds %>% pivot_longer(5:9, names_to = "Size_class", values_to = "Count")
seeds2 <- seeds1 %>% rename(DS_quad = Quad) %>% mutate(Location_m = DS_quad)
seeds2$Location_m <- str_replace(seeds2$Location_m, "m", "")

str(seeds2)
seeds2$Location_m <- as.numeric(seeds2$Location_m)
qnum2 <- seeds2 %>% select(Site, Location_m) %>% 
                       unique() %>%  add_column(Quadrat_Num = rep(1:30, times = 7))
seeds3 <- left_join(seeds2, qnum2, by = c("Site","Location_m"))

names(seeds3)
names(WP2_seeds5)
seeds4 <-rbind(WP2_seeds5, seeds3)

write.csv(seeds4, "C:/Users/cseirup/Documents/Personal/grad school/Davis project/data/Davis_quad_seedlings_long_20221209.csv", row.names = FALSE) 


