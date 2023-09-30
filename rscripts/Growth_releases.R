install.packages("tidyverse")
library(tidyverse)
library(dplR)
#library(TRADER)

# Load ring width data ----------------------------------------------------
BC <- read.tucson('../Davis_data/PIRU_best/BC_best.rwl')
BM <- read.tucson('../Davis_data/PIRU_best/BM_best.rwl')
PM <- read.tucson('../Davis_data/PIRU_best/PM_best.rwl')
OP <- read.tucson('../Davis_data/PIRU_best/OP_best.rwl')
WMA <- read.fh('../Davis_data/PIRU_best/WMALL_best_stripped.fh')
IB <- read.fh('../Davis_data/PIRU_best/IB_best_stripped.fh')
BH <- read.fh('../Davis_data/PIRU_best/BH_best_stripped.fh')
PM_TSCA <- read.fh('../Davis_data/PIRU_best/PM_TSCA_best_stripped.fh')


Cores <- read.csv('../Davis_data/Davis_tree_cores_20230405.csv')#questionable use of rings to pith from this file
Cores$Core <- str_replace(Cores$Core, "-","_")
Cores <- rename(Cores, Core_ID = Core)

#Use years to pith from CDendro attribute output
BC_attr <- read.table('../Davis_data/PIRU_best/BC_best_STRIPPED_attributes.txt', header = TRUE)#CDENDRO outputs NULL d2p/y2p as 0, change to NA + remove header manually before bringing into R.
BM_attr <- read.table('../Davis_data/PIRU_best/BM_best_STRIPPED_attributes.txt', header = TRUE)
PM_attr <- read.table('../Davis_data/PIRU_best/PM_best_STRIPPED_attributes.txt', header = TRUE)
OP_attr <- read.table('../Davis_data/PIRU_best/OP_STRIPPED_attributes.txt', header = TRUE)
WMA_attr <- read.table('../Davis_data/PIRU_best/WMALL_best_STRIPPED_attributes.txt', header = TRUE)
IB_attr <- read.table('../Davis_data/PIRU_best/IB_best_STRIPPED_attributes.txt', header = TRUE)
BH_attr <- read.table('../Davis_data/PIRU_best/BH_best_STRIPPED_attributes.txt', header = TRUE)
PM_TSCA_attr <- read.table('../Davis_data/PIRU_best/PM_TSCA_best_STRIPPED_attributes.txt', header = TRUE)

attr <- bind_rows(BC_attr, BM_attr, PM_attr, OP_attr, WMA_attr, IB_attr, BH_attr, PM_TSCA_attr)

# Create analysis dataset -------------------------------------------------
#combine ring width data for all sites
BC1 <- rownames_to_column(BC, var = "year")
BM1 <- rownames_to_column(BM, var = "year")
PM1 <- rownames_to_column(PM, var = "year")
OP1 <- rownames_to_column(OP, var = "year")
WMA1 <- rownames_to_column(WMA, var = "year")
IB1 <- rownames_to_column(IB, var = "year")
BH1 <- rownames_to_column(BH, var = "year")
PM_TSCA1 <- rownames_to_column(PM_TSCA, var = "year")


All_rw <- full_join(BC1, BM1, by = "year")
All_rw2 <- full_join(All_rw, PM1, by = "year")
All_rw3 <- full_join(All_rw2, OP1, by = "year")
All_rw4 <- full_join(All_rw3, WMA1, by = "year")
All_rw5 <- full_join(All_rw4, IB1, by = "year")
All_rw6 <- full_join(All_rw5, BH1, by = "year")
All_rw7 <- full_join(All_rw6, PM_TSCA1, by = "year")

#setting up attributes to include first/last years, length of series
rw7 <- column_to_rownames(All_rw7, var = "year")
stats <- rwl.stats(rw7)
attr2 <- stats %>% full_join(attr, by = "series") %>% select("series", "first", "last", "year", "d2pith", "years2pith")
attr2 <- rename(attr2, Core_ID = series)
attr2 <- rename(attr2, length = year)

#Reshape from wide to long
names(All_rw7)
All_long <-All_rw7 %>% pivot_longer(cols = c(BC_003:PM_213), names_to = "Core_ID",
               values_to = "Ring_width")

#All_long2 <- drop_na(All_long) #for Jay
#write.csv(All_long2, "C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_rseries_long_20221003.csv", row.names = FALSE )

#combine with core information
rgwd <- left_join(All_long, Cores, by = "Core_ID")
rgwd.a <- left_join(rgwd, attr2, by = "Core_ID")
names(rgwd.a)
rgwd.a <- rgwd.a %>% select("Site", "Core_ID", "year", "Ring_width", "Species", "first", "last", "length", "years2pith", "d2pith") %>% arrange(Core_ID, year)

# Creating lagged and leading columns -------------------------------------------------
rgwd1 <- rgwd.a %>% mutate(lag1 = lag(Ring_width, n = 1)) %>% 
                   mutate(lag2 = lag(Ring_width, n = 2)) %>% 
                    mutate(lag3 = lag(Ring_width, n = 3)) %>%
                     mutate(lag4 = lag(Ring_width, n = 4)) %>%
                      mutate(lag5 = lag(Ring_width, n = 5)) %>%
                       mutate(lag6 = lag(Ring_width, n = 6)) %>%
                        mutate(lag7 = lag(Ring_width, n = 7)) %>%
                         mutate(lag8 = lag(Ring_width, n = 8)) %>%
                          mutate(lag9 = lag(Ring_width, n = 9)) %>%
                           mutate(lag10 = lag(Ring_width, n = 10)) %>% 
                            mutate(lg10mean = (lag1 + lag2 + lag3 + lag4 + 
                                                 lag5 + lag6 + lag7 + lag8 + lag9 + lag10)/10)

rgwd2 <- rgwd1 %>% mutate(lead1 = lead(Ring_width, n = 1)) %>% 
                    mutate(lead2 = lead(Ring_width, n = 2)) %>% 
                     mutate(lead3 = lead(Ring_width, n = 3)) %>%
                      mutate(lead4 = lead(Ring_width, n = 4)) %>%
                       mutate(lead5 = lead(Ring_width, n = 5)) %>%
                        mutate(lead6 = lead(Ring_width, n = 6)) %>%
                         mutate(lead7 = lead(Ring_width, n = 7)) %>%
                          mutate(lead8 = lead(Ring_width, n = 8)) %>%
                           mutate(lead9 = lead(Ring_width, n = 9)) %>%
                            mutate(lead10 = lead(Ring_width, n = 10)) %>% 
                             mutate(ld10mean = (lead1 + lead2 + lead3 + lead4 + lead5 + 
                                                  lead6 + lead7 + lead8 + lead9 + lead10)/10)

#Need to remove averages for first and last 10 years of each core: 
#complete by removing max and min 10 years -/+ rows with NAs in average columns 
max(rgwd2$year)
min(rgwd2$year)

rgwd3 <- rgwd2 %>% mutate(lg10mean = replace(lg10mean, year<1750, NA)) %>% 
                    mutate(ld10mean = replace(ld10mean, year>2010, NA)) %>% 
                     drop_na(Ring_width)

table(complete.cases(rgwd3$years2pith)) #making sure I still have ring to pith NAs

# Calculating % increase and absolute increase ----------------------------
names(rgwd3)
rgwd4 <- rgwd3 %>% mutate(Per_inc = (ld10mean/lg10mean), Abs_inc = (ld10mean - lg10mean)) %>% 
                    select(Site, Core_ID, year, Ring_width, Species, years2pith, first, last, length, Per_inc, Abs_inc)
                    
#reshape for plotting
rgwd5 <- rgwd4 %>% pivot_longer(cols = c(Ring_width, Per_inc, Abs_inc), names_to = "RW_types", values_to = "RW_values")

# Plotting ----------------------------------------------------------------
class(rgwd5$year)
rgwd5$year <- as.numeric(rgwd5$year) #covert year to numeric value
class(rgwd5$year)
Core_IDs <- unique(rgwd5$Core_ID) #create vector of core names

gr_plot<-function(df, ID){ # df = data frame, ID = Core_ID value
    df2<- filter(df, Core_ID == ID) #filters data frame by Core_ID
    ggplot(aes(x = year, y = RW_values, colour = RW_types), data = df2)+
     geom_line(stat = "identity")+
     scale_x_continuous(breaks = seq(1750, 2020, 10), expand = expansion(add = 2))+
     scale_colour_manual(values = c(Abs_inc = "#46ac19", Per_inc = '#8479ff', Ring_width = '#bb0028'))+
     labs(x = "Year", y = "Ring width")+
     ggtitle(paste("Site =", df2$Site," ", "Core ID = ", ID, " ", "Rings to pith = ", df2$years2pith, " Years: ", df2$first, "-", df2$last, "(", df2$length, ")"))+
     geom_hline(yintercept = 2, linetype = 'longdash', colour = '#8479ff', linewidth = 1)+
     geom_hline(yintercept = .58, linetype = 'twodash', colour = '#46ac19', linewidth = 1)+
     theme_bw()
} #function to create plot; df = dataframe created above; ID = Core to plot

gr_plot(rgwd5, 'PM_213')#select a specific core to view

gr_plotlist <- map(Core_IDs, ~gr_plot(rgwd5, ID =.x)) #apply gr_plot function to each Core in the Core_IDs list

gr_plotlist[[300]] # view a plot, in numerical order by Core_ID

#Save all plots to pdf for easy viewing
pdf("gr_plots_all_sites2.pdf", height = 7, width = 10)
gr_plotlist
dev.off()

#manually tallying gap recruits and growth releases.
#will need to join with dataset of all years + count releases per year. NAs to 0
#Fancy reach task: make viz and us brush to extract values
################################################

