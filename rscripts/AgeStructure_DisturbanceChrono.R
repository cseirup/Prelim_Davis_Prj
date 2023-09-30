# Dendrochronology packages
library(dplR)
#library(treeclim)
# Data management packages
library(tidyverse)
library(forestNETN)
#install.packages("ggtext")
library(ggplot2)
library(ggtext)

source("C:/01_NETN/Forest_Health/R_Dev/Prelim_Davis_Prj/rscripts/davis_functions.R")
# Load data ---------------------------------------------------------------
dCores <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_1959_tree_cores.xlsx")

gr <- readxl::read_xlsx("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_growth_releases.xlsx")
gr2 <- select(gr, 1:3)
gr2$Core_ID <- str_replace(gr2$Core_ID, "-","_")
Cores <- read.csv("C:/01_NETN/Forest_Health/R_Dev/Davis_data/Davis_tree_cores_20230405.csv")# questionable use of rings to pith from this file
Cores$Core <- str_replace(Cores$Core, "-","_")
Cores <- rename(Cores, Core_ID = Core)

#change core IDs for cores measured on velmex so Core IDs match Cdendro output
Cores$Core_ID <- recode(Cores$Core_ID, IB_045 = "IB_045v")
Cores$Core_ID <- recode(Cores$Core_ID, IB_075 = "IB_075v")
Cores$Core_ID <- recode(Cores$Core_ID, IB_077 = "IB_077v")
Cores$Core_ID <- recode(Cores$Core_ID, IB_108 = "IB_108v")
Cores$Core_ID <- recode(Cores$Core_ID, IB_146 = "IB_146v")
Cores$Core_ID <- recode(Cores$Core_ID, IB_165 = "IB_165v")
Cores$Core_ID <- recode(Cores$Core_ID, IB_202 = "IB_202v")
Cores$Core_ID <- recode(Cores$Core_ID, WP_083 = "WP_083rs")
Cores$Core_ID <- recode(Cores$Core_ID, WP2_095 = "WP2_095x")
Cores$Core_ID <- recode(Cores$Core_ID, WP2_236 = "WP2_236v")
#Pool all Western Mtn sites
unique(Cores$Site)
Cores$Site <- recode(Cores$Site, "Western Mtn Plot" = "Western Mtn")
Cores$Site <- recode(Cores$Site, "Western Mtn Plot 2" = "Western Mtn")

# Load ring width data ----------------------------------------------------
BC <- read.tucson('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/BC_best.rwl')
BM <- read.tucson('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/BM_best.rwl')
PM <- read.tucson('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/PM_best.rwl')
OP <- read.tucson('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/OP_best.rwl')
WMA <- read.fh('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/WMALL_best_stripped.fh')
IB <- read.fh('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/IB_best_stripped.fh')
BH <- read.fh('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/BH_best_stripped.fh')
PM_TSCA <- read.fh('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/PM_TSCA_best_stripped.fh')

#Use years to pith from CDendro attribute output: Cores measured on Velmex do not have rings to pith data in attribute files
BC_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/BC_best_STRIPPED_attributes.txt', header = TRUE)#CDENDRO outputs NULL d2p/y2p as 0, change to NA + remove header manually before bringing into R.
BM_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/BM_best_STRIPPED_attributes.txt', header = TRUE)
PM_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/PM_best_STRIPPED_attributes.txt', header = TRUE)
OP_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/OP_STRIPPED_attributes.txt', header = TRUE)
WMA_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/WMALL_best_STRIPPED_attributes.txt', header = TRUE)
IB_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/IB_best_STRIPPED_attributes.txt', header = TRUE)
BH_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/BH_best_STRIPPED_attributes.txt', header = TRUE)
PM_TSCA_attr <- read.table('C:/01_NETN/Forest_Health/R_Dev/Davis_data/PIRU_best/PM_TSCA_best_STRIPPED_attributes.txt', header = TRUE)

attr <- bind_rows(BC_attr, BM_attr, PM_attr, OP_attr, WMA_attr, IB_attr, BH_attr, PM_TSCA_attr)
attr <- rename(attr, Core_ID = series)

#adding missing rings to pith to attr file
Cores2 <- attr %>% left_join(Cores, by = "Core_ID") #left join with attr to only include cross-dated cores included in analysis
names(Cores2)
Cores3 <- Cores2 %>% mutate(diff = years2pith-Rings.to.pith) #correct any discrepancies

#WM_032 .pos missing rings to pith info in attr file. Updated .pos but did not rerun cdendro. Adding the data here
Cores3$d2pith[139] <- '27.3'
Cores3$years2pith[139] <- '9'  
Cores4 <- Cores3 %>%  select(1:6)

#Sample depth ------------------------------------------------------------
names(BC)
BClong <- BC %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(BC_003:BC_271), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
BC_smdp <- BClong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Blackwoods")

names(BM)
BMlong <- BM %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(BM_001:BM_X3), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
BM_smdp <- BMlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Beech Mtn")

names(PM)
PMlong <- PM %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(PM_012:PM_203), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
PM_smdp <- PMlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Pemetic Mtn")

names(OP)
OPlong <- OP %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(OP_004:OP_345), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
OP_smdp <- OPlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Otter Point")

names(IB)
IBlong <- IB %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(IB_010:IB_202v), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
IB_smdp <- IBlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Ironbound Island")

names(BH)
BHlong <- BH %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(BH_002:BH_206), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
BH_smdp <- BHlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Bass Harbor Head")

names(WMA)
WMAlong <- WMA %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(WM_002:WP2_251), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
WMA_smdp <- WMAlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Western Mtn")

names(PM_TSCA)
PM_TSCAlong <- PM_TSCA %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(PM_010:PM_213), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
PM_TSCA_smdp <- PM_TSCAlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Pemetic Mtn TSCA")

All_smdp <- bind_rows(BC_smdp, BM_smdp, PM_smdp, OP_smdp, IB_smdp, BH_smdp, WMA_smdp, PM_TSCA_smdp)



# Compile data for age structures ------------------------------------------
BC_stats<-rwl.stats(BC) # first and last year, distribution parameters
BM_stats<-rwl.stats(BM)
OP_stats<-rwl.stats(OP)
PM_stats<-rwl.stats(PM)
WMA_stats<-rwl.stats(WMA)
IB_stats<-rwl.stats(IB)
BH_stats<-rwl.stats(BH)

stats_all <- rbind(BC_stats, BM_stats, OP_stats, PM_stats, WMA_stats, IB_stats, BH_stats)
stats_all <- rename(stats_all, Core_ID = series)

#set cores >10 rings from the pith to NA so they are not included in age structure
Cores4$years2pith <- as.numeric(Cores4$years2pith)
Cores4$d2pith <- as.numeric(Cores4$d2pith)
Cores5 <- Cores4
table(complete.cases(Cores4$years2pith))
Cores5$years2pith[Cores5$years2pith > 10] <-NA
table(complete.cases(Cores5$years2pith))

#Combine series stats with pith info. PM TSCA excluded #calc. recruit age/year, drop series w/o pith info
age <- stats_all %>% left_join(Cores5, by = c("Core_ID"))%>% mutate(r_year = first - years2pith,) %>% 
                      mutate(r_age = year + years2pith) %>% drop_na() %>% 
                      mutate(SampleEventNum = 2)      
                      
age2 <- age %>% select(-c(mean, median, stdev, skew, gini, ar1, d2pith, years2pith))

#Davis 1959 tree core data: calculate year (length of series), r_age, r_year, to match my data.No years to pith info. 
d_age <- dCores %>% mutate(year = last-first+1) %>% mutate(r_year = first) %>% 
                        mutate(r_age = year) %>% drop_na() %>% 
                        mutate(SampleEventNum = 1)


#combine 1959 and 2020s core info
ageAll <- rbind(age2, d_age)
ageAll$SampleEventNum <- as.character(ageAll$SampleEventNum)
ageAll2 <- ageAll %>% mutate(speciesEvent = str_c(Species, SampleEventNum, sep = ""))
# Plot Recruitment Age ---------------------------------------------------------
custom_breaks <- seq(1720,2010,10)
ageD<-ggplot(ageAll2, aes(x=r_year, fill = speciesEvent, ))+
  geom_bar()+
  scale_fill_manual(name = "Year/Species Cored", labels = c("1959 *Picea glauca*", '1959 *Picea rubens*', "2020s *Picea rubens*"), 
                    values = c("PIRU1" = '#b38e80', "PIRU2" = '#75abbd', "PIGL1"= '#b3b380'))+
  labs(x='Recruitment Decade', y='Number of Trees')+ 
  scale_x_binned(show.limits = TRUE, limits = c(1720, 2010), breaks = custom_breaks,
                 labels = every_nth(custom_breaks, 4, inverse = TRUE))+ 
  facet_wrap(~Site, ncol = 4)+
  theme(axis.text.x=element_text(angle=30, hjust = 1, vjust = 1.1, size = 10), # change axis label size
        axis.text.y = element_text(size = 10), 
        strip.text = element_text(size = 12), # change facet text size
        axis.title = element_text(size = 12), # change axis title size
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = c(.95,.1),
        legend.justification = c(1,0),
        panel.grid.major = element_blank(), #indented from theme_FHM()
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = '#696969', fill = 'white', size = 0.4),
        plot.background = element_blank(),
        strip.background = element_rect(color = '#696969', fill = 'grey90', size = 0.4),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "#696969", size = 0.2),
        axis.line.y = element_line(color = "#696969", size = 0.2),
        axis.ticks.y = element_line(color = "#696969", size = 0.7),
        axis.ticks.x = element_line(color = "#696969", size = 0.3),
        axis.ticks.length.x = unit(.15, "cm"))+
  theme(legend.text = element_markdown())

ageD

#ageD_forjpeg<-ggplot(ageAll2, aes(x=r_year, fill = speciesEvent, ))+
  #geom_bar()+
  #scale_fill_manual(name = "Year/Species Cored", labels = c("1959 *Picea glauca*", '1959 *Picea rubens*', "2020s *Picea rubens*"), 
                     #values = c("PIRU1" = '#b38e80', "PIRU2" = '#75abbd', "PIGL1"= '#b3b380'))+
  #labs(x='Recruitment Decade', y='Number of Trees')+ 
  #scale_x_binned(show.limits = TRUE, limits = c(1720, 2010), breaks = custom_breaks,
                # labels = every_nth(custom_breaks, 4, inverse = TRUE))+ 
  #facet_wrap(~Site, ncol = 4)+
 # theme(axis.text.x=element_text(angle=30, hjust = 1, vjust = 1.4, size = 3), # change axis label size
     #   axis.text.y = element_text(size = 4), 
      #  strip.text = element_text(size = 6), # change facet text size
       # axis.title = element_text(size = 7), # change axis title size
        #axis.title.y = element_text(margin = margin(r = 5)),
        #legend.text = element_text(size = 5),
        #legend.title = element_text(size = 5),
        #legend.position = c(1,0),
        #legend.justification = c(1,0),
         # panel.grid.major = element_blank(), #indented from theme_FHM()
          #panel.grid.minor = element_blank(),
          #panel.background = element_rect(color = '#696969', fill = 'white', size = 0.4),
          #plot.background = element_blank(),
          #strip.background = element_rect(color = '#696969', fill = 'grey90', size = 0.4),
          #legend.key = element_blank(),
          #axis.line.x = element_line(color = "#696969", size = 0.2),
          #axis.line.y = element_line(color = "#696969", size = 0.2),
          #axis.ticks.y = element_line(color = "#696969", size = 0.2),
          #axis.ticks.x = element_line(color = "#696969", size = 0.15),
          #axis.ticks.length.x = unit(.05, "cm"))+
  #theme(legend.text = element_markdown())

#ageD_forjpeg
#ggsave(file="./figures/AgeDist_comb_event_species2.jpg", ageD_forjpeg, dpi=300, width=5, height=3, units='in')

# Combine and Plot Disturbance Chronology------------------------------------------------------------
gr3 <- gr2 %>% left_join(Cores4, by = "Core_ID")
names(gr3)
gr4<-gr3 %>% group_by(Site, year, Species, Release_type) %>% summarize(count = n()) %>% drop_na()
gr4$year <- as.numeric(gr4$year)

#Calculate total sample depth per site to get % sample depth per year
Tot_sd <- Cores4 %>% group_by(Site) %>% summarise(tot_sd = n())
All_smdp$year <- as.numeric(All_smdp$year)
All_smdp2 <- All_smdp %>% left_join(Tot_sd, by = ("Site")) %>% 
                          mutate(Per_sd = (sample_depth/tot_sd)*100)

gr5 <- left_join(All_smdp2, gr4, by = c("Site", "year"))
#Express gap recruit and growth release counts as a percent of sample depth
gr6 <- gr5 %>% mutate(Per_gr = (count/sample_depth)*100)

#Releases binned decadaly
DistPlot<-ggplot(gr5)+
  geom_line(stat = "identity", linewidth = .7, color = "black", aes(x = year, y = sample_depth))+
  theme_bw()+
  geom_bar(stat = 'identity', position = 'stack', aes(fill=Release_type, x=year, y=count))+
  labs(x='Growth releases', y='Number of trees')+ 
  facet_wrap(~Site, ncol = 2, scales = "free")+
  scale_x_binned(n.breaks = 20, nice.breaks = TRUE, limits = c(1740, 2022))+
  #scale_y_continuous(limits = c(0,35))+
  scale_fill_manual(values = c( "cornflowerblue","darkorange"), name = "Release type", labels = c("Gap recruit", "Growth release"))+ 
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

DistPlot

#ggsave("../figures/GR_smdp_plot.jpg", DistPlot, height = 8, width = 8, dpi = 300)

#Just Western Mtn
DistPlotWM <- gr5 %>% filter(Site == "Western Mtn") %>% ggplot()+
  geom_line(stat = "identity", linewidth = .7, color = "black", aes(x = year, y = sample_depth))+
  theme_bw()+
  geom_bar(stat = 'identity', position = 'stack', aes(fill=Release_type, x=year, y=count))+
  labs(x='Growth releases', y='Number of trees')+ 
  facet_wrap(~Site, ncol = 1, scales = "free")+
  scale_x_binned(n.breaks = 20, nice.breaks = TRUE, limits = c(1740, 2022))+
  #scale_y_continuous(limits = c(0,35))+
  scale_fill_manual(values = c( "cornflowerblue","darkorange"), name = "Release type", labels = c("Gap recruit", "Growth release"))+ 
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

DistPlotWM

#ggsave("../figures/GR_smdp_plotWM.jpg", DistPlotWM, height = 8, width = 8, dpi = 300)

#Annual bins instead of decadel
#add grey shading when sample depth falls below 5 cores
gr7 <- gr6 %>% group_by(Site) %>% filter(sample_depth < 5) %>% 
                                  mutate(XMAX = max(year),
                                         XMIN = min(year)) %>% 
                                  select(3,10,11) %>% 
                                  unique()

gr8 <- gr6 %>% left_join(gr7, by = "Site")

##################

DistPlotAnn<-ggplot(gr8)+
  #geom_vline(aes(xintercept = year, color = "Sample Depth < 5"), data = gr7)+
  geom_rect(aes(ymax = 100, ymin = 0, xmin = XMIN, xmax = XMAX), fill="lightgrey", alpha = 0.5)+
  geom_line(stat = "identity", linetype = "dotted", linewidth = .9, aes(x = year, y = Per_sd, color = "% of Sample Depth"))+
  scale_colour_manual(values = c("black","grey")) +
  guides(colour = guide_legend(title = NULL))+
  theme_bw()+
  geom_bar(stat = 'identity', position = 'stack', aes(fill=Release_type, x=year, y=Per_gr))+
  labs(x='Year', y='Percent of trees showing a disturbance')+ 
  facet_wrap(~Site, ncol = 4, scales = "free")+
  scale_x_continuous(n.breaks = 20)+
  scale_fill_manual(values = c( "cornflowerblue","darkorange"), name = "Release type", labels = c("Gap recruit", "Growth release"))+ 
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

DistPlotAnn

#ggsave("../figures/GR_smdp_plot_annual.jpg", DistPlotAnn, height = 5, width = 10, dpi = 300)



#At temping base R plot
#set.seed(25852)                             # Create example data
#x <- rnorm(30)
#y1 <- x + rnorm(30)
#y2 <- x + rnorm(30, 5)

#names(gr8)
#grBC <- gr8 %>% filter(Site == "Blackwoods") %>% select(1,7,9)

#par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
#barplot(grBC$year, grBC$Per_gr, xlim = c(1750, 2022), ylim = c(0,100))              # Create first plot

#par(new = TRUE)                             # Add new plot
#plot(x, y2, pch = 17, col = 3,              # Create second plot without axes
    # axes = FALSE, xlab = "", ylab = "")
#axis(side = 4, at = pretty(range(y2)))      # Add second axis
#mtext("y2", side = 4, line = 3)             # Add second axis label



































# Combining Pemetic Mtn PIRU and TSCA -------------------------------------

gr3.5 <- gr2 %>% left_join(Cores4, by = "Core_ID")
gr3.5$Site <- recode(gr3.5$Site, "Pemetic Mtn TSCA" = "Pemetic Mtn")
gr4.5<-gr3.5 %>% group_by(Site, year, Release_type) %>% summarize(count = n()) %>% drop_na()
gr4.5$year <- as.numeric(gr4.5$year)

#Calculate total sample depth per site to get % sample depth per year
Cores4.5 <- Cores4
Cores4.5$Site <- recode(Cores4.5$Site, "Pemetic Mtn TSCA" = "Pemetic Mtn")
Tot_sd.5 <- Cores4.5 %>% group_by(Site) %>% summarise(tot_sd = n())

All_smdp.5 <- All_smdp
All_smdp.5$Site <- recode(All_smdp.5$Site, "Pemetic Mtn TSCA" = "Pemetic Mtn")
All_smdp.5$year <- as.numeric(All_smdp.5$year)
All_smdp1.5 <- All_smdp.5 %>% group_by(year, Site) %>% summarise(sample_depthPM = sum(sample_depth))
All_smdp2.5 <- All_smdp1.5 %>% left_join(Tot_sd.5, by = ("Site")) %>% 
                                mutate(Per_sd = (sample_depthPM/tot_sd)*100)

gr5.5 <- left_join(All_smdp2.5, gr4.5, by = c("Site", "year"))
#Express gap recruit and growth release counts as a percent of sample depth
gr6.5 <- gr5.5 %>% mutate(Per_gr = (count/sample_depthPM)*100)

#Annual bins instead of decadal
gr7.5 <- gr6.5 %>% group_by(Site) %>% filter(sample_depthPM < 5) %>% 
                                      mutate(XMAX = max(year),
                                             XMIN = min(year)) %>% 
                                     select(2,9,10) %>% 
                                     unique()

gr8.5 <- gr6.5 %>% left_join(gr7.5, by = "Site")

DistAnnPMpool <-ggplot(gr8.5)+
  geom_rect(aes(ymax = 100, ymin = 0, xmin = XMIN, xmax = XMAX), fill="lightgrey", alpha = 0.5)+
  geom_line(stat = "identity", linetype = "dotted", linewidth = .9, aes(x = year, y = Per_sd, color = "% of Sample Depth"))+
  scale_colour_manual(values = c("black","grey")) +
  guides(colour = guide_legend(title = NULL))+
  theme_bw()+
  geom_bar(stat = 'identity', position = 'stack', aes(fill=Release_type, x=year, y=Per_gr))+
  labs(x='Year', y='Percent of trees showing a disturbance')+ 
  facet_wrap(~Site, ncol = 4, scales = "free")+
  scale_x_continuous(n.breaks = 20)+
  scale_fill_manual(values = c( "cornflowerblue","darkorange"), name = "Release type", labels = c("Gap recruit", "Growth release"))+ 
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(.95,.15),
        legend.justification = c(1,0)) 

DistAnnPMpool


ggsave("./figures/GR_annual_pooledPM.jpg", DistAnnPMpool, height = 5, width = 10, dpi = 300)


