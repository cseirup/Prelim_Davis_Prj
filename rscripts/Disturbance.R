# Dendrochronology packages
#install.packages("treeclim")
library(dplR)
#library(treeclim)
# Data management packages
library(tidyverse)


# Load data ---------------------------------------------------------------

gr <- readxl::read_xlsx("../data/Davis_growth_releases.xlsx")
gr2 <- select(gr, 1:4)
Cores <- read.csv('../data/Davis_tree_cores_20210118.csv')# questionable use of rings to pith from this file
Cores$Core <- str_replace(Cores$Core, "-","_")
Cores <- rename(Cores, Core_ID = Core)

# Load ring width data ----------------------------------------------------
BC <- read.tucson('../data/PIRU_best/BC_best.rwl')
BM <- read.tucson('../data/PIRU_best/BM_best.rwl')
PM <- read.tucson('../data/PIRU_best/PM_best.rwl')
OP <- read.tucson('../data/PIRU_best/OP_best.rwl')
WM <- read.tucson('../data/PIRU_best/WM_best.rwl')


# Sample depth ------------------------------------------------------------
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

names(WM)
WMlong <- WM %>% rownames_to_column(var = "year") %>% pivot_longer(cols = c(WM_002:WM_047), names_to = "core_ID", values_to = "ring_width") %>% drop_na()
WM_smdp <- WMlong %>% group_by(year) %>% summarise(sample_depth = n()) %>% add_column(Site = "Western Mtn")

All_smdp <- bind_rows(BC_smdp, BM_smdp, PM_smdp, OP_smdp, WM_smdp)

# Combine and Plot ------------------------------------------------------------

gr3 <- gr2 %>% left_join(Cores, by = "Core_ID")  %>%  select(1:10)
names(gr3)
gr4<-gr3 %>% group_by(Site, year, Release_type) %>% summarize(count = n()) %>% drop_na()
gr4$year <- as.numeric(gr4$year)
All_smdp$year <- as.numeric(All_smdp$year)
gr5 <- left_join(All_smdp, gr4, by = c("Site", "year"))


DistPlot<-ggplot(gr5)+
  geom_line(stat = "identity", size = .7, color = "black", aes(x = year, y = sample_depth))+
  theme_bw()+
  geom_bar(stat = 'identity', position = 'stack', aes(fill=Release_type, x=year, y=count))+
  labs(x='Growth releases', y='Number of trees')+ 
  facet_wrap(~Site, ncol = 1, scales = "free")+
  scale_x_binned(n.breaks = 20, nice.breaks = TRUE, limits = c(1740, 2020))+
  scale_y_continuous(limits = c(0,35))+
  scale_fill_manual(values = c( "cornflowerblue","darkorange"), name = "Release type", labels = c("Gap recruit", "Growth release"))+ 
  theme(axis.text.x=element_text(angle=45,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

DistPlot

ggsave("./figures/Dist_smdp_plot2.jpg", DistPlot, height = 8, width = 8, dpi = 300)
