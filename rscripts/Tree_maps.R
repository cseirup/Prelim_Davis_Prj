# Exploratory plotting of Davis data

# Load packages -----------------------------------------------------------

library(tidyverse)
library(forestMIDNarch)

# Load data ---------------------------------------------------------------

saps <- read.csv("../Davis_data/Davis_saps_20211101.csv")
trees <- read.csv("../Davis_data/Davis_trees_20211101.csv")
#gr <- readxl::read_xlsx("../data/Davis_growth_releases.xlsx")
#CWD <- readxl::read_xlsx("../data/CWD_data.xlsx")

# Plotting ----------------------------------------------------------------
as.factor(trees$Site)
trees$Site <- ordered(trees$Site, levels = c("BM", "BC", "OP", "PM")) #order sites for consistency
levels(trees$Site)

site_names <- c('BM' = 'Beech Mtn', 'BC' = 'Blackwoods','OP' = 'Otter Point', 'PM' = "Pemetic Mtn")
Tree_layout <- ggplot(data = trees, aes(x = cY, y = cX))+
        geom_point(aes(color = Status, size = DBH))+
        #coord_equal()+# makes x and y the same scale
        scale_color_discrete(name = "Status", labels = c("Dead", "Live"))+
        scale_size_continuous(range = c(.5, 3))+
        facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names))+
        labs(x = "Y (meters)", y = "X (meters)")+ 
        theme(axis.text = element_text(size = 9), 
              strip.text = element_text(size = 12), #facet wrap text size
              axis.title = element_text(size = 12),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line = element_line(color = "#696969", size = 0.01),
              axis.ticks = element_line(color = "#696969", size = 0.5),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              axis.title.y = element_text(margin = margin(r = 5)))+
        theme_FVM()
Tree_layout
        
ggsave("./figures/Tree_layout15.jpg", Tree_layout, dpi = 300, height = 7, width = 12)

#plot individual site tree mapp
site_names <- c('BM' = 'Beech Mtn', 'BC' = 'Blackwoods','OP' = 'Otter Point', 'PM' = "Pemetic Mtn")
BM_trees <- trees %>% filter(Site == 'BM') %>% filter(cY >= 100)
tree_map <- ggplot(data = BM_trees, aes(x = cY, y = cX))+
        geom_point(aes(color = Species, shape = Status, size = DBH))+
        #coord_equal()+# makes x and y the same scale
        #scale_color_discrete(name = "Status", labels = c("Dead", "Live"))+
        scale_shape_manual(values = c(17,16))+
        scale_size_continuous(range = c(1, 6))+
        scale_color_manual(values = c("#c8598f", "#6882cf","#8a71c9", "#6ca74d", "#cc5a43", "#b68f40", "#49adad"))+
        facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names))+
        labs(x = "Y (meters)", y = "X (meters)")+ 
        theme(axis.text = element_text(size = 9), 
              strip.text = element_text(size = 12), #facet wrap text size
              axis.title = element_text(size = 12),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line = element_line(color = "#696969", size = 0.01),
              axis.ticks = element_line(color = "#696969", size = 0.5),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              axis.title.y = element_text(margin = margin(r = 5)))+
        geom_hline(yintercept = 5, linetype = 2)+
        #geom_vline(xintercept = 100, linetype = 2)+
        geom_text(aes(label=Tag), size = 3)+ #add tag # labels
        coord_flip()+
        theme_FVM()
tree_map

ggsave("./figures/BM_tree_map_north.jpg", tree_map, dpi = 300, height = 11, width = 3)


# Spatial tree cores ------------------------------------------------------
names(trees)
trees$Tag <- str_pad(trees$Tag, 3, side = "left", pad = "0")
trees_coreID <- trees %>% mutate(Core_ID = paste0(Site, "_", Tag))
head(trees_coreID)

names(gr)
gr2 <- select(gr, 1:4)
Sp_gr <- trees_coreID %>% left_join(gr2, by = "Core_ID") %>% distinct(Core_ID, .keep_all = TRUE) %>%  select(1:10,12:17)
names(Sp_gr)

OP_gr <- Sp_gr %>% filter(Site == 'OP') %>% replace_na(replace = list(Spatial_interest = 'E'))
# A= did count 2000s
# B = met threshold, didn't count 2000s
# C = didn't meet threshold, didn't count
# D = no 2000s bump
# E = no tree core

site_names <- c('BC' = 'Blackwoods', 'BM' = 'Beech Mtn', 'OP' = 'Otter Point', 'PM' = "Pemetic Mtn")
OPcore_layout <- ggplot(data = OP_gr, aes(x = cY, y = cX))+
        geom_point(aes(color = Spatial_interest, size = DBH))+
        scale_color_manual(values = c("blue", 'green', 'orange', 'red', 'grey'), labels = c("2000s GR counted",  "met threshold, didn't count", "didn't meet threshold", "no 2000s bump", "no tree core"))+
        #coord_equal()+# makes x and y the same scale+
        scale_size_continuous(range = c(.5, 3))+
        facet_wrap(~Site, ncol = 1, labeller = as_labeller(site_names))+
        labs(x = "Y (meters)", y = "X (meters)")+ 
        theme(legend.position = "bottom")+
        theme(axis.text = element_text(size = 6), 
              strip.text = element_text(size = 8), 
              axis.title = element_text(size = 9),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line = element_line(color = "#696969", size = 0.01),
              axis.ticks = element_line(color = "#696969", size = 0.5),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              axis.title.y = element_text(margin = margin(r = 5)))+
        theme_FVM()
OPcore_layout

ggsave("./figures/OPcore_layout.jpg", OPcore_layout, height = 3, width = 12, dpi = 300)


# CWD spatial -------------------------------------------------------------
names(CWD)
CWD <- select(CWD, 1:7)
OPcwd <- filter(CWD, Site == "OP")
OPcwd2 <- OPcwd %>% group_by(Point, DC) %>% summarise(Num_cwd = length(DC))
OPcwd2$Point <-as.character(OPcwd2$Point)
OPcwd2$DC <-as.character(OPcwd2$DC)
OPcwd2$Point <- ordered(OPcwd2$Point, levels = c("0", '25', '50', '75', '100', '125', '150', '175', '200'))


OPcwdPlot<-ggplot(OPcwd2, aes(x=Point,  y= Num_cwd, fill = DC))+
        geom_bar(stat = "identity")+#, position = position_dodge()) +
        scale_fill_brewer()+
        scale_y_continuous(breaks = seq(0,10,1))+
        labs(x='Turkey Foot Location', y='Number of CWD pieces', title = 'OP: CWD by turkey foot and decay class')+ 
        theme_bw() + #scale_x_binned(show.limits = TRUE)+ 
        theme(axis.text.x=element_text(angle=30,hjust=1, size = 6), 
              axis.ticks = element_line(),
              axis.text.y = element_text(size = 6),
              axis.title.y = element_text(size = 8),
              axis.title.x = element_text(size = 8),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) 

OPcwdPlot
ggsave(file="./figures/OP_cwd_plot.jpg", OPcwdPlot, dpi=300, width=5, height=3, units='in')



