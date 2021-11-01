#Davis project: signal sample depth, standardization, and chronologies

# Dendrochronology packages
#install.packages("TRADER")
library(dplR)
library(TRADER)
install.packages("na_if")
#library(treeclim)
# Data managment packages
library(tidyverse)

# Load ring width data ----------------------------------------------------
BC <- read.tucson('../data/PIRU_best/BC_best.rwl')
BM <- read.tucson('../data/PIRU_best/BM_best.rwl')
PM <- read.tucson('../data/PIRU_best/PM_best.rwl')
OP <- read.tucson('../data/PIRU_best/OP_best.rwl')
WM <- read.tucson('../data/PIRU_best/WM_best.rwl')

BC_attr <- read.table('../data/PIRU_best/BC_best_STRIPPED_attributes.txt', header = TRUE)#CDENDRO outputs NULL d2p/y2p as 0, change to NA + remove header manually before bringing into R.
BM_attr <- read.table('../data/PIRU_best/BM_best_STRIPPED_attributes.txt', header = TRUE)
PM_attr <- read.table('../data/PIRU_best/PM_best_STRIPPED_attributes.txt', header = TRUE)
OP_attr <- read.table('../data/PIRU_best/OP_STRIPPED_attributes.txt', header = TRUE)
WM_attr <- read.table('../data/PIRU_best/WM_STRIPPED_attributes.txt', header = TRUE)

# Summary Stats -----------------------------------------------------------
#spag.plot(BC, plot.type="spag") #plot as stacked spaghetti 
BC_stats<-rwl.stats(BC) # first and last year, distribution parameters
#spag.plot(BM, plot.type="spag") #plot as stacked spaghetti 
BM_stats<-rwl.stats(BM)
#spag.plot(OP, plot.type="spag") #plot as stacked spaghetti 
OP_stats<-rwl.stats(OP)
#spag.plot(PM, plot.type="spag") #plot as stacked spaghetti 
PM_stats<-rwl.stats(PM)
#spag.plot(WM, plot.type="spag") #plot as stacked spaghetti 
WM_stats<-rwl.stats(WM)

#Changing years to pith to NA if greater than 10
BC_attr$years2pith[BC_attr$years2pith > 10] <-NA
BM_attr$years2pith[BM_attr$years2pith > 10] <-NA
PM_attr$years2pith[PM_attr$years2pith > 10] <-NA
OP_attr$years2pith[OP_attr$years2pith > 10] <-NA
WM_attr$years2pith[WM_attr$years2pith > 10] <-NA

#Combine series stats with pith info
BC_age <- BC_stats %>% full_join(BC_attr, by = c("series"))%>% mutate(r_year = first - years2pith,) %>% 
                        mutate(r_age = year + years2pith) %>% drop_na() #calc. recruit age/year, drop series w/o pith info
BM_age <- BM_stats %>% full_join(BM_attr, by = c("series")) %>% mutate(r_year = first - years2pith) %>% 
                        mutate(r_age = year + years2pith) %>% drop_na()
PM_age <- PM_stats %>% full_join(PM_attr, by = c("series")) %>% mutate(r_year = first - years2pith) %>% 
                          mutate(r_age = year + years2pith) %>% drop_na()
OP_age <- OP_stats %>% full_join(OP_attr, by = c("series")) %>% mutate(r_year = first - years2pith) %>% 
                          mutate(r_age = year + years2pith) %>% drop_na()
WM_age <- WM_stats %>% full_join(WM_attr, by = c("series")) %>% mutate(r_year = first - years2pith) %>% 
  mutate(r_age = year + years2pith) %>% drop_na()

# Plot Recruitment Age ---------------------------------------------------------

ageD<-ggplot(BC_age, aes(x=r_year))+
  geom_bar(fill = '#006400') +
  #geom_line(aes(y = "samp.depth")) + 
  labs(x='Recruitment Decade', y='Number of Trees', title = 'Recruitment age of red spruce at Blackwoods')+ 
  theme_bw()+scale_x_binned(show.limits = TRUE, limits = c(1810, 2020), n.breaks = 21)+ 
  theme(axis.text.x=element_text(angle=30,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

ageD
ggsave(file="./figures/BC_age_dist2.jpg", ageD, dpi=300, width=5, height=3, units='in')

ageD<-ggplot(BM_age, aes(x=r_year))+
  geom_bar(fill = '#006400') +
  #geom_line(aes(y = "samp.depth")) + 
  labs(x='Recruitment Decade', y='Number of Trees', title = 'Recruitment age of red spruce at Beech Mtn')+ 
  theme_bw() + scale_x_binned(show.limits = TRUE, limits = c(1810, 2020), n.breaks = 21)+ 
  theme(axis.text.x=element_text(angle=30,hjust=1, size = 6), 
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

ageD
ggsave(file="./figures/BM_age_dist2.jpg", ageD, dpi=300, width=5, height=3, units='in')

ageD<-ggplot(PM_age, aes(x=r_year))+
  geom_bar(fill = '#006400') +
  #geom_line(aes(y = "samp.depth")) + 
  labs(x='Recruitment Decade', y='Number of Trees', title = 'Recruitment age of red spruce at Pemetic Mtn')+ 
  theme_bw()+scale_x_binned(show.limits = TRUE, limits = c(1810, 2020), n.breaks = 21)+ 
  theme(axis.text.x=element_text(angle=30,hjust=1, size = 6),
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

ageD
ggsave(file="./figures/PM_age_dist2.jpg", ageD, dpi=300, width=5, height=3, units='in')

ageD<-ggplot(OP_age, aes(x=r_year))+
  geom_bar(fill = '#006400') +
  #geom_line(aes(y = "samp.depth")) + 
  labs(x='Recruitment Decade', y='Number of Trees', title = 'Recruitment age of red spruce at Otter Point')+ 
  theme_bw()+scale_x_binned(show.limits = TRUE, limits = c(1810, 2020), n.breaks = 21)+ 
  theme(axis.text.x=element_text(angle=30,hjust=1, size = 6),
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

ageD
ggsave(file="./figures/OP_age_dist2.jpg", ageD, dpi=300, width=5, height=3, units='in')

ageD<-ggplot(WM_age, aes(x=r_year))+
  geom_bar(fill = '#006400') +
  #geom_line(aes(y = "samp.depth")) + 
  labs(x='Recruitment Decade', y='Number of Trees', title = 'Recruitment age of red spruce at Bernard Mtn')+ 
  theme_bw()+scale_x_binned(show.limits = TRUE, limits = c(1740, 2020), n.breaks = 21)+ 
  theme(axis.text.x=element_text(angle=30,hjust=1, size = 6),
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

ageD
ggsave(file="./figures/WM_age_dist2.jpg", ageD, dpi=300, width=5, height=3, units='in')

# Detrending --------------------------------------------------------------
#i.detrend() #Detrend one by one (recommended) or detrend(df, method = c("Friedman") to do all at once
#BC_detrend<- i.detrend(BC) #to detrend each series individually
BC_Fr<-detrend(BC, method = c("Friedman")) #for climate growth relationship
BC_M<-detrend(BC, method = c("Mean")) #for stand dynamics

OP_Fr<-detrend(OP, method = c("Friedman")) #for climate growth relationship
OP_M<-detrend(OP, method = c("Mean")) #for stand dynamics

BM_Fr<-detrend(BM, method = c("Friedman")) #for climate growth relationship
BM_M<-detrend(BM, method = c("Mean")) #for stand dynamics

PM_Fr<-detrend(PM, method = c("Friedman")) #for climate growth relationship
PM_M<-detrend(PM, method = c("Mean")) #for stand dynamics

WM_Fr<-detrend(WM, method = c("Friedman")) #for climate growth relationship
WM_M<-detrend(WM, method = c("Mean")) #for stand dynamics
# Sample Signal Strength --------------------------------------------------
#Run subsample signal strength analysis to see where the chronology should be truncated due to to small sample depth.
#sss(df_rwi; df_ids) need to get a table with IDs (col 1 = tree ID, col 2 = core ID)

plot_sss<-function(d.df, df.rwl){
  site<-deparse(substitute(df.rwl))
  fig_name<-paste0("./chrono_graphs/", site, "_", "sss", '.jpg')
  df.sss<-sss(d.df)
  yr<-time(df.rwl)
  ppi<-300
  jpeg(file = fig_name, units='px', width=10*ppi, height=7*ppi, res=300) 
  plot(yr, df.sss,type="p", ylim=c(0.4,1), xaxt="none",
       col="blue", lwd=2, xlab="year", ylab="SSS") +
    axis(3, seq(1600, 2020, 20)) +
    abline(h=0.85, col="red") 
  dev.off()
} # function to make SSS generation faster
#need folder chrono_graphs in directory
#d.df = detrend df; df.rwl = imported rwl file 

plot_sss(BC_Fr, BC) #BC should truncate at 1880s

plot_sss(BM_Fr, BM) #BM should truncate at 1880s

plot_sss(OP_Fr, OP) #OP should truncate at 1960s

plot_sss(PM_Fr, PM) #PM should truncate at 1880s

plot_sss(WM_Fr, WM) #WM should truncate at 1785s
# Chronologies------------------------------------------------------------
#generate chronology file, format .crn
BC_Fr_chr <- chron(BC_Fr, prefix = "BCF")
BC_M_chr <- chron(BC_M, prefix = "BCM")

BM_Fr_chr <- chron(BM_Fr, prefix = "BMF")
BM_M_chr <- chron(BM_M, prefix = "BMM")

OP_Fr_chr <- chron(OP_Fr, prefix = "OPF")
OP_M_chr <- chron(OP_M, prefix = "OPM")

PM_Fr_chr <- chron(PM_Fr, prefix = "PMF")
PM_M_chr <- chron(PM_M, prefix = "PMM")

WM_Fr_chr <- chron(WM_Fr, prefix = "WMF")
WM_M_chr <- chron(WM_M, prefix = "WMM")

#plot chronology

#jpeg(file = "./chrono_graphs/BC_Fr_chr.jpg", units='px', width=10*ppi, height=7*ppi, res=300)
#plot.crn(BC_Fr_chr, add.spline=TRUE, nyrs=20) # add xlim = c(1880, 2018), ylim = c(0,2) to limit x to years specified by SSS or time period of interest
#dev.off()

print_crn<-function(df_chr){
  site<-deparse(substitute(df_chr))
  fig_name<-paste0("./chrono_graphs/", site, '.jpg')
  ppi<-300
  jpeg(file = fig_name, units='px', width=10*ppi, height=7*ppi, res=300) 
  df.plot<-plot.crn(df_chr, add.spline=TRUE, nyrs=20, lab = c(7,7,7))
  dev.off()
}# function to print a chronology to specified folder; lab used to specify number of axis tick marks

print_crn(BC_Fr_chr)
print_crn(BC_M_chr)

print_crn(BM_Fr_chr)
print_crn(BM_M_chr)

print_crn(OP_Fr_chr)
print_crn(OP_M_chr)

print_crn(PM_Fr_chr)
print_crn(PM_M_chr)

print_crn(WM_Fr_chr)
print_crn(WM_M_chr)
# Disturbance Reconstruction ----------------------------------------------
#use absolute increase method for determining growth release. PIRU threshold .58 (Shawn data)
#Sub par compared to Shawn's SAS code, how to ID gap recruits? Do I use RWI or raw rign widths?


