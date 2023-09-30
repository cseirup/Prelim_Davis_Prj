#Davis project: signal sample depth, standardization, and chronologies


#############################
#NEEDS TO BE UPDATED WITH CODE FROM Standardize_Chrono and AgeStructure_Disturbance_Chrono
##################################

# Dendrochronology packages
#install.packages("treeclim")
library(dplR)
library(TRADER)
library(treeclim)
# Data management packages
library(tidyverse)

# Load ring width data ----------------------------------------------------
BC <- read.tucson('../Davis_data/PIRU_best/BC_best.rwl')
BM <- read.tucson('../Davis_data/PIRU_best/BM_best.rwl')
PM <- read.tucson('../Davis_data/PIRU_best/PM_best.rwl')
OP <- read.tucson('../Davis_data/PIRU_best/OP_best.rwl')
WMA <- read.fh('../Davis_data/PIRU_best/WMALL_best_stripped.fh')
IB <- read.fh('../Davis_data/PIRU_best/IB_best_stripped.fh')
BH <- read.fh('../Davis_data/PIRU_best/BH_best_stripped.fh')
PM_TSCA <- read.fh('../Davis_data/PIRU_best/PM_TSCA_best_stripped.fh')

#Use years to pith from CDendro attribute output: Cores measured on Velmex do not have rings to pith data in attribute files
BC_attr <- read.table('../Davis_data/PIRU_best/BC_best_STRIPPED_attributes.txt', header = TRUE)#CDENDRO outputs NULL d2p/y2p as 0, change to NA + remove header manually before bringing into R.
BM_attr <- read.table('../Davis_data/PIRU_best/BM_best_STRIPPED_attributes.txt', header = TRUE)
PM_attr <- read.table('../Davis_data/PIRU_best/PM_best_STRIPPED_attributes.txt', header = TRUE)
OP_attr <- read.table('../Davis_data/PIRU_best/OP_STRIPPED_attributes.txt', header = TRUE)
WMA_attr <- read.table('../Davis_data/PIRU_best/WMALL_best_STRIPPED_attributes.txt', header = TRUE)
IB_attr <- read.table('../Davis_data/PIRU_best/IB_best_STRIPPED_attributes.txt', header = TRUE)
BH_attr <- read.table('../Davis_data/PIRU_best/BH_best_STRIPPED_attributes.txt', header = TRUE)
PM_TSCA_attr <- read.table('../Davis_data/PIRU_best/PM_TSCA_best_STRIPPED_attributes.txt', header = TRUE)

attr <- bind_rows(BC_attr, BM_attr, PM_attr, OP_attr, WMA_attr, IB_attr, BH_attr, PM_TSCA_attr)
attr <- rename(attr, Core_ID = series)

#adding missing rings to pith to attr file
a2 <- attr %>% left_join(Cores, by = "Core_ID") #left join with attr to only include cross-dated cores included in analysis
names(Cores2)
Cores3 <- Cores2 %>% mutate(diff = years2pith-Rings.to.pith) #correct any discrepancies
Cores4 <- Cores3 %>%  select(1:6)

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
if(!dir.exists("../chrono_graphs")){dir.create('../chrono_graphs')} # creating a folder for graphs

plot_sss(BC_Fr, BC) #BC should truncate at 1880s

plot_sss(BM_Fr, BM) #BM should truncate at 1880s

plot_sss(OP_Fr, OP) #OP should truncate at 1960s

plot_sss(PM_Fr, PM) #PM should truncate at 1880s

plot_sss(WM_Fr, WM) #WM should truncate at 1785s
# Chronologies------------------------------------------------------------
#generate chronology file, format .crn
BC_Fr_chr <- chron(BC_Fr, prefix = "BCF")
BC_M_chr <- chron(BC_M, prefix = "BCM")
BC_chr <- chron(BC, prefix = "BC")

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
print_crn(BC_chr)

print_crn(BM_Fr_chr)
print_crn(BM_M_chr)

print_crn(OP_Fr_chr)
print_crn(OP_M_chr)

print_crn(PM_Fr_chr)
print_crn(PM_M_chr)

print_crn(WM_Fr_chr)
print_crn(WM_M_chr)

# Climate ----------------------------------------------

# Load Climate Data -------------------------------------------------------
clim <-read.csv("../data/Prism.csv",header=TRUE)
spei <- read.csv2("../data/SPEI_BH.csv")
cld <- read.csv("../data/BH_clouds.csv", header = TRUE)

clim1<- separate(clim, Date, into = c("year", 'month'), sep = "[^[:alnum:]]+")
clim2<-rename(clim1, ppt = ppt..mm., tmin = tmin..degrees.C., tmean = tmean..degrees.C., tmax = tmax..degrees.C., 
              vpdmin = vpdmin..hPa., vpdmax = vpdmax..hPa.)
clim2$year<-as.numeric(clim2$year)
clim2$month<-as.numeric(clim2$month) # clim2 starts at 1895, no spei
clim.tmax <-select(clim2, year, month, tmax)
clim.tmin <-select(clim2, year, month, tmin)
clim.tmean <-select(clim2, year, month, tmean)
clim.ppt <-select(clim2, year, month, ppt)
clim.vpdmin <- select(clim2, year, month, vpdmin)
clim.vpdmax <- select(clim2, year, month, vpdmax)

colnames(spei)
spei$days.since.1900.1.1<-as.character(spei$days.since.1900.1.1)
spei1 <- separate(spei,days.since.1900.1.1, into = c("year","month"), sep = "[^[:alnum:]]+")
spei1$year<-as.numeric(spei1$year)
spei1$month<-as.numeric(spei1$month)
spei1$spei<-as.numeric(spei1$spei)
clim.spei<-spei1
str(clim.spei)

cld1 <- cld[,c(1:13)]
cld2 <- rename(cld1, "year" = Year, '1' = Jan, '2' = Feb, '3' = Mar, '4' = Apr, '5' = May, '6' = Jun, '7' = Jul, '8' = Aug, 
               '9' = Sep, '10' = Oct, '11' = Nov, '12' = Dec)
clim.cld <- gather(cld2, '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',  
                   key = "month", value = "cld_cov") # cloud cover by month,1948-2018
clim.cld$month<-as.numeric(clim.cld$month)

clim3 <- clim2 %>% full_join(clim.spei, by = c("year", "month"))
clim_all <- clim3 %>% full_join(clim.cld, by = c("year", "month"))

# Truncate chronos to .85 sss--------------------------------------------
BM_Fr_chr_tr <- BM_Fr_chr %>% rownames_to_column(var = "year") %>% mutate(sss = sss(BM_Fr)) %>% filter(sss>=.85) %>% 
                               column_to_rownames(var = "year")
BC_Fr_chr_tr <- BC_Fr_chr %>% rownames_to_column(var = "year") %>% mutate(sss = sss(BC_Fr)) %>% filter(sss>=.85) %>%   
                               column_to_rownames(var = "year")
OP_Fr_chr_tr <- OP_Fr_chr %>% rownames_to_column(var = "year") %>% mutate(sss = sss(OP_Fr)) %>% filter(sss>=.85) %>% 
                               column_to_rownames(var = "year")
PM_Fr_chr_tr <- PM_Fr_chr %>% rownames_to_column(var = "year") %>% mutate(sss = sss(PM_Fr)) %>% filter(sss>=.85) %>% 
                               column_to_rownames(var = "year")
WM_Fr_chr_tr <- WM_Fr_chr %>% rownames_to_column(var = "year") %>% mutate(sss = sss(WM_Fr)) %>% filter(sss>=.85) %>% 
                               column_to_rownames(var = "year")

# Response function analysis ----------------------------------------------
if(!dir.exists("../climate_graphs")){dir.create('../climate_graphs')} # creating a data folder for datasets when created
Cdata_list = list(clim.cld, clim.spei, clim.ppt, clim.tmax, clim.tmean, clim.tmin, clim.vpdmax, clim.vpdmin)

# simplifying function for purrr; going to try all climate for one site at a time
#BC - timespan limited by climate
fBCresp<-function(climdf){ #climdf = list of climate variable datasets
  title <- paste0("Blackwoods", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df<- dcc(chrono = BC_Fr_chr_tr, climate = climdf, selection = -6:9, method='response', 
              dynamic = "static", win_size = 25, win_offset = 1,
              start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
              boot = "stationary", sb = TRUE)
  
  df_2 <- data.frame(df$coef$month,df$coef$coef,df$coef$significant)
  colnames(df_2) <- c('month','coef','significant')
  df_2$mon=1:16
  
  plot.df2<-ggplot(data=df_2, aes(x=mon,y=coef,fill=significant))+
    geom_bar(stat='identity')+theme_classic(base_size=12)+
    ylab('Coefficients\n')+ggtitle(title)+
    scale_fill_manual(values=c("lightsteelblue","red3"))+ theme(axis.title.x=element_blank())+
    scale_x_continuous(breaks=c(1:16),labels=df_2$month)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot.df2
}
fBCresp(clim.spei)
BC_climplots <- map(Cdata_list, fBCresp)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/BC_climplots.pdf", height = 7, width = 10)
BC_climplots
dev.off()

#BM - timespan limited by climate
fBMresp<-function(climdf){ #clim = list of climate variable datasets
  #para <- deparse(substitute(climdf))
  title <- paste0("Beech Mtn", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df<- dcc(chrono = BM_Fr_chr_tr, climate = climdf, selection = -6:9, method='response', 
           dynamic = "static", win_size = 25, win_offset = 1,
           start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
           boot = "stationary", sb = TRUE)
  
  df_2 <- data.frame(df$coef$month,df$coef$coef,df$coef$significant)
  colnames(df_2) <- c('month','coef','significant')
  df_2$mon=1:16
  
  plot.df2<-ggplot(data=df_2, aes(x=mon,y=coef,fill=significant))+
    geom_bar(stat='identity')+theme_classic(base_size=12)+
    ylab('Coefficients\n')+ggtitle(title)+
    scale_fill_manual(values=c("lightsteelblue","red3"))+ theme(axis.title.x=element_blank())+
    scale_x_continuous(breaks=c(1:16),labels=df_2$month)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot.df2
}
fBMresp(clim.spei)
BM_climplots <- map(Cdata_list, fBMresp)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/BM_climplots.pdf", height = 7, width = 10)
BM_climplots
dev.off()

#OP - timespan limited by chrono, title manually entered
fOPresp<-function(climdf){ #clim = list of climate variable datasets
  #para <- deparse(substitute(climdf))
  title <- paste0("Otter Point",  " ", names(climdf[3])," ", "1949", "-", "2018")
  
  df<- dcc(chrono = OP_Fr_chr_tr, climate = climdf, selection = -6:9, method='response', 
           dynamic = "static", win_size = 25, win_offset = 1,
           start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
           boot = "stationary", sb = TRUE)
  
  df_2 <- data.frame(df$coef$month,df$coef$coef,df$coef$significant)
  colnames(df_2) <- c('month','coef','significant')
  df_2$mon=1:16
  
  plot.df2<-ggplot(data=df_2, aes(x=mon,y=coef,fill=significant))+
    geom_bar(stat='identity')+theme_classic(base_size=12)+
    ylab('Coefficients\n')+ggtitle(title)+
    scale_fill_manual(values=c("lightsteelblue","red3"))+ theme(axis.title.x=element_blank())+
    scale_x_continuous(breaks=c(1:16),labels=df_2$month)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot.df2
}
fOPresp(clim.spei)
OP_climplots <- map(Cdata_list, fOPresp)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/OP_climplots.pdf", height = 7, width = 10)
OP_climplots
dev.off()

#PM
fPMresp<-function(climdf){ #clim = list of climate variable datasets
  #para <- deparse(substitute(climdf))
  title <- paste0("Pemetic Mtn", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df<- dcc(chrono = PM_Fr_chr_tr, climate = climdf, selection = -6:9, method='response', 
           dynamic = "static", win_size = 25, win_offset = 1,
           start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
           boot = "stationary", sb = TRUE)
  
  df_2 <- data.frame(df$coef$month,df$coef$coef,df$coef$significant)
  colnames(df_2) <- c('month','coef','significant')
  df_2$mon=1:16
  
  plot.df2<-ggplot(data=df_2, aes(x=mon,y=coef,fill=significant))+
    geom_bar(stat='identity')+theme_classic(base_size=12)+
    ylab('Coefficients\n')+ggtitle(title)+
    scale_fill_manual(values=c("lightsteelblue","red3"))+ theme(axis.title.x=element_blank())+
    scale_x_continuous(breaks=c(1:16),labels=df_2$month)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot.df2
}
fPMresp(clim.spei)
PM_climplots <- map(Cdata_list, fPMresp)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/PM_climplots.pdf", height = 7, width = 10)
PM_climplots
dev.off()

#WM
fWMresp<-function(climdf){ #clim = list of climate variable datasets
  #para <- deparse(substitute(climdf))
  title <- paste0("Bernard Mtn", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df<- dcc(chrono = WM_Fr_chr_tr, climate = climdf, selection = -6:9, method='response', 
           dynamic = "static", win_size = 25, win_offset = 1,
           start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
           boot = "stationary", sb = TRUE)
  
  df_2 <- data.frame(df$coef$month,df$coef$coef,df$coef$significant)
  colnames(df_2) <- c('month','coef','significant')
  df_2$mon=1:16
  
  plot.df2<-ggplot(data=df_2, aes(x=mon,y=coef,fill=significant))+
    geom_bar(stat='identity')+theme_classic(base_size=12)+
    ylab('Coefficients\n')+ggtitle(title)+
    scale_fill_manual(values=c("lightsteelblue","red3"))+ theme(axis.title.x=element_blank())+
    scale_x_continuous(breaks=c(1:16),labels=df_2$month)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot.df2
}
fWMresp(clim.spei)
WM_climplots <- map(Cdata_list, fWMresp)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/WM_climplots.pdf", height = 7, width = 10)
WM_climplots
dev.off()

# Response graphs points --------------------------------------------------
#BC
fBCresp_pt<-function(climdf){ #climdf = list of climate variable datasets
  title <- paste0("Blackwoods", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))

df <- dcc(BC_Fr_chr_tr, climdf, selection = -6:9, method = 'response', dynamic = "static", win_size = 25, win_offset = 1,
          start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
          boot = "stationary", sb = TRUE)
df_2 <- data.frame(df$coef$month, df$coef$coef, df$coef$significant, df$coef$ci_lower, df$coef$ci_upper)
colnames(df_2) <- c('month', 'coef', 'significant', 'ci_lower', 'ci_upper')
df_2$mon <- 1:16
df_2$significant <- as.factor(df_2$significant)

plot_pts <- ggplot(df_2, aes(x = mon, y = coef))+
            geom_point(aes(colour = significant))+
            scale_color_manual(values = c("grey", '#CD5C5C'))+
            theme_classic()+ ylab('Coefficients\n')+ ggtitle(title)+
            geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, colour = significant))+
            scale_x_continuous(breaks=c(1:16), labels = df_2$month)+
            geom_hline(yintercept = 0, lty = 2, lwd =1, color = "darkgrey")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))  
plot_pts

}
fBCresp_pt(clim.spei)
BC_climplot_point <- map(Cdata_list, fBCresp_pt)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/BC_climplot_point.pdf", height = 7, width = 10)
BC_climplot_point
dev.off()

#BM
fBMresp_pt<-function(climdf){ #climdf = list of climate variable datasets
  title <- paste0("Beech Mtn", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df <- dcc(BM_Fr_chr_tr, climdf, selection = -6:9, method = 'response', dynamic = "static", win_size = 25, win_offset = 1,
            start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
            boot = "stationary", sb = TRUE)
  df_2 <- data.frame(df$coef$month, df$coef$coef, df$coef$significant, df$coef$ci_lower, df$coef$ci_upper)
  colnames(df_2) <- c('month', 'coef', 'significant', 'ci_lower', 'ci_upper')
  df_2$mon <- 1:16
  df_2$significant <- as.factor(df_2$significant)
  
  plot_pts <- ggplot(df_2, aes(x = mon, y = coef))+
    geom_point(aes(colour = significant))+
    scale_color_manual(values = c("grey", '#CD5C5C'))+
    theme_classic()+ ylab('Coefficients\n')+ ggtitle(title)+
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, colour = significant))+
    scale_x_continuous(breaks=c(1:16), labels = df_2$month)+
    geom_hline(yintercept = 0, lty = 2, lwd =1, color = "darkgrey")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  plot_pts
  
}
fBMresp_pt(clim.spei)
BM_climplot_point <- map(Cdata_list, fBMresp_pt)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/BM_climplot_point.pdf", height = 7, width = 10)
BM_climplot_point
dev.off()

#OP
fOPresp_pt<-function(climdf){ #climdf = list of climate variable datasets
  title <- paste0("Otter Point", " ", names(climdf[3])," ", "1949", "-", "2018")
  
  df <- dcc(OP_Fr_chr_tr, climdf, selection = -6:9, method = 'response', dynamic = "static", win_size = 25, win_offset = 1,
            start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
            boot = "stationary", sb = TRUE)
  df_2 <- data.frame(df$coef$month, df$coef$coef, df$coef$significant, df$coef$ci_lower, df$coef$ci_upper)
  colnames(df_2) <- c('month', 'coef', 'significant', 'ci_lower', 'ci_upper')
  df_2$mon <- 1:16
  df_2$significant <- as.factor(df_2$significant)
  
  plot_pts <- ggplot(df_2, aes(x = mon, y = coef))+
    geom_point(aes(colour = significant))+
    scale_color_manual(values = c("grey", '#CD5C5C'))+
    theme_classic()+ ylab('Coefficients\n')+ ggtitle(title)+
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, colour = significant))+
    scale_x_continuous(breaks=c(1:16), labels = df_2$month)+
    geom_hline(yintercept = 0, lty = 2, lwd =1, color = "darkgrey")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  plot_pts
  
}
fOPresp_pt(clim.spei)
OP_climplot_point <- map(Cdata_list, fOPresp_pt)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/OP_climplot_point.pdf", height = 7, width = 10)
OP_climplot_point
dev.off()

#PM
fPMresp_pt<-function(climdf){ #climdf = list of climate variable datasets
  title <- paste0("Pemetic Mtn", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df <- dcc(PM_Fr_chr_tr, climdf, selection = -6:9, method = 'response', dynamic = "static", win_size = 25, win_offset = 1,
            start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
            boot = "stationary", sb = TRUE)
  df_2 <- data.frame(df$coef$month, df$coef$coef, df$coef$significant, df$coef$ci_lower, df$coef$ci_upper)
  colnames(df_2) <- c('month', 'coef', 'significant', 'ci_lower', 'ci_upper')
  df_2$mon <- 1:16
  df_2$significant <- as.factor(df_2$significant)
  
  plot_pts <- ggplot(df_2, aes(x = mon, y = coef))+
    geom_point(aes(colour = significant))+
    scale_color_manual(values = c("grey", '#CD5C5C'))+
    theme_classic()+ ylab('Coefficients\n')+ ggtitle(title)+
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, colour = significant))+
    scale_x_continuous(breaks=c(1:16), labels = df_2$month)+
    geom_hline(yintercept = 0, lty = 2, lwd =1, color = "darkgrey")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  plot_pts
  
}
fPMresp_pt(clim.spei)
PM_climplot_point <- map(Cdata_list, fPMresp_pt)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/PM_climplot_point.pdf", height = 7, width = 10)
PM_climplot_point
dev.off()

BC_pt_ppt <- fBCresp_pt(clim.ppt)
jpeg(BC_pt_ppt,"./climate_graphs/indiv/BC_pt_ppt.jpeg")

#WM
fWMresp_pt<-function(climdf){ #climdf = list of climate variable datasets
  title <- paste0("Western Mtn", " ", names(climdf[3])," ", min(climdf[1]), "-", max(climdf[1]))
  
  df <- dcc(WM_Fr_chr_tr, climdf, selection = -6:9, method = 'response', dynamic = "static", win_size = 25, win_offset = 1,
            start_last = TRUE, timespan = NULL, var_names = NULL, ci = 0.05,
            boot = "stationary", sb = TRUE)
  df_2 <- data.frame(df$coef$month, df$coef$coef, df$coef$significant, df$coef$ci_lower, df$coef$ci_upper)
  colnames(df_2) <- c('month', 'coef', 'significant', 'ci_lower', 'ci_upper')
  df_2$mon <- 1:16
  df_2$significant <- as.factor(df_2$significant)
  
  plot_pts <- ggplot(df_2, aes(x = mon, y = coef))+
    geom_point(aes(colour = significant))+
    scale_color_manual(values = c("grey", '#CD5C5C'))+
    theme_classic()+ ylab('Coefficients\n')+ ggtitle(title)+
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, colour = significant))+
    scale_x_continuous(breaks=c(1:16), labels = df_2$month)+
    geom_hline(yintercept = 0, lty = 2, lwd =1, color = "darkgrey")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  plot_pts
  
}
fWMresp_pt(clim.spei)
WM_climplot_point <- map(Cdata_list, fWMresp_pt)
#Save all plots to pdf for easy viewing
pdf("./climate_graphs/WM_climplot_point.pdf", height = 7, width = 10)
WM_climplot_point
dev.off()

# printing individual figs
BC_pt_ppt <- fBCresp_pt(clim.ppt)
BC_pt_ppt
ggsave(file = "./climate_graphs/indiv/BC_pt_ppt.jpeg", BC_pt_ppt, dpi=300, width=5, height=3, units='in')
dev.off()

BM_pt_ppt <- fBMresp_pt(clim.ppt)
BM_pt_ppt
ggsave(file = "./climate_graphs/indiv/BM_pt_ppt.jpeg", BM_pt_ppt, dpi=300, width=5, height=3, units='in')
dev.off()

OP_pt_ppt <- fOPresp_pt(clim.ppt)
OP_pt_ppt
ggsave(file = "./climate_graphs/indiv/OP_pt_ppt.jpeg", OP_pt_ppt, dpi=300, width=5, height=3, units='in')
dev.off()

PM_pt_ppt <- fPMresp_pt(clim.ppt)
PM_pt_ppt
ggsave(file = "./climate_graphs/indiv/PM_pt_ppt.jpeg", PM_pt_ppt, dpi=300, width=5, height=3, units='in')
dev.off()

WM_pt_ppt <- fWMresp_pt(clim.ppt)
WM_pt_ppt
ggsave(file = "./climate_graphs/indiv/WM_pt_ppt.jpeg", WM_pt_ppt, dpi=300, width=5, height=3, units='in')
dev.off()

# Variable growth response ------------------------------------------------

#### Note: Win_size works at 30, but errors if < 30
#### Note: .mean(3:5) = mean values months 3-5, etc 

#BC - first moving analysis, then evolving
vBC_ppt_mv <- dcc(chrono=BC_Fr_chr_tr, climate=clim.ppt,
          selection = -6:9,
          method = 'response', dynamic = "moving", win_size = 30, 
          sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vBC_ppt_mv_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vBC_ppt_mv_plot <- plot(vBC_ppt_mv)
vBC_ppt_mv_plot
dev.off()

vBC_ppt_ev <- dcc(chrono=BC_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "evolving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vBC_ppt_ev_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vBC_ppt_ev_plot <- plot(vBC_ppt_ev)
vBC_ppt_ev_plot
dev.off()

#BM - first moving analysis, then evolving
vBM_ppt_mv <- dcc(chrono=BM_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "moving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vBM_ppt_mv_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vBM_ppt_mv_plot <- plot(vBM_ppt_mv)
vBM_ppt_mv_plot
dev.off()

vBM_ppt_ev <- dcc(chrono=BM_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "evolving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vBM_ppt_ev_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vBM_ppt_ev_plot <- plot(vBM_ppt_ev)
vBM_ppt_ev_plot
dev.off()

#OP - first moving analysis, then evolving
vOP_ppt_mv <- dcc(chrono=OP_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "moving", win_size = 20, 
                  sb=FALSE, win_offset = 10)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vOP_ppt_mv_plot2.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vOP_ppt_mv_plot <- plot(vOP_ppt_mv)
vOP_ppt_mv_plot
dev.off()

vOP_ppt_ev <- dcc(chrono=OP_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "evolving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vOP_ppt_ev_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vOP_ppt_ev_plot <- plot(vOP_ppt_ev)
vOP_ppt_ev_plot
dev.off()

#PM - first moving analysis, then evolving
vPM_ppt_mv <- dcc(chrono=PM_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "moving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vPM_ppt_mv_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vPM_ppt_mv_plot <- plot(vPM_ppt_mv)
vPM_ppt_mv_plot
dev.off()

vPM_ppt_ev <- dcc(chrono=PM_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "evolving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vPM_ppt_ev_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vPM_ppt_ev_plot <- plot(vPM_ppt_ev)
vPM_ppt_ev_plot
dev.off()

#WM - first moving analysis, then evolving
vWM_ppt_mv <- dcc(chrono=WM_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "moving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vWM_ppt_mv_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vWM_ppt_mv_plot <- plot(vWM_ppt_mv)
vWM_ppt_mv_plot
dev.off()

vWM_ppt_ev <- dcc(chrono=WM_Fr_chr_tr, climate=clim.ppt,
                  selection = -6:9,
                  method = 'response', dynamic = "evolving", win_size = 30, 
                  sb=FALSE, win_offset = 15)
ppi<-300
jpeg(file = "./climate_graphs/indiv/vWM_ppt_ev_plot.jpg", units='px', width=10*ppi, height=7*ppi, res=300) 
vWM_ppt_ev_plot <- plot(vWM_ppt_ev)
vWM_ppt_ev_plot
dev.off()