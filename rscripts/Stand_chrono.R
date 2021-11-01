#Davis project: signal sample depth, standardization, and chronologies

# Dendrochronology packages
#install.packages("dplR")
library(dplR)
#library(treeclim)
# Data managment packages
library(ggplot2)
library(gdata)
library(dplyr)
library(tidyr)
library(tibble)

# Load ring width data ----------------------------------------------------
BC <- read.tucson('../data/PIRU_best/BC_best.rwl')
BM <- read.tucson('../data/PIRU_best/BM_best.rwl')
PM <- read.tucson('../data/PIRU_best/PM_best.rwl')
OP <- read.tucson('../data/PIRU_best/OP_best.rwl')

# Summary Stats -----------------------------------------------------------
spag.plot(BC, plot.type="spag") #plot as stacked spaghetti 
BC_stats<-rwl.stats(BC)

spag.plot(BM, plot.type="spag") #plot as stacked spaghetti 
BM_stats<-rwl.stats(BM)

spag.plot(OP, plot.type="spag") #plot as stacked spaghetti 
OP_stats<-rwl.stats(OP)

spag.plot(PM, plot.type="spag") #plot as stacked spaghetti 
PM_stats<-rwl.stats(PM)

# Detrending --------------------------------------------------------------
#i.detrend() #Detrend one by one (recommended) or detrend(df, method = c("Friedman") to do all at once
BC_Fr<-detrend(BC, method = c("Friedman")) #for climate growth relationship
BC_M<-detrend(BC, method = c("Mean")) #for stand dynamics

OP_Fr<-detrend(OP, method = c("Friedman")) #for climate growth relationship
OP_M<-detrend(OP, method = c("Mean")) #for stand dynamics

BM_Fr<-detrend(BM, method = c("Friedman")) #for climate growth relationship
BM_M<-detrend(BM, method = c("Mean")) #for stand dynamics

PM_Fr<-detrend(PM, method = c("Friedman")) #for climate growth relationship
PM_M<-detrend(PM, method = c("Mean")) #for stand dynamics

# Chronologies ------------------------------------------------------------
BC_Fr_chr <- chron(BC_Fr, prefix = "BCF")
BC_M_chr <- chron(BC_M, prefix = "BCM")


