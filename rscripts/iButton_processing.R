#Code written by Jay Wason to process first batch of ibutton data for SCS Davis project. Modified by Camilla Seirup 10/24/22 to combine second year of data
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

library(tidyverse)

# Processing 2020-2021 data -----------------------------------------------


setwd("C:/01_NETN/Forest_Health/R_Dev/Davis_data/ibuttons/Acadia_ibuttons_06_23_21/raw")
ibutton_list = c(list.files())

dat.out = data.frame()
for(i in ibutton_list){
  #i = "bernard_original_north_soilT.csv"
  #if thermochron (soil) then change number of rows to skip at read in
  my.skip = 19
  if(substrRight(i, 8) == "oilT.csv"){my.skip = 14}
  #read in the ibutton data
  temp = read.csv(file = i, header = F, skip = 0)#pull out serial number
  serial = as.character(substrRight(as.character(temp$V1[2]), 16))
  dat = read.csv(file = i, header = T, skip = my.skip)
  dat$file = i
  dat$serial = as.character(serial)
dat.out = rbind(dat.out, dat)
}

setwd("C:/01_NETN/Forest_Health/R_Dev/Davis_data/ibuttons")
str(dat.out)


#I think I want three columns for each ibutton location. AirT, AirH, and SoilT
#so, merge in the ibutton data first.
ibutton_info = read.csv("iButton_info_2021_06_22.csv", header = T)
ibutton_info2 = ibutton_info[c("iButton_name", "serial", "site", "site_code", "NSEW", "air_or_soil")]
dat.out2 = right_join(ibutton_info2, dat.out, by = "serial")
unique(substrRight(as.character(dat.out2$Date.Time), 8)) #substring to view minutes and seconds for measurements
## all times are either on the hour or one minute past
#so we can just delete the minutes and seconds and that would round to the nearest even hour
#testing converting this
as.POSIXct(strptime("7/29/20 2:00:01 PM", "%m/%d/%y %I:%M:%S %p"))
dat.out2$Date.Time2 = as.POSIXct(strptime(as.character(dat.out2$Date.Time), "%m/%d/%y %I:%M:%S %p"))
str(dat.out2)

dat.out2$Date.TimeR = as.POSIXct(round(dat.out2$Date.Time2, units = "hours"))
#okay, should be able to use this to merge things together
temp = subset(dat.out2, site_code == "BEMO")
temp2 = subset(temp, NSEW == "North")
temp3 = subset(temp2, air_or_soil == "soil")

unique(dat.out2$site_code)
unique(ibutton_info2$site_code)


dat.combined = data.frame()
for(i in unique(dat.out2$site_code)){
  #i = "BEMO"
  temp = subset(dat.out2, site_code == i)
  
  for(x in unique(levels(as.factor(as.character(temp$NSEW))))){
    #x = "North"
    temp2 = subset(temp, NSEW == x)

    temp.air = subset(temp2, air_or_soil == "air")
    temp.air.t = subset(temp.air, Unit == "C") #that is our base data
    #trim out some garbage columns
    temp.air.t$airT = temp.air.t$Value
    temp.air.t = temp.air.t[c("site", "site_code", "NSEW", "Date.TimeR", "airT")]
    
    #isolate humidity
    temp.air.h = subset(temp.air, Unit == "%RH")[c("Date.TimeR", "Value")]
    temp.air.h$airH = temp.air.h$Value
    temp.air.h = temp.air.h[c("Date.TimeR", "airH")]
    temp.air.out = merge(temp.air.t, temp.air.h, by = "Date.TimeR")
    
    #isolate soil temp
    temp.soil.t = subset(temp2, air_or_soil == "soil")
    temp.soil.t$soilT = temp.soil.t$Value
    temp.soil.t = temp.soil.t[c("Date.TimeR", "soilT")]
    temp.out = merge(temp.air.out, temp.soil.t, by = "Date.TimeR", all = T)
    
    dat.combined = rbind(dat.combined, temp.out)
  }#NSEW loop
}#site loop end

str(dat.combined)

#these were deployed July 30, 2020.
#lets look at a few plots to confirm this.
temp  = subset(dat.combined, site_code == "BLWO")
temp2 = subset(temp, NSEW = "North")
temp3 = subset(temp2, Date.TimeR <= "2020-08-03 16:00:00")
plot(airT ~ Date.TimeR, data = temp3)
abline(v = as.POSIXct("2020-07-30 22:00:00"))
#yep. so lets trim so August 1st is first day.
dat.combined2 = subset(dat.combined, Date.TimeR >= "2020-08-01 00:00:00")

#now lets do some quality control
hist(dat.combined2$airT)
hist(dat.combined2$airH)
#will need to decide how to handle values above 100
hist(dat.combined2$soilT)
#looks like some warm soilT days that should be checked too.

#lets do this by site
dat.combined2$site_NSEW = as.factor(paste(dat.combined2$site_code, dat.combined2$NSEW, sep = "_"))
str(dat.combined2)


table(dat.combined2$site_code, dat.combined2$NSEW)
plot(dat.combined2$soilT ~ dat.combined2$site_NSEW)




#################################################################
#lets try to correct for saturation drift.

dat.out = data.frame()

for(i in unique(levels(dat.combined2$site_NSEW))){
  #https://datasheets.maximintegrated.com/en/ds/DS1923.pdf
  #i = "BECL_North"
  temp = subset(dat.combined2, site_NSEW == i)
  temp$ov70 = ifelse(temp$airH >= 70, 2, 0) #hours over 70 since last reading
  temp$ov70na = ifelse(temp$airH >= 70, 2, NA) #mark out the times not over 70 with NA
  #calculate the cummulative hours of 70 for each event.
  temp$ov70cum = with(temp, ave(ov70,cumsum(is.na(ov70na)),FUN=cumsum))
  
  #calculate the partial correction factor
  temp$hcor.par = (0.0156 * temp$airH * 2.54^(-0.3502 * temp$ov70cum)) / (1 + (temp$airT - 25)/100)
  #make it 0 if not during a period over 70
  temp$hcor.par2 = ifelse(temp$ov70cum > 0, temp$hcor.par, 0)
  
  #now sum those up for each event.
  temp$hcor.cum = with(temp, ave(hcor.par2,cumsum(is.na(ov70na)),FUN=cumsum))
  temp$airHc = temp$airH - temp$hcor.cum
  
  #still a few over 100 so then move all those to just 100%.
  temp$airHc[temp$airHc >= 100] <- 100 
  hist(temp$airHc)
  
  #remove old columns
  temp2 = subset (temp, select = -c(ov70, ov70na, ov70cum, hcor.par, hcor.par2, hcor.cum))
  dat.out = rbind(dat.out, temp2)   
  
}
########################################################

dat.out$soilT = as.numeric(dat.out$soilT)

#now we need to check data for each site

for(i in unique(levels(dat.out$site_NSEW))){
  # i = "BEMO_North"
  png(filename = paste("2020_2021_iButtons", paste(i,".png", sep = ""), sep = "_"), width = 9, height = 6, units = "in", res = 300)
  temp = subset(dat.out, site_NSEW == i)
  
  par(oma=c(3,3,3,1), mar=c(2,3,1,0.5), mgp = c(3, 1, 0)  )
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow  = T))
  layout.show(n=6)
  #air temp examination
  hist(temp$airT, main = "")
  mtext(side = 1, line = 2, "Air Temp (C)")
  
  plot(airT ~ Date.TimeR, data = temp, type = "l")
  abline(h = 0)
  mtext(side = 2, line = 2, "Air Temp (C)")
  

  #humidity examination
  hist(temp$airHc, main = "")
  mtext(side = 1, line = 2, "Relative Humidity (%)")
  
  plot(airHc ~ Date.TimeR, data = temp, type = "l")
  mtext(side = 2, line = 2, "Relative Humidity (%)")
  
  #soil temp
  if(length(na.omit(temp$soilT)) >= 2){
  hist(temp$soilT, main = "")
  mtext(side = 1, line = 2.5, "Soil Temp (C)")
  
  temp.soil = temp[!is.na(temp$soilT), ]
  plot(soilT ~ Date.TimeR, data = temp.soil, type = "l")
  abline(h = 0)
  mtext(side = 2, line = 2, "Soil Temp (C)")
    mtext(side = 1, line = 2.5, "Date")

  }
  mtext(side = 3, line = 1, i, outer = T)
  mtext(side = 2, line = 1, "Frequency", outer = T)
  
  
  dev.off()
}

###############################
#now calculate and output the daily and monthly datasets
dat.daily = data.frame()
dat.monthly = data.frame()

for(i in unique(levels(dat.out$site_NSEW))){
  # i = "BEMO_North"
  temp = subset(dat.out, site_NSEW == i)
  temp$day = as.POSIXct(round(as.POSIXct(temp$Date.TimeR), "day"))
  
  #air temperature
    #min
    dailymin = tapply(temp$airT, temp$day, min) #minimum temperature each day
    dailyrecord = data.frame(day=names(dailymin), airTmin=dailymin)
    #max#
    dailymax = tapply(temp$airT, temp$day, max)
    dailyrecord$airTmax = dailymax
    #average#
    dailyave = tapply(temp$airT, temp$day, mean)
    dailyrecord$airTave = dailyave
    
  #air humidity
    #min
    dailymin = tapply(temp$airHc, temp$day, min) #minimum humidity each day
    dailyrecord$airHmin = dailymin
    #max#
    dailymax = tapply(temp$airHc, temp$day, max)
    dailyrecord$airHmax = dailymax
    #average#
    dailyave = tapply(temp$airHc, temp$day, mean)
    dailyrecord$airHave = dailyave
    
  #VPD calculations
    temp$satvap = exp((16.78*temp$airT - 116.9) / (temp$airT + 237.3))
    temp$actvap = temp$satvap * (temp$airHc/100)
    temp$vpd = temp$satvap - temp$actvap#units are kilopascals
    temp$vpd = ifelse(temp$airT > 0, temp$vpd, NA)
    
    #min
    dailymin = tapply(temp$vpd, temp$day, min, na.rm = T) #minimum humidity each day
    dailyrecord$airVPDmin = dailymin    
    dailyrecord$airVPDmin = ifelse(is.infinite(dailyrecord$airVPDmin), NA, dailyrecord$airVPDmin)

    #max#
    dailymax = tapply(temp$vpd, temp$day, max, na.rm =T)
    dailyrecord$airVPDmax = dailymax
    dailyrecord$airVPDmax = ifelse(is.infinite(dailyrecord$airVPDmax), NA, dailyrecord$airVPDmax)
    
    #average#
    dailyave = tapply(temp$vpd, temp$day, mean, na.rm = T)
    dailyrecord$airVPDave = dailyave
    dailyrecord$airVPDave = ifelse(is.infinite(dailyrecord$airVPDave), NA, dailyrecord$airVPDave)
    
  #soil temperature
    #min
    dailymin = tapply(temp$soilT, temp$day, min, na.rm = T) #minimum temperature each day
    dailyrecord$soilTmin = dailymin
    #max#
    dailymax = tapply(temp$soilT, temp$day, max, na.rm = T)
    dailyrecord$soilTmax = dailymax
    #average#
    dailyave = tapply(temp$soilT, temp$day, mean,na.rm = T)
    dailyrecord$soilTave = dailyave
  
  dailyrecord$site_NSEW = i
  
  ##############
  #now for monthly calculations
  dailyrecord$month = strftime(dailyrecord$day, format = "%Y-%m")
  monthlyairTave = tapply(dailyrecord$airTave, dailyrecord$month, mean)
  monthlyrecord = data.frame(month=names(monthlyairTave), airTave = monthlyairTave)
  monthlyrecord$airTmin = tapply(dailyrecord$airTmin, dailyrecord$month, mean)
  monthlyrecord$airTmax = tapply(dailyrecord$airTmax, dailyrecord$month, mean)
  
  monthlyrecord$airHave = tapply(dailyrecord$airHave, dailyrecord$month, mean)
  monthlyrecord$airHmin = tapply(dailyrecord$airHmin, dailyrecord$month, mean)
  monthlyrecord$airHmax = tapply(dailyrecord$airHmax, dailyrecord$month, mean)
  
  monthlyrecord$airVPDave = tapply(dailyrecord$airVPDave, dailyrecord$month, mean)
  monthlyrecord$airVPDmin = tapply(dailyrecord$airVPDmin, dailyrecord$month, mean)
  monthlyrecord$airVPDmax = tapply(dailyrecord$airVPDmax, dailyrecord$month, mean)
  
  monthlyrecord$soilTave = tapply(dailyrecord$soilTave, dailyrecord$month, mean)
  monthlyrecord$soilTmin = tapply(dailyrecord$soilTmin, dailyrecord$month, mean)
  monthlyrecord$soilTmax = tapply(dailyrecord$soilTmax, dailyrecord$month, mean)
  
  monthlyrecord$length = tapply(dailyrecord$day, dailyrecord$month, function(x) length(unique(na.omit(x))))
  
  monthlyrecord$site_NSEW = i
  
  
  #output
  dat.daily = rbind(dat.daily, dailyrecord)
  dat.monthly = rbind(dat.monthly, monthlyrecord)
}

dat.daily20 <- dat.daily
dat.monthly20 <- dat.monthly

today = Sys.Date()
write.csv(dat.daily, file = paste(paste("2020_2021_dat_daily", today, sep = "_"), ".csv", sep = ""),
          row.names = F)
write.csv(dat.monthly, file = paste(paste("2020_2021_dat_monthly", today, sep = "_"), ".csv", sep = ""),
          row.names = F)


# Processing 2021-2022 data -----------------------------------------------

setwd("C:/01_NETN/Forest_Health/R_Dev/Davis_data/ibuttons/Acadia_ibuttons_2022_08_18/raw")
ibutton_list2 = c(list.files())

dat.out = data.frame()
for(i in ibutton_list2){
  #i = "bernard_original_north_soilT.csv"
  #if thermochron (soil) then change number of rows to skip at read in
  my.skip = 19
  if(substrRight(i, 6) == "ST.csv"){my.skip = 14}
  #read in the ibutton data
  temp = read.csv(file = i, header = F, skip = 0)#pull out serial number
  serial = as.character(substrRight(as.character(temp$V1[2]), 16))
  dat = read.csv(file = i, header = T, skip = my.skip)
  dat$file = i
  dat$serial = as.character(serial)
  dat.out = rbind(dat.out, dat)
}

setwd("C:/01_NETN/Forest_Health/R_Dev/Davis_data/ibuttons")
str(dat.out)


#I think I want three columns for each ibutton location. AirT, AirH, and SoilT
#so, merge in the ibutton data first.
ibutton_info = read.csv("iButton_info_2023_10_24.csv", header = T)
ibutton_info2 = ibutton_info[c("iButton_name", "serial", "site", "site_code", "NSEW", "air_or_soil")]
ibutton_info2$serial <-recode(ibutton_info2$serial, '6.60E+130' = '6600000055E00121')#had to hard code because kept coming in as scietific notation
ibutton_info2$serial <-recode(ibutton_info2$serial, '7.90E+15' = '7900000078656141')#typo (extra space) in info csv
dat.out2 = right_join(ibutton_info2, dat.out, by = "serial")
table(complete.cases(dat.out2))
unique(substrRight(as.character(dat.out2$Date.Time), 8)) #substring to view minutes and seconds for measurements
## all times are either on the hour or one minute past
#so we can just delete the minutes and seconds and that would round to the nearest even hour
#testing converting this
as.POSIXct(strptime("7/29/20 2:00:01 PM", "%m/%d/%y %I:%M:%S %p"))
dat.out2$Date.Time2 = as.POSIXct(strptime(as.character(dat.out2$Date.Time), "%m/%d/%y %I:%M:%S %p"))
str(dat.out2)

dat.out2$Date.TimeR = as.POSIXct(round(dat.out2$Date.Time2, units = "hours"))
#okay, should be able to use this to merge things together
temp = subset(dat.out2, site_code == "BEMO")
temp2 = subset(temp, NSEW == "North")
temp3 = subset(temp2, air_or_soil == "soil")

test <- unique(dat.out2$site_code)

dat.combined = data.frame()
for(i in unique(dat.out2$site_code)){
  #i = "BEMO"
  temp = subset(dat.out2, site_code == i)
  
  for(x in unique(levels(as.factor(as.character(temp$NSEW))))){
    #x = "North"
    temp2 = subset(temp, NSEW == x)
    
    temp.air = subset(temp2, air_or_soil == "air")
    temp.air.t = subset(temp.air, Unit == "C") #that is our base data
    #trim out some garbage columns
    temp.air.t$airT = temp.air.t$Value
    temp.air.t = temp.air.t[c("site", "site_code", "NSEW", "Date.TimeR", "airT")]
    
    #isolate humidity
    temp.air.h = subset(temp.air, Unit == "%RH")[c("Date.TimeR", "Value")]
    temp.air.h$airH = temp.air.h$Value
    temp.air.h = temp.air.h[c("Date.TimeR", "airH")]
    temp.air.out = merge(temp.air.t, temp.air.h, by = "Date.TimeR")
    
    #isolate soil temp
    temp.soil.t = subset(temp2, air_or_soil == "soil")
    temp.soil.t$soilT = temp.soil.t$Value
    temp.soil.t = temp.soil.t[c("Date.TimeR", "soilT")]
    temp.out = merge(temp.air.out, temp.soil.t, by = "Date.TimeR", all = T)
    
    dat.combined = rbind(dat.combined, temp.out)
  }#NSEW loop
}#site loop end

str(dat.combined)
unique(dat.combined$site_code)

#BAHA and WEMO were deployed June 21, 2021. The rest on 6/23/21?
#lets look at a few plots to confirm this.
temp  = subset(dat.combined, site_code == "BLWO")
temp2 = subset(temp, NSEW == "North")
temp3 = subset(temp2, Date.TimeR <= "2021-08-03 16:00:00")
plot(airT ~ Date.TimeR, data = temp3)
abline(v = as.POSIXct("2021-06-23 22:00:00"))
#let's trim everything to  so 6/24/21 is the first day
dat.combined2 = subset(dat.combined, Date.TimeR >= "2021-06-24 00:00:00")

#now lets do some quality control
hist(dat.combined2$airT)
hist(dat.combined2$airH)
#will need to decide how to handle values above 100
hist(dat.combined2$soilT)
#looks like some warm soilT days that should be checked too.

#lets do this by site
dat.combined2$site_NSEW = as.factor(paste(dat.combined2$site_code, dat.combined2$NSEW, sep = "_"))
str(dat.combined2)


table(dat.combined2$site_code, dat.combined2$NSEW)
plot(dat.combined2$soilT ~ dat.combined2$site_NSEW)




#################################################################
#lets try to correct for saturation drift.

dat.out = data.frame()

for(i in unique(levels(dat.combined2$site_NSEW))){
  #https://datasheets.maximintegrated.com/en/ds/DS1923.pdf
  #i = "BECL_North"
  temp = subset(dat.combined2, site_NSEW == i)
  temp$ov70 = ifelse(temp$airH >= 70, 2, 0) #hours over 70 since last reading
  temp$ov70na = ifelse(temp$airH >= 70, 2, NA) #mark out the times not over 70 with NA
  #calculate the cummulative hours of 70 for each event.
  temp$ov70cum = with(temp, ave(ov70,cumsum(is.na(ov70na)),FUN=cumsum))
  
  #calculate the partial correction factor
  temp$hcor.par = (0.0156 * temp$airH * 2.54^(-0.3502 * temp$ov70cum)) / (1 + (temp$airT - 25)/100)
  #make it 0 if not during a period over 70
  temp$hcor.par2 = ifelse(temp$ov70cum > 0, temp$hcor.par, 0)
  
  #now sum those up for each event.
  temp$hcor.cum = with(temp, ave(hcor.par2,cumsum(is.na(ov70na)),FUN=cumsum))
  temp$airHc = temp$airH - temp$hcor.cum
  
  #still a few over 100 so then move all those to just 100%.
  temp$airHc[temp$airHc >= 100] <- 100 
  hist(temp$airHc)
  
  #remove old columns
  temp2 = subset (temp, select = -c(ov70, ov70na, ov70cum, hcor.par, hcor.par2, hcor.cum))
  dat.out = rbind(dat.out, temp2)   
  
}
########################################################

dat.out$soilT = as.numeric(dat.out$soilT)

#now we need to check data for each site

for(i in unique(levels(dat.out$site_NSEW))){
  # i = "BEMO_North"
  png(filename = paste("2021_2022_iButtons", paste(i,".png", sep = ""), sep = "_"), width = 9, height = 6, units = "in", res = 300)
  temp = subset(dat.out, site_NSEW == i)
  
  par(oma=c(3,3,3,1), mar=c(2,3,1,0.5), mgp = c(3, 1, 0)  )
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow  = T))
  layout.show(n=6)
  #air temp examination
  hist(temp$airT, main = "")
  mtext(side = 1, line = 2, "Air Temp (C)")
  
  plot(airT ~ Date.TimeR, data = temp, type = "l")
  abline(h = 0)
  mtext(side = 2, line = 2, "Air Temp (C)")
  
  
  #humidity examination
  hist(temp$airHc, main = "")
  mtext(side = 1, line = 2, "Relative Humidity (%)")
  
  plot(airHc ~ Date.TimeR, data = temp, type = "l")
  mtext(side = 2, line = 2, "Relative Humidity (%)")
  
  #soil temp
  if(length(na.omit(temp$soilT)) >= 2){
    hist(temp$soilT, main = "")
    mtext(side = 1, line = 2.5, "Soil Temp (C)")
    
    temp.soil = temp[!is.na(temp$soilT), ]
    plot(soilT ~ Date.TimeR, data = temp.soil, type = "l")
    abline(h = 0)
    mtext(side = 2, line = 2, "Soil Temp (C)")
    mtext(side = 1, line = 2.5, "Date")
    
  }
  mtext(side = 3, line = 1, i, outer = T)
  mtext(side = 2, line = 1, "Frequency", outer = T)
  
  
  dev.off()
}

###############################
#now calculate and output the daily and monthly datasets
dat.daily = data.frame()
dat.monthly = data.frame()

for(i in unique(levels(dat.out$site_NSEW))){
  # i = "BEMO_North"
  temp = subset(dat.out, site_NSEW == i)
  temp$day = as.POSIXct(round(as.POSIXct(temp$Date.TimeR), "day"))
  
  #air temperature
  #min
  dailymin = tapply(temp$airT, temp$day, min) #minimum temperature each day
  dailyrecord = data.frame(day=names(dailymin), airTmin=dailymin)
  #max#
  dailymax = tapply(temp$airT, temp$day, max)
  dailyrecord$airTmax = dailymax
  #average#
  dailyave = tapply(temp$airT, temp$day, mean)
  dailyrecord$airTave = dailyave
  
  #air humidity
  #min
  dailymin = tapply(temp$airHc, temp$day, min) #minimum humidity each day
  dailyrecord$airHmin = dailymin
  #max#
  dailymax = tapply(temp$airHc, temp$day, max)
  dailyrecord$airHmax = dailymax
  #average#
  dailyave = tapply(temp$airHc, temp$day, mean)
  dailyrecord$airHave = dailyave
  
  #VPD calculations
  temp$satvap = exp((16.78*temp$airT - 116.9) / (temp$airT + 237.3))
  temp$actvap = temp$satvap * (temp$airHc/100)
  temp$vpd = temp$satvap - temp$actvap#units are kilopascals
  temp$vpd = ifelse(temp$airT > 0, temp$vpd, NA)
  
  #min
  dailymin = tapply(temp$vpd, temp$day, min, na.rm = T) #minimum humidity each day
  dailyrecord$airVPDmin = dailymin    
  dailyrecord$airVPDmin = ifelse(is.infinite(dailyrecord$airVPDmin), NA, dailyrecord$airVPDmin)
  
  #max#
  dailymax = tapply(temp$vpd, temp$day, max, na.rm =T)
  dailyrecord$airVPDmax = dailymax
  dailyrecord$airVPDmax = ifelse(is.infinite(dailyrecord$airVPDmax), NA, dailyrecord$airVPDmax)
  
  #average#
  dailyave = tapply(temp$vpd, temp$day, mean, na.rm = T)
  dailyrecord$airVPDave = dailyave
  dailyrecord$airVPDave = ifelse(is.infinite(dailyrecord$airVPDave), NA, dailyrecord$airVPDave)
  
  #soil temperature
  #min
  dailymin = tapply(temp$soilT, temp$day, min, na.rm = T) #minimum temperature each day
  dailyrecord$soilTmin = dailymin
  #max#
  dailymax = tapply(temp$soilT, temp$day, max, na.rm = T)
  dailyrecord$soilTmax = dailymax
  #average#
  dailyave = tapply(temp$soilT, temp$day, mean,na.rm = T)
  dailyrecord$soilTave = dailyave
  
  dailyrecord$site_NSEW = i
  
  ##############
  #now for monthly calculations
  dailyrecord$month = strftime(dailyrecord$day, format = "%Y-%m")
  monthlyairTave = tapply(dailyrecord$airTave, dailyrecord$month, mean)
  monthlyrecord = data.frame(month=names(monthlyairTave), airTave = monthlyairTave)
  monthlyrecord$airTmin = tapply(dailyrecord$airTmin, dailyrecord$month, mean)
  monthlyrecord$airTmax = tapply(dailyrecord$airTmax, dailyrecord$month, mean)
  
  monthlyrecord$airHave = tapply(dailyrecord$airHave, dailyrecord$month, mean)
  monthlyrecord$airHmin = tapply(dailyrecord$airHmin, dailyrecord$month, mean)
  monthlyrecord$airHmax = tapply(dailyrecord$airHmax, dailyrecord$month, mean)
  
  monthlyrecord$airVPDave = tapply(dailyrecord$airVPDave, dailyrecord$month, mean)
  monthlyrecord$airVPDmin = tapply(dailyrecord$airVPDmin, dailyrecord$month, mean)
  monthlyrecord$airVPDmax = tapply(dailyrecord$airVPDmax, dailyrecord$month, mean)
  
  monthlyrecord$soilTave = tapply(dailyrecord$soilTave, dailyrecord$month, mean)
  monthlyrecord$soilTmin = tapply(dailyrecord$soilTmin, dailyrecord$month, mean)
  monthlyrecord$soilTmax = tapply(dailyrecord$soilTmax, dailyrecord$month, mean)
  
  monthlyrecord$length = tapply(dailyrecord$day, dailyrecord$month, function(x) length(unique(na.omit(x))))
  
  monthlyrecord$site_NSEW = i
  
  
  #output
  dat.daily = rbind(dat.daily, dailyrecord)
  dat.monthly = rbind(dat.monthly, monthlyrecord)
  
  
}

dat.daily21 <- dat.daily
dat.monthly21 <-dat.monthly

today = Sys.Date()
write.csv(dat.daily, file = paste(paste("2021_2022_dat_daily", today, sep = "_"), ".csv", sep = ""),
          row.names = F)
write.csv(dat.monthly, file = paste(paste("2021_2022_dat_monthly", today, sep = "_"), ".csv", sep = ""),
          row.names = F)


# Combining daily and monthly datasets -----------------------------------
dat.dailyFull <- rbind(dat.daily20, dat.daily21)
dat.monthlyFull <- rbind(dat.monthly20, dat.monthly21)

write.csv(dat.dailyFull, file = paste(paste("2020-2022_dat_daily", today, sep = "_"), ".csv", sep = ""),
          row.names = F)
write.csv(dat.monthlyFull, file = paste(paste("2020_2022_dat_monthly", today, sep = "_"), ".csv", sep = ""),
          row.names = F)

# Data Exploration --------------------------------------------------------
#explore some plots of the summary data
dat.monthly$month.p = as.Date(paste(as.character(dat.monthly$month), "01", sep = "-"),
                              format='%Y-%m-%d')
dat.daily$day.p = as.Date(as.character(dat.daily$day),
                              format='%Y-%m-%d')

hist(dat.monthly$airTave)
hist(dat.monthly$soilTave)
hist(dat.monthly$airVPDmax)
hist(dat.daily$airVPDmax)

str(dat.monthly)
plot(airTmin ~ month.p, data = dat.monthly)
plot(airTmax ~ month.p, data = dat.monthly)
plot(soilTmin ~ month.p, data = dat.monthly)
plot(airVPDmin ~ month.p, data = dat.monthly)


temp = subset(dat.daily, site_NSEW == "BEMO_North")

plot(airTmin ~ day.p, data = temp, type = "l", col = "blue")
points(airTmax ~ day.p, data = temp, type = "l", col = "red")
points(airTave ~ day.p, data = temp, type = "l", col = "green")
points(soilTave ~ day.p, data = temp, type = "l", col = "brown", lwd = 3)

temp2 = subset(temp, day.p < "2020-09-01")
plot(airTmin ~ day.p, data = temp2, type = "l", col = "blue", ylim = c(8, 32))
points(airTmax ~ day.p, data = temp2, type = "l", col = "red")
points(airTave ~ day.p, data = temp2, type = "l", col = "green")
points(soilTave ~ day.p, data = temp2, type = "l", col = "brown", lwd = 3)


temp = subset(dat.monthly, site_NSEW == "BEMO_North")

plot(airTmin ~ month.p, data = temp, type = "l", col = "blue")
points(airTmax ~ month.p, data = temp, type = "l", col = "red")
points(soilTmax ~ month.p, data = temp, type = "l", col = "black", lwd = 3)
points(soilTave ~ month.p, data = temp, type = "l", col = "brown", lwd = 3)
