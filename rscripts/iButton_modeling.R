library(doBy)
sumfun <- function(x, ...){
  c(m=mean(x, ...), sd=sd(x, ...), l=length(x),
    sd.up = mean(x, ...) + sd(x, ...),
    sd.down = mean(x, ...) - sd(x, ...),
    se.up = mean(x, ...) + (sd(x, ...) / sqrt(length(x))),
    se.down = mean(x, ...) - (sd(x, ...) / sqrt(length(x))))
}

setwd("C:\\Users\\Jay Wason\\OneDrive - University of Maine System\\UMaine\\Research\\Grants\\Second Century 2019\\ibuttons\\")
dat.d = read.csv("dat_daily_2021-10-29.csv", header = T)

head(dat.d)
hist(dat.d$airVPDmax)
table(dat.d$site_NSEW)
temp = subset(dat.d, airVPDmax >= 2)
table(temp$site_NSEW)

#so, lets see if we can find anything unique about the climate of western mountain

plot(soilTmax ~ as.Date(day), data = dat.d, cex = 0.5)
points(soilTmax ~ as.Date(day), data = subset(dat.d, site_NSEW == "BEMO_South"),
       col = 'red', pch = 19, cex = 0.5, type = "l")

temp = subset(dat.d, month == "2020-08")
table(temp$site_NSEW)

plot(airVPDmax ~ as.Date(day), data = temp, cex = 0.5)
points(airVPDmax ~ as.Date(day), data = subset(temp, site_NSEW == "BEMO_South"),
       col = 'red', pch = 19, cex = 0.5, type = "l")
points(airVPDmax ~ as.Date(day), data = subset(temp, site_NSEW == "BEMO_North"),
       col = 'red', pch = 19, cex = 0.5, type = "l")


plot(airTmax   ~ as.Date(day), data = temp, cex = 0.5)
points(airTmax   ~ as.Date(day), data = subset(temp, site_NSEW == "BEMO_South"),
       col = 'red', pch = 19, cex = 0.5, type = "l")
points(airTmax   ~ as.Date(day), data = subset(temp, site_NSEW == "BEMO_North"),
       col = 'red', pch = 19, cex = 0.5, type = "l")

temp.sum = summaryBy(airTmax + airTave + airTmin + airVPDmax ~ day, data = temp, FUN = sumfun)



png(file = "temp_vpd.png", width = 7, height = 7, units = 'in', res = 300)
layout(matrix(c(1,2), 2, 1))
par(oma=c(3.5,3.5,0.5,0.5), mar=c(1,0.5,1,0.5), mgp = c(3, 0.6, 0)  )

plot(airTmax.m ~ as.Date(day), data = temp.sum, pch = NA,
     ylim = c(10,30), ylab = "Air temp (C)", xlab = "Day")
polygon(x = c(as.Date(temp.sum$day), rev(as.Date(temp.sum$day))), 
        y = c(temp.sum$airTmax.m, rev(temp.sum$airTmin.m)),
        col = "tomato3")
points(airTave.m ~ as.Date(day), data = temp.sum, type = "l")
mtext(side = 2, line = 2, "Air temp. (C)")

plot(airVPDmax.m ~ as.Date(day), data = temp.sum, pch = NA,
     ylim = c(0,3), ylab = "VPD (kPa)", xlab = "Day")
points(airVPDmax.m ~ as.Date(day), data = temp.sum, type = "l")
polygon(x = c(as.Date("2020-08-01"), as.Date(temp.sum$day), as.Date("2020-08-31")),
        y = c(0, temp.sum$airVPDmax.m, 0), col = "royalblue2")
mtext(side = 2, line = 2, "VPD (kPa)")
mtext(side = 1, line = 2, "Date")

dev.off()
