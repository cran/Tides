###Typical time series
require(Tides)
pdf(file="./projecten/artikel-GGG2/Rpackage/doc/TideFig1.pdf",width=12,height=6)
plot(waterlevels$time,waterlevels$h,ylab="water level [mTAW]",type="l")
lines(waterlevels$time,waterlevels$h0,col="grey")
dev.off()

### Detail of typical time series and visualistion of algorithm
pdf(file="./projecten/artikel-GGG2/Rpackage/doc/TideFig2.pdf",width=12,height=6)
TCwl <- TidalCharacteristics(waterlevels)
N1 <- 5
N2 <- 12
date1 <- TCwl$h$time[match(N1,TCwl$h$N)]
date2 <- TCwl$h$time[match(N2+1,TCwl$h$N)]
wssub <- subset(waterlevels,time>date1&time<date2)
HLsub <- subset(TCwl$HL,time>date1&time<date2)
Hsub <- subset(HLsub,HL=="H")
as.POSIXct(diff(Hsub$time))
plot(TCwl$h$time,TCwl$h$h,ylab="water level [mTAW]",type="l",xlim=c(date1,date2),ylim=c(280,360),axes=F)
axis(2)
days <- as.POSIXct(strptime(paste("2007-03-",1:31),format="%F"))
axis.POSIXct(1,waterlevels$time,format="%a",at=days)
box()

hH <- subset(TCwl$h,HL=="H")
points(hH$time,hH$h,col="blue",pch=20)

T2 <- 5*60*60
pt <- wssub[95,]
pt1 <- data.frame(time=0,h=0)
pt1$time <- pt$time-T2
pt1$h <- approx(wssub$time,wssub$h,xout=pt1$time)$y
pt2 <- data.frame(time=0,h=0)
pt2$time <- pt$time+T2
pt2$h <- approx(wssub$time,wssub$h,xout=pt2$time)$y
points(c(pt1$time,pt$time,pt2$time),c(pt1$h,pt$h,pt2$h),pch=20,type="b")
text(c(pt1$time,pt$time-60*60,pt2$time+2*60*60),c(pt1$h-5,pt$h,pt2$h),c("h(t-T/2)","h(t)","h(t+T/2)"))




#lines(waterlevels$time,waterlevels$h0,col="grey")
#text()
#abline(v=subset(HL,HL=="H")$time,col="red")
dev.off()

pdf(file="./projecten/artikel-GGG2/Rpackage/doc/TideFig3.pdf",width=12,height=6)
plot(TidalCharacteristics(waterlevels))
dev.off()

