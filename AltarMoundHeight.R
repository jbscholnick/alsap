library(ineq)

#read the data files-----------------------------------------
rm(list=ls())
setwd('~/Documents/AltarProject/Analysis')
elev <- read.csv2(file='~/Documents/AltarProject/Data/Elevs2.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
elev$Elev2 <- as.numeric(elev$Elev2)
base.area <- read.csv2(file='~/Documents/AltarProject/Data/MoundBaseArea.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
base.area$SHAPE_Area <- as.numeric(base.area$SHAPE_Area)


#plot histogram
pdf("histElevs.pdf")
hist(elev$Elev2, breaks = 10, col = "light grey", xlab="Mound Height (Meters)", main = "Distribution of Mound Heights", xlim=c(0,14))
dev.off()

pdf("histDensityElevs.pdf")
hist(elev$Elev2, breaks = 10, freq = F, xlab="Mound Height (Meters)", main = "Density of Mound Heights", xlim=c(0,14), ylim=c(0,0.3))
lines(density(elev$Elev2))
dev.off()

#density function plot
d <- density(elev$Elev2)
pdf("kernelDensityElevs.pdf")
plot(d, main = "Kernel Density of Mound Heights", col = "dark blue")
polygon(d, col = "dark blue")
dev.off()

#logged distribution
logd <- density(log(elev$Elev2))
pdf("logkernelDensityElevs.pdf")
plot(logd, main = "Log-tranformed Kernel Density of Mound Heights", col = "dark blue")
polygon(logd, col = "dark blue")
dev.off()

#ecdf
pdf("cdfElevs.pdf")
plot(ecdf(elev$Elev2), verticals=TRUE, pch=46)
dev.off()

#subset Harvard mounds
harvard_mounds <- elev[ which(elev$OBJECTID <= 81), ]
#subset house mounds
house_mounds <- elev[ which(elev$Location == "house_mound"), ]

#plot Lorenz curve 
pdf("lorenz_curve_mound_height.pdf")
par(mar = c(4.5, 4.5, 3, 3), xpd = F, bg="white")
plot(Lc(elev$Elev2), col = "red3", xlab="cumulative % of mounds", ylab="cumulative % of height", main = NULL, lwd=3, las=1, cex.axis = 1.25, cex.lab = 1.5)
dev.off()

#plot Lorenz curve 
pdf("lorenz_curve_mound_area.pdf")
par(mar = c(4, 4, 3, 3), xpd = F, bg="white")
plot(Lc(base.area$SHAPE_Area), xlab="cumulative % of mounds", ylab="cumulative % of area", main = NULL, lwd=2, las=1, cex = 1.5)
dev.off()

#harvard mounds only
plot(Lc(harvard_mounds$Elev2))
#house mounds only
plot(Lc(house_mounds$Elev2))

#gini index for mound height
gBaseArea <- ineq(base.area$SHAPE_Area, type = "Gini")
#gini index for mound height
gMound <- ineq(elev$Elev2)
#gini index for Harvard mound height
gHarvardMound <- ineq(harvard_mounds$Elev2)
#gini index for Harvard mound height
gHouseMound <- ineq(house_mounds$Elev2)
