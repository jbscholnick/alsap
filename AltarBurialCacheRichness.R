library(dplyr)
library(vegan)
library(vioplot)
library(RODBC)

#Variables for Altar Project burial and cache analysis
#burial goods (richness, diversity measure, and perhaps an index)
#cache goods (richness, diversity measure, and perhaps an index)
#biophysical markers (cranial modification, dental decoration)
#aggregate by time period (phase or period)
rm(list=ls())
setwd('~/Documents/AltarProject/Analysis')


#open connection to the ALSAP database
als <- odbcConnect(dsn="alsap_fmp_db", uid = "ADMIN", pwd = "")

#create dataframes from tables
tblFeatureArtifact <- sqlFetch(als, "tblFeatureArtifact")
tblFeature <- sqlFetch(als, "tblFeature")
tblHB <- sqlFetch(als, "tblHB_Analysis 2")
tblHBiso <- sqlFetch(als, "tblHB_Isotopes")
tblRV <- sqlFetch(als, "tblVS_Reconstructible 2")

#close connections to the odbc data source
odbcClose()

#read the data files-----------------------------------------

ca <- read.csv(file='~/Documents/AltarProject/Data/CA_Analysis.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
hb <- read.csv(file='~/Documents/AltarProject/Data/HB_Analysis.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
cahb <- read.csv(file='~/Documents/AltarProject/Data/CAHB.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
phase <- read.csv(file='~/Documents/AltarProject/Data/ALS_Phase.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
period <- c("Middle Preclassic", "Late Preclassic", "Protoclassic", "Early Classic", "Late Classic", "Terminal Classic")

#remove phase = unassigned caches and burials
hb <- hb[which(hb$Phase != "Unassigned"),]
ca <- ca[which(ca$Phase != "Unassigned"),]
cahb <- cahb[which(cahb$Phase != "Unassigned"),]

ca <- merge(ca, phase, by = "Phase")
hb <- merge(hb, phase, by = "Phase")
cahb <- merge(cahb, phase, by = "Phase")
  
#make contingency tables by phase
ca_phase_mat <- xtabs(Frequency~Phase2+Material, data = ca)
ca_period_mat <- xtabs(Frequency~Period+Material, data = ca)
hb_phase_mat <- xtabs(Frequency~Phase2+Material, data = hb)
hb_period_mat <- xtabs(Frequency~Period+Material, data = hb)

#transform tables to binary
ca_phase_mat[ca_phase_mat > 0.0] <- 1
hb_phase_mat[hb_phase_mat > 0.0] <- 1
ca_period_mat[ca_period_mat > 0.0] <- 1
hb_period_mat[hb_period_mat > 0.0] <- 1

ca_phase_rich <- margin.table(ca_phase_mat[1:14,],1)
hb_phase_rich <- margin.table(hb_phase_mat[1:19,], 1)
ca_rich <- margin.table(ca_period_mat[1:5,],1)
hb_rich <- margin.table(hb_period_mat[1:5,], 1)

#ceramics per cache/burial
ca_vs <- ca %>%
  filter(Material == "VS", Period != "") %>%
  group_by(CacheID,Period) %>%
  summarise(vs_ct = sum(Frequency)) %>%
  arrange(Period) %>%
  group_by(Period) %>%
  summarise(ca_mean_vs = mean(vs_ct), ca_max_vs = max(vs_ct))
hb_vs <- hb %>%
  filter(Material == "VS", Period != "") %>%
  group_by(BurialID,Period) %>%
  summarise(vs_ct = sum(Frequency)) %>%
  arrange(Period) %>%
  group_by(Period) %>%
  summarise(hb_mean_vs = mean(vs_ct), hb_max_vs = max(vs_ct))

#cache versus burial
ca_ct <- cahb%>%
  filter(Type=="cache", Period != "")%>%
  group_by(Context, Period)%>%
  summarize(ct=n())%>%
  group_by(Period)%>%
  summarise(caches=n())
  
hb_ct <- cahb%>%
  filter(Type=="burial", Period != "")%>%
  group_by(Context, Period)%>%
  summarize(ct=n())%>%
  group_by(Period)%>%
  summarise(burials=n())

cahb_ct <- merge(ca_ct, hb_ct, by = "Period", all=T)
cahb_ct[is.na(cahb_ct)] <- 0
cahb_ct$cache_index <- cahb_ct$caches/(cahb_ct$caches+cahb_ct$burials)

#merge the cache and burial vessel count tables
cahb_vs <- merge(ca_vs, hb_vs, by = "Period", all=T)

#merge the cache/burial count and the vessel counts
cahb_analysis <- merge(cahb_ct, cahb_vs, by ="Period", all=T)

#append the richness by Period
cahb_analysis$ca_rich <- margin.table(ca_period_mat[1:5,],1)
cahb_analysis$hb_rich <- margin.table(hb_period_mat[1:5,],1)

#make boxplots of A) ceramic vessel count, B) richness by period and by deposit type
cahb_tab <- xtabs(~Context+Material, data = cahb)
cahb_tab <- cahb_tab[,1:17]

cahb_sample <- cahb %>%
  select(Context, Period, Type) %>%
  distinct(Context, Type, Period) %>%
  group_by(Period, Type) %>%
  summarise(count = n())

write.csv2(cahb_sample, file = "burial+cache_period.csv")

cahb_context <- cahb %>%
  select(Context, Type, Period, Material) %>%
  distinct(Context, Type, Period, Material) %>%
  group_by(Context, Type, Period) %>%
  summarise(rich = n())

cahb_context_vs <- cahb %>%
  select(Context, Period, Material) %>%
  filter(Material == "VS") %>%
  group_by(Context) %>%
  summarise(vs_ct = n())

cahb_context <- merge(cahb_context, cahb_context_vs, by = "Context")

cahb_context<- cahb_context[which(cahb_context$Period != ""),]

boxplot(cahb_context$rich~cahb_context$Period+cahb_context$Type, xlab = "richness")
boxplot(cahb_context$vs_ct~cahb_context$Period+cahb_context$Type, xlab = "ceramic counts")


pdf("richness.period.equal.pdf")
par(mar = c(4, 7, 3, 3), xpd = F, bg="white")
boxplot(cahb_context$rich~cahb_context$Period, ylim=c(0,10), horizontal = T, outline =T, varwidth=F, yaxt = "n")
axis(side = 2, at = 1:7, las = 2, labels = period, cex.axis = .8, srt=45)
mtext("richness", side = 1, line = 2, cex = 1)
dev.off()

x2 <- cahb_context$rich[cahb_context$Period == "2LatePreclassic"]
x3 <- cahb_context$rich[cahb_context$Period == "3Protoclassic"]
x4 <- cahb_context$rich[cahb_context$Period == "4EarlyClassic"]
x5 <- cahb_context$rich[cahb_context$Period == "5LateClassic"]
x6 <- cahb_context$rich[cahb_context$Period == "6TerminalClassic"]

pdf("richness.violin.pdf")
par(mar = c(5, 5, 4, 2), xpd = F, bg="white")
vioplot(x2, x3, x4, x5, x6, names = period, col = "grey", ylim=c(1,10))
mtext("richness", side = 2, line = 3, cex = 1)
dev.off()

pdf("ceramic.count.period.pdf")
par(mar = c(4, 7, 3, 3), xpd = F, bg="white")
boxplot(cahb_context$vs_ct~cahb_context$Period, xlab = "ceramic counts", horizontal = T, outline =T,  ylim=c(0,10), varwidth =T, yaxt="n")
axis(side = 2, at = 1:7, las = 2, labels = period, cex.axis = .8)
dev.off()

write.csv(cahb_analysis, file="cahb_analysis.csv")
