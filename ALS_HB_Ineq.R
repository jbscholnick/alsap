library(dplyr)
library(vegan)
library(vioplot)
library(RODBC)
library(tidyr)
library(ineq)
library(ggplot2)

#Variables for Altar Project burial and cache analysis
#burial goods (richness, diversity measure, and perhaps an index)
#cache goods (richness, diversity measure, and perhaps an index)
#biophysical markers (cranial modification, dental decoration)
#aggregate by time period (phase or period)


#open connection to the ALSAP database
als <- odbcConnect(dsn="alsap_fmp_db", uid = "ADMIN", pwd = "")

#create dataframes from tables
tblFeatureArtifact <- sqlFetch(als, "tblFeatureArtifact", stringsAsFactors=FALSE)
tblFeature <- sqlQuery(als, "SELECT * FROM tblFeature WHERE featureType = 'burial' and featurePhase != 'Unassigned'", stringsAsFactors=FALSE)
tblHB <- sqlFetch(als, "tblHB_Analysis 2", stringsAsFactors=FALSE)
#tblHBiso <- sqlFetch(als, "tblHB_Isotopes", stringsAsFactors=FALSE)
tblRV <- sqlFetch(als, "tblVS_Reconstructible 2", stringsAsFactors=FALSE)
tblVS_TypeDef <- sqlFetch(als, "tblVS_TypeDef", stringsAsFactors=FALSE)
tblContext <- sqlFetch(als, "tblContext", stringsAsFactors=FALSE)
tblContextFeature <- sqlFetch(als, "tblContextFeature", stringsAsFactors=FALSE)

#close connections to the odbc data source
odbcClose(als)

#read the data files-----------------------------------------
phase <- read.csv(file='~/Documents/AltarProject/Data/ALS_Phase.csv', header=T, sep = ",", na.strings = "NA", stringsAsFactors = F)
period <- c("Middle Preclassic", "Late Preclassic", "Protoclassic", "Early Classic", "Late Classic", "Terminal Classic")

#aggregate tblFeatureArtifact for artifact totals and artifact classes
tblFeatureArtifact <- tblFeatureArtifact %>%
  mutate(prov = case_when(materialClass == "JD" ~ "nonlocal",
                          materialClass == "OB" ~ "nonlocal",
                          materialClass == "SH" ~ "nonlocal",
                          materialClass == "AS" ~ "nonlocal",
                          materialClass == "SS" ~ "nonlocal",
                          TRUE ~ "local"))

nonlocal <- tblFeatureArtifact %>%
  group_by(featureNumber, prov) %>%
  summarise(nonlocal = sum(frequency)) %>%
  filter(prov == "nonlocal")

artifacts <- tblFeatureArtifact %>%
  group_by(featureNumber) %>%
  summarize(artifact_rich = n_distinct(materialClass), artifact_ct = sum(frequency))

#Relative abundance index of non-local burial goods = sum of jade, obsidian, shell (spondylus and olivella only), 
#and stingray spine / numerator + ceramic vessels


#rename tblRV$BurialCacheID
tblRV$featureNumber <- tblRV$BurialCacheID

#aggregate RVs
rvs <- tblRV %>%
  group_by(featureNumber, Form) %>%
  summarize(RVct = n()) %>%
  summarize(vsForm_rich = n_distinct(Form), vs_ct = sum(RVct))

rvs$vsPresence <- ifelse(rvs$vs_ct > 0, c(1), c(0)) 

#drop unnecessary variables in tblFeature and tblContext
varsFea <- c("featureNumber","featurePhase")
feature <- tblFeature[varsFea]
varsCxt <- c("contextID","op")
context <- tblContext[varsCxt]

context$housemound <- ifelse(context$op =="D", c("house mound"), c("ceremonial core")) 
feature <- left_join(feature, tblContextFeature, by = "featureNumber")
feature <- left_join(feature, context, by = "contextID")

#join feature, artifact, and RV tables
burialArtifacts <- left_join(feature, artifacts, by = "featureNumber")

burialArtifacts <- left_join(burialArtifacts, rvs, by = "featureNumber")

#remove phase = unassigned burials
HB <- full_join(burialArtifacts, tblHB[which(tblHB$Phase != "Unassigned"),], by = c("featureNumber" = "BurialID"))
HB <- inner_join(HB, phase, by = c("featurePhase"="Phase"))
HB$artifact_rich[is.na(HB$artifact_rich)] <- 0
HB$artifact_ct[is.na(HB$artifact_ct)] <- 0
HB$vsForm_rich[is.na(HB$vsForm_rich)] <- 0
HB$vs_ct[is.na(HB$vs_ct)] <- 0
HB$vsPresence[is.na(HB$vsPresence)] <- 0

#recode the cranial modification
HB$CranModScore <- NA
HB$CranModScore[HB$CranMod=="yes"] <- 100
HB$CranModScore[HB$CranMod=="Yes"] <- 100
HB$CranModScore[HB$CranMod=="no"] <- 0
HB$CranModScore[HB$CranMod=="No"] <- 0

HB$DentalDecPA <- NA
HB$DentalDecPA <- recode(as.factor(HB$DentalDec), "FI" = 1, "JI" = 1, "PI" = 1, "NN" = 0)
HB <- mutate(HB, DentalDecScore = ifelse(FinalAge != "child" & FinalAge != "infant",  100*(DentalDecQty/32), 0))

#recode the burial type
HB$GraveTypeScore <- recode(HB$GraveType, "Simple burial" = 1, "Urn burial" = 25, "Cyst burial" = 50, "Crypt burial" = 100)

#richness, sample size, cranial/dental decoration
HB <- mutate(HB, artVS_rich = artifact_rich + vsPresence, artVS_ct = artifact_ct + vs_ct, cranial_dental = CranModScore + DentalDecScore) 

#make contingency table for sex by period
periodHB_sex <- HB %>%
  group_by(Period, FinalSex) %>%
  summarise(ct = n()) %>%
  spread(FinalSex, ct)

#make contingency table for age by period
periodHB_age <- HB %>%
  group_by(Period, FinalAge) %>%
  summarise(ct = n()) %>%
  spread(FinalAge, ct)

#make contingency table for location by period
periodHB_location <- HB %>%
  group_by(Period, housemound) %>%
  summarise(ct = n()) %>%
  spread(housemound, ct)

period_HB_age_sex <- inner_join(periodHB_age, periodHB_sex, by = "Period")
period_HB_age_sex_loc <- inner_join(period_HB_age_sex, periodHB_location, by = "Period")
write.csv(period_HB_age_sex_loc, file="period_HB_age_sex_loc.csv")

#####Calculate Health Index (Steckel et al. 2002)
#Exclude infant, child and adolescent burials from cranial mod, dental dec, and health measures
#Use FinalAge
#Teeth - must have at least 8 teeth (or tooth loss greater than 24)

#Enamel Hypoplasia (0 - 4 recoded as 100, 75, 50, 25, 0)
HB$EnamelHypoplasiaS <- recode(HB$EnamelHypoplasia, '0'=100, '1'=75, '2'=50, '3'=25, '4'=0)
#Porotic Hypertosis/Anemia (0 - 4 recoded as 100, 75, 50, 25, 0)
HB$PoroticHyperOstosisS <- recode(HB$PoroticHyperOstosis, '0'=100, '1'=75, '2'=50, '3'=25, '4'=0)
#Dental Abscesses (recode values 0, 1, 2+ as 100, 50, 0)
HB$DentalAbsecessesS <- recode(HB$DentalAbsecesses, '0'=100, '1'=50, '2'=0, '3'=0)
#Osteitis (0 - 4 recoded as 100, 67, 33, 0, 0)
HB$OsteitisS <- recode(HB$Osteitis, '0'=100, '1'=75, '2'=50, '3'=25, '4'=0)
#Arthritis (0 - 1 recoded as 100, 0)
HB$ArthritisS <- recode(HB$Arthritis, '0'=100, '1'=0)
#Trauma/Healed Fractures  (0 - 1 recoded as 100, 0)
HB$HealedFracturesS <- recode(HB$HealedFractures, '0'=100, '1'=0)

#Dental health = 1 - (# caries + # premortem tooth loss)/(# teeth)
# number of teeth = 32
HB <- HB %>%
  mutate(DentCompleteness = 100*(1 - (Caries + PremortemToothLoss)/32))

HB <- mutate(HB, HealthIndex = 
               ifelse(is.na(EnamelHypoplasiaS), 0, EnamelHypoplasiaS) 
             + ifelse(is.na(PoroticHyperOstosisS), 0, PoroticHyperOstosisS)
             + ifelse(is.na(OsteitisS), 0, OsteitisS) 
             + ifelse(is.na(ArthritisS), 0, ArthritisS) 
             + ifelse(is.na(HealedFracturesS), 0, HealedFracturesS) 
             + ifelse(PremortemToothLoss < 24 | is.na(PremortemToothLoss), (ifelse(is.na(DentCompleteness), 0, .75*DentCompleteness) 
                                                + ifelse(is.na(DentalAbsecessesS), 0, .25*DentalAbsecessesS)), 0))

HB <- mutate(HB, PossHI = 
               ifelse(is.na(EnamelHypoplasiaS), 0, 100) 
             + ifelse(is.na(PoroticHyperOstosisS), 0, 100) 
             + ifelse(is.na(OsteitisS), 0, 100) 
             + ifelse(is.na(ArthritisS), 0, 100) 
             + ifelse(is.na(HealedFracturesS), 0, 100) 
             + ifelse(PremortemToothLoss < 24 | is.na(PremortemToothLoss), (ifelse(is.na(DentCompleteness), 0, 75) 
                                                                              + ifelse(is.na(DentalAbsecessesS), 0, 25)), 0))
HB <- mutate(HB, AdjHealthIndex = 
               ifelse(PossHI > 0, 100*HealthIndex/PossHI, NA)
             )

HB$Period=as.factor(HB$Period)
HB$Period=factor(HB$Period, levels=levels(HB$Period)[c(1,2,3,4)])

ineqArtifacts <- HB %>%
  #filter(FinalAge != c("child", "infant")) %>%
  group_by(Period) %>%
  summarize(N_burials = n(), gArtifactRich = ineq(artVS_rich, na.rm=TRUE), meanArtifactRich = mean(artVS_rich), 
            sdArtifactRich = sd(artVS_rich), medianArtifactRich = median(artVS_rich), maxArtifactRich = max(artVS_rich), 
            minArtifactRich = min(artVS_rich),gArtifactCount = ineq(artVS_ct, na.rm=TRUE), meanArtifactCount = mean(artVS_ct), 
            sdArtifactCount = sd(artVS_ct), medianArtifactCount = median(artVS_ct), maxArtifactCount = max(artVS_ct), 
            minArtifactCount = min(artVS_ct), gVSForm = ineq(vsForm_rich, na.rm=TRUE), meanVSFormRich = mean(vsForm_rich),
            sdVSFormRich = sd(vsForm_rich), medianVSFormRich = median(vsForm_rich), maxVSFormRich = max(vsForm_rich), 
            minVSFormRich = min(vsForm_rich), gVS_ct = ineq(vs_ct, na.rm=TRUE), meanVS_ct = mean(vs_ct), sdVS_ct = sd(vs_ct),
            medianVS_ct = median(vs_ct), maxVS_ct = max(vs_ct), minVS_ct = min(vs_ct), gCranialDental = ineq(cranial_dental, na.rm=TRUE), 
            meanCranialDental = mean(cranial_dental, na.rm=TRUE), sdCranialDental = sd(cranial_dental, na.rm=TRUE), 
            medianCranialDental = median(cranial_dental, na.rm=TRUE), minCranialDental = min(cranial_dental, na.rm=TRUE),
            maxCranialDental = max(cranial_dental, na.rm=TRUE), gGraveTypeScore = ineq(GraveTypeScore, na.rm=TRUE),
            meanGraveTypeScore = mean(GraveTypeScore, na.rm=TRUE), sdGraveTypeScore = sd(GraveTypeScore, na.rm=TRUE), 
            medianGraveTypeScore = median(GraveTypeScore, na.rm=TRUE), minGraveTypeScore = min(GraveTypeScore, na.rm=TRUE),
            maxGraveTypeScore = max(GraveTypeScore, na.rm=TRUE),
            gCranial = ineq(CranModScore, na.rm=TRUE), gDental = ineq(DentalDecScore, na.rm=TRUE))

ineqArtifactsHousemounds <- HB %>%
  group_by(housemound) %>%
  summarize(N_burials = n(), gArtifactRich = ineq(artVS_rich, na.rm=TRUE), meanArtifactRich = mean(artVS_rich), 
            sdArtifactRich = sd(artVS_rich), medianArtifactRich = median(artVS_rich), maxArtifactRich = max(artVS_rich), 
            minArtifactRich = min(artVS_rich),gArtifactCount = ineq(artVS_ct, na.rm=TRUE), meanArtifactCount = mean(artVS_ct), 
            sdArtifactCount = sd(artVS_ct), medianArtifactCount = median(artVS_ct), maxArtifactCount = max(artVS_ct), 
            minArtifactCount = min(artVS_ct), gVSForm = ineq(vsForm_rich, na.rm=TRUE), meanVSFormRich = mean(vsForm_rich),
            sdVSFormRich = sd(vsForm_rich), medianVSFormRich = median(vsForm_rich), maxVSFormRich = max(vsForm_rich), 
            minVSFormRich = min(vsForm_rich), gVS_ct = ineq(vs_ct, na.rm=TRUE), meanVS_ct = mean(vs_ct), sdVS_ct = sd(vs_ct),
            medianVS_ct = median(vs_ct), maxVS_ct = max(vs_ct), minVS_ct = min(vs_ct), gCranialDental = ineq(cranial_dental, na.rm=TRUE), 
            gGraveTypeScore = ineq(GraveTypeScore, na.rm=TRUE), gCranial = ineq(CranModScore, na.rm=TRUE), 
            gDental = ineq(DentalDecScore, na.rm=TRUE))

ineqHealth <- HB %>%
  filter(FinalAge != "child" & FinalAge != "infant" & AdjHealthIndex > 0) %>%
  group_by(Period) %>%
  summarize(N_Burials = n(), gAdjHealthIndex = ineq(AdjHealthIndex), mean(AdjHealthIndex), 
            sd(AdjHealthIndex), median(AdjHealthIndex), min(AdjHealthIndex), max(AdjHealthIndex))

ineqHealthHousemound <- HB %>%
  filter(FinalAge != "child" & FinalAge != "infant" & AdjHealthIndex > 0) %>%
  group_by(housemound) %>%
  summarize(N_Burials = n(), gAdjHealthIndex = ineq(AdjHealthIndex), mean(AdjHealthIndex), 
            sd(AdjHealthIndex), median(AdjHealthIndex), min(AdjHealthIndex), max(AdjHealthIndex))

#HBiso <- inner_join(HB, tblHBiso, by = c("featureNumber" = "BurialID"))
#ineqIso <- HBiso %>%
  #filter(FinalAge != c("child", "infant")) %>%
  #group_by(Period) %>%
  #summarize(N_Isotopes = n(), gSigmaC13 = ineq(abs(SigmaC13), na.rm=TRUE), meanSigmaC13 = mean(SigmaC13, na.rm=TRUE), sdSigmaC13 = sd(SigmaC13, na.rm=TRUE), medianSigmaC13 = median(SigmaC13, na.rm=TRUE), maxSigmaC13 = max(SigmaC13, na.rm=TRUE), minSigmaC13 = min(SigmaC13, na.rm=TRUE),
  #         gSigmaN15 = ineq(SigmaN15, na.rm=TRUE), meanSigmaN15 = mean(SigmaN15, na.rm=TRUE), sdSigmaN15 = sd(SigmaN15, na.rm=TRUE), medianSigmaN15 = median(SigmaN15, na.rm=TRUE), maxSigmaN15 = max(SigmaN15, na.rm=TRUE), minSigmaN15 = min(SigmaN15, na.rm=TRUE))

#HB_ineq <- inner_join(ineqArtifacts, ineqIso, by = "Period")


write.csv(ineqArtifacts, file="ineqArtifacts.csv")
write.csv(ineqHealth, file="ineqHealth.csv")
write.csv(ineqArtifactsHousemounds, file="ineqArtifactsHousemounds.csv")
write.csv(ineqHealthHousemound, file="ineqHealthHousemounds.csv")
write.csv(HB, file="HB.csv")

phases <- HB %>%
  group_by(Phase, featurePhase) %>%
  summarise(n = n())

#plot line graph of Material (artifact count), Relational (artifact class richness), and embodied weatlh (cranial/dental)
ineqArtifacts <- arrange(ineqArtifacts, desc(ineqArtifacts$Period))
pdf("inequality_by_period.pdf")
par(mar = c(8, 5, 3, 3), cex = 1.25, xpd = F, bg="white")
plot(ineqArtifacts$gArtifactCount, type="o", col="dark blue", pch = 1, axes=FALSE, ann=FALSE, ylim = c(0,1), lwd = 1.5, lty =1)
axis(1, at=c(1:4), lab=F)
text(x=1:4, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=c("Preclassic", "Early Classic", "Late Classic", "Terminal Classic"), srt=45, adj=1, xpd=TRUE)
axis(2, las=1, at=seq(from = 0, to = 1, by = 0.25))
box()
lines(ineqArtifacts$gArtifactRich, type="o", col = "dark green", pch =2, lwd = 1.5, lty =1)
lines(ineqArtifacts$gCranialDental, type="o", col = "black", pch =3, lwd = 1.5, lty =1)
lines(ineqHealth$gAdjHealthIndex, type = "o", col = "dark orange", pch= 4, lwd = 1.5, lty =1)
title(ylab = "Gini Coefficient")
legend("bottomright", c("Material wealth","Relational wealth", "Embodied wealth", "Health index"), cex=0.8, col=c("dark blue","dark green", "black", "dark orange"), pch=1:4, lwd=2, lty=1)
dev.off()

pdf("inequality_by_period_detail.pdf")
par(mar = c(8, 5, 3, 3), cex = 1.25, xpd = F, bg="white")
plot(ineqArtifacts$gArtifactCount, type="o", col="red3", pch = 15, axes=FALSE, ann=FALSE, ylim = c(0,1), lwd = 1.5, lty =1)
axis(1, at=c(1:4), lab=F)
text(x=1:4, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=c("Preclassic", "Early Classic", "Late Classic", "Terminal Classic"), srt=45, adj=1, xpd=TRUE)
axis(2, las=1, at=seq(from = 0, to = 1, by = 0.25))
box()
lines(ineqArtifacts$gVS_ct, type="o", col = "red3", pch =16, lwd = 1.5, lty =1)
lines(ineqArtifacts$gArtifactRich, type="o", col = "#6699CC", pch =17, lwd = 1.5, lty =1)
lines(ineqArtifacts$gVSForm, type="o", col = "#6699CC", pch =18, lwd = 1.5, lty =1)
lines(ineqArtifacts$gCranialDental, type="o", col = "#336633", pch =3, lwd = 1.5, lty =1)
lines(ineqHealth$gAdjHealthIndex, type="o", col = "#336633", pch =4, lwd = 1.5, lty =1)
title(ylab = "Gini Coefficient")
legend("topleft", c("Burial goods", "Ceramic vessels", "Material rich", "Vessel form rich", "Body mod. index", "Health index"), 
       col=c("red3", "red3", "#6699CC", "#6699CC", "#336633", "#336633"), pch=c(15:18,3:4), lwd=2, lty=1, cex = 0.7)
dev.off()

#join health index and artifact ineq measures
ineqHousemound <- inner_join(ineqArtifactsHousemounds, ineqHealthHousemound, by = "housemound")
ineqHousemound <- ineqHousemound %>%
  select(housemound, gArtifactCount, gVS_ct, gArtifactRich, gVSForm, gCranialDental, gAdjHealthIndex)

#transpose the df
ineqHousemound <- gather(ineqHousemound,  measure, gini, -housemound)

#readjust the order
#x$name <- factor(x$name, levels = x$name[order(x$val)])
ineqHousemound <- ineqHousemound[c(2,1,3)]
ineqHousemound$measure <- as.factor(ineqHousemound$measure)
ineqHousemound <- ineqHousemound %>%
  mutate(measure = case_when(measure == "gArtifactCount" ~ "Burial goods",
                   measure == "gVS_ct" ~ "Ceramic vessels",
                   measure == "gArtifactRich" ~ "Material rich",
                   measure == "gVSForm" ~ "Vessel form",
                   measure == "gCranialDental" ~ "Body mod. index",
                   measure == "gAdjHealthIndex" ~ "Health index",
                   TRUE                      ~  "other"
                   ))
ineqHousemound$measure=factor(ineqHousemound$measure, levels=c("Burial goods", "Ceramic vessels", "Material rich", "Vessel form", "Body mod. index", "Health index"))

#barplot
pdf("inequality_by_location.pdf")
ggplot(ineqHousemound, aes(x = measure, y= gini, fill = housemound)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  scale_fill_manual(values=c("#999999", "white"), labels=c("Ceremonial core", "House mound")) +
  xlab("") + ylab("Gini Coefficient") +
  theme_bw()+
  theme(axis.title.x = element_text(face="bold", colour="#990000"),
        axis.text.x  = element_text(angle= 90, vjust=0.5, size = 12),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        legend.title=element_blank())
dev.off()

#Make a boxplot of the health index
HB$Period=factor(HB$Period, levels=levels(HB$Period)[c(4,3,2,1)])
pdf("health_index_by_period.pdf")
par(mar = c(8, 5, 3, 3), cex = 1.25, xpd = F, bg="white")
boxplot(AdjHealthIndex~Period, data = HB, col="#6699CC", axes=FALSE, ann=FALSE, ylab="Health Index")
axis(1, at=c(1:4), lab=F)
box()
text(x=1:4, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=c("Preclassic", "Early Classic", "Late Classic", "Terminal Classic"), srt=45, adj=1, xpd=TRUE)
axis(2, las=1, at=seq(from = 0, to = 100, by = 20))
dev.off()

########___old_code___###############################################################################
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
