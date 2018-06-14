library(RODBC)
library(dplyr)
library(tidyr)

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
tblArtifactLog <- sqlFetch(als, "tblArtifactLog", stringsAsFactors=FALSE)
tblVS_Analysis <- sqlFetch(als, "tblVS_Analysis", stringsAsFactors=FALSE)

#close connections to the odbc data source
odbcClose(als)

diagnostics <- tblVS_Analysis %>%
  filter(type < 999991) %>%
  filter(rim > 0 | base > 0)

diagnostics <- inner_join(diagnostics, tblVS_TypeDef)

diagnostics <- select(diagnostics, contextID, bagNo, CRNo, muestra, wareName, groupName, typeName, freq, rim, base)

write.csv(diagnostics, file = "diagnostic_CRNo.csv")
