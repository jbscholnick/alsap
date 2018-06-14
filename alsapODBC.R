library(RODBC)
odbcDataSources()

#open connection to the ALSAP database
als <- odbcConnect(dsn="alsap_fmp_db", uid = "ADMIN", pwd = "")

#create dataframes from tables
cxt <- sqlFetch(als, "tblContext")
baglist <- sqlFetch(als, "tblArtifactLog")

#close connections to the odbc data source
odbcClose()


