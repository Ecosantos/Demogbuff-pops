rm(list=ls())
install.packages("tidyverse")
install.packages("terra")
install.packages("openxlsx")

library(tidyverse)
library(terra)
library(openxlsx)

dir.create(file.path("Data/ChelsacrutsData/Tmax"), showWarnings = FALSE)
dir.create(file.path("Data/ChelsacrutsData/Tmin"), showWarnings = FALSE)
dir.create(file.path("Data/ChelsacrutsData/Prec"), showWarnings = FALSE)


# LOAD POPULATION TO HAVE DATA EXTRACTED
linkMetadata<- "https://github.com/Ecosantos/Demogbuff-pops/raw/refs/heads/main/Data/CleanData.RDS"
CleanData<- readRDS(url(linkMetadata, method="libcurl"))


Metadata<-CleanData$Metadata
MetadataClean<-CleanData$MetadataClean

#Reduce (more) data to improve redability
MetadataFinal<-MetadataClean%>%select(-c(lambda,Ecoregion,Binomial))

MetadataFinal<-MetadataFinal%>%
  left_join(.,
            Metadata%>%select(ID,StudyStart, StudyDuration, StudyEnd)%>%distinct(),
            by="ID")


ClimaticData<-Metadata%>%
  dplyr::select(ID,Lat,Lon,SpeciesAccepted,StudyStart,StudyEnd)%>%
  distinct(ID,.keep_all=TRUE)%>%
  filter(!is.na(Lat)&  between(Lat, -90, 90))%>%
  filter(!is.na(Lon) &  between(Lon, -180, 180))

ClimaticData%>%
  arrange(Lon)

#Hlcn.72_716  - Exemplo do brasil +-20ºC

ClimaticData<-ClimaticData%>%
  #filter(SpeciesAccepted=="Clinocottus analis")%>%
  mutate(Lat=case_when(
    SpeciesAccepted=="Clinocottus analis" ~ 	32.82162,
    SpeciesAccepted=="Ursus maritimus" ~    	69.35431,
    SpeciesAccepted=="Nephtys incisa" ~     	41.26105,
    SpeciesAccepted=="Amphimedon compressa" ~ 21.09022,
    TRUE~Lat),
    Lon=case_when(
      SpeciesAccepted=="Clinocottus analis" ~  -117.27783,
      SpeciesAccepted=="Ursus maritimus" ~     -135.44041,
      SpeciesAccepted=="Nephtys incisa" ~      -72.879102,
      SpeciesAccepted=="Amphimedon compressa" ~ 70.792523,
      TRUE~Lon))



#=====================================================================================================
# LOAD A SAMPLING
#=====================================================================================================
sample_link<-"https://os.zhdk.cloud.switch.ch/chelsav1/chelsa_cruts/tmin/CHELSAcruts_tmin_9_1990_V.1.0.tif"
download.file(sample_link,destfile = "Data/ChelsacrutsData/Chelsa sample.tif",mode="wb")

sampleraster<-terra::rast("Data/ChelsacrutsData/Chelsa sample.tif")

Data_point<-terra::vect(data.frame(ClimaticData[,c(1,3,2)]), geom=c("Lon", "Lat"), crs=crs(sampleraster), keepgeom=TRUE)
Data_point<-terra::project(Data_point, crs(sampleraster))



#=====================================================================================================
# Check extraction
#=====================================================================================================

Data_point_example<-terra::extract(sampleraster, Data_point, na.rm=TRUE,bind=T)


#Fix population ID
Climate_dataID<-data.frame(
  ID=Data_point_example$ID,
  IDindex=terra::extract(sampleraster, Data_point, na.rm=TRUE,bind=F)[,1])
  

# Once checked, remove file to avoid problems with GITHUB
# Github doesn't work well with files close to 100mb in a standard account
file.remove("Data/ChelsacrutsData/Chelsa sample.tif")

#=====================================================================================================
# CLIMATIC DATA
#=====================================================================================================
LinksChelsa<-"https://github.com/Ecosantos/Demogbuff-pops/raw/refs/heads/IntegratingGoogleCollab/Data/ChelsacrutsData/ChelsacrutsDataV2.xlsx"


# LOAD LIST OF CLIMATIC DATA TO BE EXTRACTED
Tmin<- openxlsx::read.xlsx(xlsxFile = LinksChelsa, sheet = "Tmin", colNames = TRUE,skipEmptyRows = FALSE)
Tmax<- openxlsx::read.xlsx(xlsxFile = LinksChelsa, sheet = "Tmax", colNames = TRUE,skipEmptyRows = FALSE)
Prec<- openxlsx::read.xlsx(xlsxFile = LinksChelsa, sheet = "Preciptation", colNames = TRUE,skipEmptyRows = FALSE)



Tmin_df<-Tmin%>%
  mutate(param=gsub(Link, pattern=".*/CHELSAcruts_tmin_", replace="\\"))%>%
  mutate(Month=str_split(param, "_", n = 3, simplify = TRUE)[,1])%>%
  mutate(Year=str_split(param, "_", n = 3, simplify = TRUE)[,2])

Tmax_df<-Tmax%>%
  mutate(param=gsub(Link, pattern=".*/CHELSAcruts_tmax_", replace="\\"))%>%
  mutate(Month=str_split(param, "_", n = 3, simplify = TRUE)[,1])%>%
  mutate(Year=str_split(param, "_", n = 3, simplify = TRUE)[,2])

Prec_df<-Prec%>%
  mutate(param=gsub(Link, pattern=".*/CHELSAcruts_prec_", replace="\\"))%>%
  mutate(Month=str_split(param, "_", n = 3, simplify = TRUE)[,1])%>%
  mutate(Year=str_split(param, "_", n = 3, simplify = TRUE)[,2])

#---------------------------------------------------------------------------
# Create data frames to receive extracted information
#---------------------------------------------------------------------------

Prec_df_out<-Tmax_df_out<-Tmin_df_out<-data.frame(matrix(NA, ncol=length(Data_point), nrow=dim(Tmin_df)[1]))
chelsarast<-NULL
destfile<-NULL

#---------------------------------------------------------------------------
# MINIMUM TEMPERATURE
#---------------------------------------------------------------------------

for( i in 1:dim(Tmin_df)[1]){
  #Determine a name for raster file
  destfile = paste0("Data/ChelsacrutsData/Tmin/CHELSA_tmin_m",
                    Tmin_df$Month[i],"y",Tmin_df$Year[i],".tif")
  #Download raster
  download.file(Tmin_df$Link[i],destfile,mode="wb")
  #Routine
  chelsarast<-terra::rast(destfile)# Open raster
  Tmin_df_out[i,]<-terra::extract(chelsarast, Data_point, na.rm=TRUE)[,2] #Extract values
  #Delete raster  
  file.remove(destfile)
}

colnames(Tmin_df_out)<-Climate_dataID$ID

#saveRDS(Tmin_df_out,file="Data/ChelsacrutsData/MinTemperatureChelsa.rds")


#----------------------------------------------------------------------------------------
# Tmax
#----------------------------------------------------------------------------------------

chelsarast<-destfile<-NULL

for( i in 1:dim(Tmax_df)[1]){
  #Determine a name for raster file
  destfile = paste0("Data/ChelsacrutsData/CHELSA_tmax_m",
                    Tmax_df$Month[i],"y",Tmax_df$Year[i],".tif")
  #Download raster
  download.file(    Tmax_df$Link[i],destfile,mode="wb")
  #Routine
  #  chelsarast<-terra::rast(destfile)# Open raster
  #  Tmin_df_out[i,]<-terra::extract(chelsarast, Data_point, na.rm=TRUE)[,2] #Extract values
  print(i)
  #Delete raster  
  file.remove(destfile)
}

colnames(Tmax_df_out)<-Climate_dataID$ID
saveRDS(Tmax_df_out,file="Data/ChelsacrutsData/MaxTemperatureChelsa.rds")


#----------------------------------------------------------------------------------------
# Preciptation
#----------------------------------------------------------------------------------------
chelsarast<-destfile<-NULL

for( i in 1:dim(Prec_df)[1]){
  #Determine a name for raster file
  destfile = paste0("Data/ChelsacrutsData/CHELSA_prec_m",
                    Prec_df$Month[i],"y",Prec_df$Year[i],".tif")
  #Download raster
  download.file(    Prec_df$Link[i],destfile,mode="wb")
  #Routine
  chelsarast<-terra::rast(destfile)# Open raster
  Prec_df_out[i,]<-terra::extract(chelsarast, Data_point, na.rm=TRUE)[,2] #Extract values
  #Delete raster  
  file.remove(destfile)
}

colnames(Prec_df_out)<-Climate_dataID$ID
saveRDS(Prec_df_out,file="Data/ChelsacrutsData/PrecChelsa.rds")


TESTE<-readRDS(file.choose())

colnames(TESTE)<-Data_point$ID

TESTE[1:5,1:10]

cbind(Prec_df,TESTE)%>%
  pivot_longer(!c(Link:Year),names_to="ID",values_to="Prec")



hist(filter(
  cbind(Prec_df,TESTE)%>%
  pivot_longer(!c(Link:Year),names_to="ID",values_to="Prec"),
  ID=="Hlcn.72_716")$Prec)
