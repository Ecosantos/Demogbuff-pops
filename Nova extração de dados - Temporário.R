
rm(list=ls())

library(tidyverse)
library(terra)


# LOAD POPULATION TO HAVE DATA EXTRACTED
CleanData<-readRDS("Data/CleanData.RDS")

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

sampleraster<-terra::rast("Data/ChelsacrutsData/Chelsa sample.tif")
Data_point<-terra::vect(data.frame(ClimaticData[,c(1,3,2)]), geom=c("Lon", "Lat"), crs=crs(sampleraster), keepgeom=TRUE)
Data_point<-terra::project(Data_point, crs(sampleraster))


Data_point
plot(sampleraster)
plot(Data_point,add=T)
plot(Data_point)


#=====================================================================================================
# Check extraction
#=====================================================================================================

terra::extract(sampleraster, Data_point, na.rm=TRUE)[,2]


#=====================================================================================================
# CLIMATIC DATA
#=====================================================================================================

# LOAD LIST OF CLIMATIC DATA TO BE EXTRACTED
Tmin<- openxlsx::read.xlsx(xlsxFile = "Data/ChelsacrutsData/ChelsacrutsDataV2.xlsx", sheet = "Tmin", colNames = TRUE,skipEmptyRows = FALSE)
Tmax<- openxlsx::read.xlsx(xlsxFile = "Data/ChelsacrutsData/ChelsacrutsDataV2.xlsx", sheet = "Tmax", colNames = TRUE,skipEmptyRows = FALSE)
Prec<- openxlsx::read.xlsx(xlsxFile = "Data/ChelsacrutsData/ChelsacrutsDataV2.xlsx", sheet = "Preciptation", colNames = TRUE,skipEmptyRows = FALSE)



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



Tmin_df%>%head()
Tmax_df%>%head()
Prec_df%>%head()


#---------------------------------------------------------------------------
# Create data frames to receive extracted information
#---------------------------------------------------------------------------

Prec_df_out<-Tmax_df_out<-Tmin_df_out<-data.frame(matrix(NA, ncol=length(Data_point), nrow=dim(Tmin_df)[1]))
chelsarast<-NULL
destfile<-NULL

#---------------------------------------------------------------------------
# MINIMUM TEMPERATURE
#---------------------------------------------------------------------------

for( i in 1:3){
#Determine a name for raster file
    destfile = paste0("Data/ChelsacrutsData/tmin/CHELSA_tmin_m",
                    Tmin_df$Month[i],"y",Tmin_df$Year[i],".tif")
#Download raster
    download.file(    Tmin_df$Link[i],destfile,mode="wb")
#Routine
    chelsarast<-terra::rast(destfile)# Open raster
    Tmin_df_out[i,]<-terra::extract(chelsarast, Data_point, na.rm=TRUE)[,2] #Extract values
#Delete raster  
    file.remove(destfile)
    }

#saveRDS(Tmin_df_out,file="MinTemperatureChelsa_complete_04Jun23.rds")

dir("Data/ChelsacrutsData")


#----------------------------------------------------------------------------------------
# Tmax
#----------------------------------------------------------------------------------------

chelsarast<-destfile<-NULL

for( i in 1:3){
  #Determine a name for raster file
  destfile = paste0("Data/ChelsacrutsData/tmax/CHELSA_tmax_m",
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

#saveRDS(Tmin_df_out,file="MinTemperatureChelsa_complete_04Jun23.rds")



#----------------------------------------------------------------------------------------
# Preciptation
#----------------------------------------------------------------------------------------
chelsarast<-destfile<-NULL

for( i in 1:dim(Prec_df)[1]){
  #Determine a name for raster file
  destfile = paste0("Data/ChelsacrutsData/prec/CHELSA_prec_m",
                    Prec_df$Month[i],"y",Prec_df$Year[i],".tif")
  #Download raster
  download.file(    Tmin_df$Link[i],destfile,mode="wb")
  #Routine
  chelsarast<-terra::rast(destfile)# Open raster
  Tmin_df_out[i,]<-terra::extract(chelsarast, Data_point, na.rm=TRUE)[,2] #Extract values
  #Delete raster  
  file.remove(destfile)
}

#saveRDS(Tmin_df_out,file="MinTemperatureChelsa_complete_04Jun23.rds")


