#============================================================================
# 		CLIMATIC DATA DOWNLOAD
# https://chelsa-climate.org/timeseries/
# Check also: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V1%2Fchelsa_cruts
#============================================================================
#----------------------------------------------------------------------------
#			TMIN
#----------------------------------------------------------------------------
Tmin<- openxlsx::read.xlsx(xlsxFile = "./ChelsacrutsData/ChelsacrutsDataV1.xlsx", sheet = "Tmin", colNames = TRUE,skipEmptyRows = FALSE)
Tmin[1,]

Tmin_df<-Tmin%>%
mutate(param=gsub(Link, pattern=".*/CHELSAcruts_tmin_", replace="\\"))%>%
mutate(Month=str_split(param, "_", n = 3, simplify = TRUE)[,1])%>%
mutate(Year=str_split(param, "_", n = 3, simplify = TRUE)[,2])

#for( i in 1:dim(Tmin_df)[1]){
#download.file(
#Tmax_df$Link[i],
#destfile = paste0("ChelsacrutsData/tmin/CHELSA_tmax_m",Tmin_df$Month[i],"y",Tmin_df$Year[i],".tif"),mode="wb")
#print(i)}

#----------------------------------------------------------------------------
#			TMAX
#----------------------------------------------------------------------------
Tmax<- openxlsx::read.xlsx(xlsxFile = "./ChelsacrutsData/ChelsacrutsDataV1.xlsx", sheet = "Tmax", colNames = TRUE,skipEmptyRows = FALSE)
Tmax[1,]

Tmax_df<-Tmax%>%
mutate(param=gsub(Link, pattern=".*/CHELSAcruts_tmax_", replace="\\"))%>%
mutate(Month=str_split(param, "_", n = 3, simplify = TRUE)[,1])%>%
mutate(Year=str_split(param, "_", n = 3, simplify = TRUE)[,2])


#for( i in 1:dim(Tmax_df)[1]){
#download.file(
#Tmax_df$Link[i],
#destfile = paste0("ChelsacrutsData/tmax/CHELSA_tmax_m",Tmax_df$Month[i],"y",Tmax_df$Year[i],".tif"),mode="wb")
#print(i)}

#----------------------------------------------------------------------------
#			PRECIPTATION
#----------------------------------------------------------------------------
Prec<- openxlsx::read.xlsx(xlsxFile = "./ChelsacrutsData/ChelsacrutsDataV1.xlsx", sheet = "Preciptation", colNames = TRUE,skipEmptyRows = FALSE)
Prec[1,]


Tmax_df<-Tmax%>%
mutate(param=gsub(Link, pattern=".*/CHELSAcruts_tmax_", replace="\\"))%>%
mutate(Month=str_split(param, "_", n = 3, simplify = TRUE)[,1])%>%
mutate(Year=str_split(param, "_", n = 3, simplify = TRUE)[,2])


#for( i in 1:dim(Tmax_df)[1]){
#download.file(
#Tmax_df$Link[i],
#destfile = paste0("ChelsacrutsData/tmax/CHELSA_tmax_m",Tmax_df$Month[i],"y",Tmax_df$Year[i],".tif"),mode="wb")
#print(i)}

#==============================================================================================
#			EXTRACT DATA AND STORE CLIMATIC DATA
#==============================================================================================

ClimaticData<-Metadata%>%
dplyr::select(ID,Lat,Lon,SpeciesAccepted,StudyStart,StudyEnd)%>%
distinct(ID,.keep_all=TRUE)

PrecChelsa<-TmaxChelsa<-TminChelsa<-NULL

climDir<-"C:/Artigos e resumos publicados submetidos ideias/3 - Em desenvolvimento/Demographic buffering continuum - Plants and animals/Data and script/ChelsacrutsData"

raster_data_tmax<-list.files(path=paste0(climDir,"/tmax"), pattern=".tif") 
raster_data_tmin<-list.files(path=paste0(climDir,"/tmin"), pattern=".tif") 
raster_data_prec<-list.files(path=paste0(climDir,"/prec"), pattern=".tif") 

#---------------------------------------------------------------------------------
#			AMMENDMENT TO REDUCE CONTACT WITH WATER
# Water superficience in CHELSA is -32768, tranform this to NA is very time consuming
# Alternatively, I moved the coordinates a bit to avoid touch the water surface. 
# It was necessary only for 4 animal as plants has enough data already
#---------------------------------------------------------------------------------
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


#---------------------------------------------------------------------------------
#			MAX TEMPERATURE
#---------------------------------------------------------------------------------
#TMAX
#TmaxChelsa
#for(i in 1:length(raster_data_tmax)){
#X<-raster::raster(paste0(paste0(climDir,"/tmax/"),raster_data_tmax[i]))
#TmaxChelsa[[i]]<-raster::extract(X,
#	ClimaticData[,c(3,2)],buffer=100, fun=mean)
#gc()
#print(i)}

Time<-str_replace(
		str_replace(raster_data_tmax,"CHELSA_tmax_",""),
					".tif","")

TmaxChelsa<-do.call(cbind,TmaxChelsa)%>%data.frame(.)
colnames(TmaxChelsa)<-Time

TmaxChelsa<-cbind(ClimaticData,TmaxChelsa)%>%
pivot_longer(!c(ID:StudyEnd),values_to="TMax",names_to="Time")

saveRDS(TmaxChelsa,file="MaxTemperatureChelsa_01Aug23.rds")


#TmaxChelsa%>%
#mutate(Year=str_split(Time,"y",simplify=T)[,2])%>%
#mutate(Month=str_replace(str_split(Time,"y",simplify=T)[,1],"m",""))

#---------------------------------------------------------------------------------
#			MIN TEMPERATURE
#---------------------------------------------------------------------------------
rm(Time)
Time<-str_replace(
		str_replace(raster_data_tmin,"CHELSA_tmin_",""),
					".tif","")


#TMIN
#for(i in 1:length(raster_data_tmin)){
#X<-raster::raster(paste0(paste0(climDir,"/tmin/"),raster_data_tmin[i]))
#TminChelsa[[i]]<-raster::extract(X,
#	ClimaticData[,c(3,2)],buffer=100, fun=mean)
#gc()
#print(i)}

TminChelsa<-do.call(cbind,TminChelsa)%>%data.frame(.)
colnames(TminChelsa)<-Time

TminChelsa<-cbind(ClimaticData,TminChelsa)%>%
pivot_longer(!c(ID:StudyEnd),values_to="TMin",names_to="Time")

saveRDS(TminChelsa,file="MinTemperatureChelsa_01Aug23.rds")


#---------------------------------------------------------------------------------
#			PRECIPTATION
#---------------------------------------------------------------------------------
rm(Time)
Time<-str_replace(
		str_replace(raster_data_prec,"CHELSA_prec_",""),
					".tif","")

#PRECIPTATION
#for(i in 1:length(raster_data_prec)){
#X<-raster::raster(paste0(paste0(climDir,"/prec/"),raster_data_prec[i]))
#PrecChelsa[[i]]<-raster::extract(X,
#	ClimaticData[,c(3,2)],buffer=100, fun=mean)
#gc()
#print(i)}



PrecChelsa<-do.call(cbind,PrecChelsa)%>%data.frame(.)
colnames(PrecChelsa)<-Time

PrecChelsa<-cbind(ClimaticData,PrecChelsa)%>%
pivot_longer(!c(ID:StudyEnd),values_to="Prec",names_to="Time")

saveRDS(PrecChelsa,file="PreciptationChelsa_01Aug23.rds")



#Check validity using Terra package
X1<-terra::rast(paste0(paste0(climDir,"/prec/"),raster_data_prec[1]))
Data_point<-terra::vect(data.frame(ClimaticData[,c(3,2)]), geom=c("Lon", "Lat"), crs=crs(Xtem), keepgeom=FALSE)
point_buf<-buffer(Data_point,15)

data.frame(ClimaticData,
MEAN=extract(X1, Data_point, buffer=1000,na.rm=TRUE,fun=mean,method="simple")[,2],
MAX=extract(X1, Data_point, buffer=1000,na.rm=TRUE,fun=max,method="simple")[,2],
MIN=extract(X1, Data_point, buffer=1000,na.rm=TRUE,fun=min,method="simple")[,2])

X1%>%plot()
Data_point%>%plot(.,add=T)
#plot(point_buf,add=T)
#plot(buffer(Data_point,5),add=T)
lines(buffer(Data_point,300000),lwd=.01,col="red")


#==============================================================================================
#			STORE DATA
#==============================================================================================
order_m<-c(10:12,1:9)

month_ord<-paste0("m",order_m)
year_ord<-paste0("y",1901:2016)

nyear<-year_ord%>%length()

yearMx<-data.frame(matrix(1,12,nyear))

dim(yearMx)[1]*dim(yearMx)[2]

colnames(yearMx)<-year_ord

climatic_ord_date<-data.frame(month_ord,yearMx)%>%
pivot_longer(!month_ord,names_to="year",values_to="values")%>%
mutate(values=paste0(month_ord,year))

climatic_ord_date

do.call(rbind,TmaxChelsa)%>%dim()
%>%
data.frame(.,Time=climatic_ord_date$values)%>%
column_to_rownames("Time")%>%t()%>%as_tibble()%>%
mutate(ID=ClimaticData$ID)%>%
pivot_longer(!ID,names_to="Date",values_to="Tmax")


TESTE<-do.call(cbind,TmaxChelsa)
colnames(TESTE)<-climatic_ord_date$values


TESTE[1,]%>%plot(,type="l")
TESTE%>%head()
TESTE%>%dim()
climatic_ord_date$values%>%length()
colnames(climatic_ord_date$values)
data.frame(
ClimaticData,
do.call(cbind,TmaxChelsa))%>%head()

tibble(ClimaticData$ID)

#==============================================================================
#  PARTIAL DATA EXTRACTED
#==============================================================================
#TMAX
#saveRDS(TmaxChelsa,file="MaxTemperatureChelsa_partial_27July23.rds")
# Problem with raster_file 283 . Completed until 1215

#TMIN
#saveRDS(TminChelsa,file="MinTemperatureChelsa_partial_27July23.rds")
# Problem with raster_file 83 . Completed until 1132

#saveRDS(PrecChelsa,file="PreciptationChelsa_partial_27July23.rds")
# No problems . Completed until 851
#===============================================================================
# FINISH DATA EXTRACTION
#===============================================================================
#-------------------------------------------------------------------------------
#		PRECIPTATION
#-------------------------------------------------------------------------------
PrecChelsa<-readRDS("./ChelsacrutsData/PreciptationChelsa_partial_01Jun23.rds")

for(i in length(PrecChelsa):length(raster_data_prec)){
X<-raster::raster(paste0("D:/Gabriel/prec/",raster_data_prec[i]))
PrecChelsa[[i]]<-raster::extract(X,
	ClimaticData[,c(3,2)],buffer=100, fun=mean)
gc()
print(i)}

#saveRDS(PrecChelsa,file="PreciptationChelsa_complete_01Jun23.rds")


#-------------------------------------------------------------------------------
#			TMIN
#-------------------------------------------------------------------------------
TminChelsa<-readRDS(file="./ChelsacrutsData/MinTemperatureChelsa_partial_01Jun23.rds")

for(i in length(TminChelsa):length(raster_data_tmin)){
X<-raster::raster(paste0("D:/Gabriel/Tmin/",raster_data_tmin[i]))
TminChelsa[[i]]<-raster::extract(X,
	ClimaticData[,c(3,2)],buffer=100, fun=mean)
gc()
print(i)}


unlist(lapply(TminChelsa,is.null))%>%length()
#saveRDS(TminChelsa,file="MinTemperatureChelsa_complete_04Jun23.rds")


#-------------------------------------------------------------------------------
#			TMAX
#-------------------------------------------------------------------------------
TmaxChelsa<-readRDS(file="./ChelsacrutsData/MaxTemperatureChelsa_partial_01Jun23.rds")

for(i in 1258:length(raster_data_tmax)){
X<-raster::raster(paste0("D:/Gabriel/Tmax/",raster_data_tmax[i]))
TmaxChelsa[[i]]<-raster::extract(X,
	ClimaticData[,c(3,2)],buffer=100, fun=mean)
gc()
print(i)}

#saveRDS(TmaxChelsa,file="MaxTemperatureChelsa_complete_04Jun23.rds")



#=================================================================
#	CHECK VALIDITY OF DATA EXTRACTION: PROBLEMS WITH BUFFERING
#=================================================================
AllAnimalsFitted<-MetadataClean%>%left_join(.,Metadata,by="ID")%>%
distinct(ID,.keep_all=TRUE)%>%
filter(Kingdom.x=="Animalia")%>%
select(ID,SpeciesAccepted.x,Lat,Lon)

TESTE<-climate_data_final%>%
filter(Prec<0)%>%
distinct(ID,.keep_all=TRUE)%>%
left_join(.,
Metadata%>%select(ID,Lat,Lon),
by="ID")%>%
distinct(ID,.keep_all=TRUE)

X<-terra::rast(file.choose())
plot(X)
points(AllAnimalsFitted%>%select(Lon,Lat),pch=21,col="blue",cex=1.1)
points(TESTE%>%select(Lon,Lat),pch=16,col="red")
