#'##########################################################
#	  	 CLIMATIC VARIABLES CALCULATION 
# 			a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 		contact by ssantos.gabriel@gmail.com
#			    25 Aug 2023
#'##########################################################

# 		CLIMATIC DATA -----
rm(list=ls())
library(tidyverse)

#'============================================================================
#		LOADING METADATA -------
#'============================================================================
CleanData<-readRDS("Data/CleanData.RDS")
Metadata<-CleanData$Metadata
MetadataClean<-CleanData$MetadataClean

#Reduce data to improve redability
MedatadaFinal<-MetadataClean%>%select(-c(lambda,Ecoregion,Binomial))
MedatadaFinal<-MedatadaFinal%>%
  left_join(.,
            Metadata%>%select(ID,StudyStart, StudyDuration, StudyEnd)%>%distinct(),
            by="ID")

rm(CleanData)	#Remove non-used data to improve memory usage



# Select only the necessary information to extract climatic data
ClimaticData<-Metadata%>%
  dplyr::select(ID,Lat,Lon,StudyStart,StudyEnd)%>%
  distinct(.)

#'============================================================================
#		Opening climatic data -------
#'============================================================================
TmaxChelsa<-readRDS(file="./Data/ChelsacrutsData/MaxTemperatureChelsa.rds")
TminChelsa<-readRDS(file="./Data/ChelsacrutsData/MinTemperatureChelsa.rds")
PrecipChelsa<-readRDS("./Data/ChelsacrutsData/PrecChelsa.rds")

# Convert to long format - 
TmaxChelsa<-TmaxChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="TMax")%>%as_tibble()
TminChelsa<-TminChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="TMin")%>%as_tibble()
PrecipChelsa<-PrecipChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="Prec")%>%as_tibble()

# Merge multiple climatic data in a single dataset
dataclim<-left_join(TmaxChelsa,TminChelsa,by=c("param","Month","Year","ID"))%>%
  left_join(.,PrecipChelsa,by=c("param","Month","Year","ID"))%>%relocate(.,"ID",1)%>%
  mutate(
    Month=as.numeric(Month),
    Year=as.numeric(Year))



climate_data_final<-ClimaticData%>%left_join(.,dataclim,by="ID")%>%
  filter(Year>=StudyStart & Year<=StudyEnd)%>%
  filter(complete.cases(.))%>%
  filter(TMax>-3000)%>%
  filter(Prec>-3000)%>%
  filter(TMin>-3000)%>%
  arrange(ID,Year,Month)

climate_data_final%>%glimpse()

#'========================================================================================
#		CALCULATING CLIMATIC VARIABLES   ------
#'========================================================================================
# Mean_trend = Mean value of trend
#Stoch_noisesize= STOCHASTICITY in the Time series; Magnitude of Stochasticity/Noise ; var(decompose(ts)$random)	#Here, I'll standardize the value using the Coefficient of variation!
#Stoch_rel = STOCHASTICITY in the Time series; (sd/mean)*100	#Coef. of variation of the timeseries
#Var_trend = STOCHASTICITY of the trend; (sd/mean)*100		#Coef. of variation of the trend	
#Ampli_season = (max(decompose(ts)$season) - min(decompose(ts)$season)/mean(trend)	#Proportional Amplitude 
#Ampli_trend = (max(decompose(ts)$trend) - min(decompose(ts)$trend))/mean(trend)	#Proportional Amplitude 

#'============================================================================
# 				CREATING AN AUTOMATIC FUNCTION
#             My_tsvars ------
#'============================================================================
My_tsvars<-function(x,freq=freq){
  out<-decomp<-tempts<-NULL
  tempts<-ts(x,frequency=freq)
  if(!any(is.na(tempts))){
    print("imputeTS NOT used")
    decomp<-decompose(tempts)
  }
  else{
    print("imputeTS used")
    decomp<-decompose(imputeTS::na_interpolation(tempts, option = "spline"))
  } 
  out$Mean_trend<-mean(na.omit(decomp$trend))
  out$Stoch_noisesize<-sd(na.omit(decomp$random))
  out$Stoch_rel<-(sd(na.omit(decomp$x))/mean(na.omit(decomp$x)))*100
  out$Var_trend<-(sd(na.omit(decomp$trend))/mean(na.omit(decomp$trend)))*100
  out$Ampli_season<-abs((max(na.omit(decomp$seasonal))-min(na.omit(decomp$seasonal)))/mean(na.omit(decomp$trend)))
  out$Ampli_trend<-abs(max(na.omit(decomp$trend))-min(na.omit(decomp$trend)))/mean(na.omit(decomp$trend))
  return(out)}

#Usage
#My_tsvars(filter(climate_data_final,ID=="Acrs.330")$TMax,freq=12)
#-----------------------------------------------------------------------------------------

tempout<-tempset<-climate_vars<-NULL
climate_vars$ID<-unique(climate_data_final$ID)

#TEMP MAXIMUM
for(i in 1:length(climate_vars$ID)){
  climate_vars$ID2[i]<-climate_vars$ID[i]	#ID2 is created as a quality check. At the end ID must be ID=ID2
  tempset<-subset(climate_data_final,ID==climate_vars$ID2[i])
  #Jump timeseries that doesnt fit the condition of two years complete
  #if(length(ts(tempset$TMax))< 2* 12) next		#Check if there is at least two years complete. This line make sure we are working with two years that might not hold in case when there is missing data in the last December - This is probably a temporary issue
  tempout<-My_tsvars(tempset$TMax,freq=12)
  climate_vars$Mean_trend_TMax[i]<-tempout$Mean_trend
  climate_vars$Stoch_noisesize_TMax[i]<-tempout$Stoch_noisesize#
  climate_vars$Ampli_season_TMax[i]<-tempout$Ampli_season
  climate_vars$Ampli_trend_TMax[i]<-tempout$Ampli_trend			#Removed given high collinearity given the amplitude of seasons
  #add verbose 
  if (i == 1 || i%%25 == 0) {
    message("Calculating mean Matrices", 
            i)
  }
}


unique(climate_vars$ID2)
climate_vars%>%do.call(data.frame,.)%>%head()

#TEMP MINIMUM
for(i in 1:length(climate_vars$ID)){
  climate_vars$ID2[i]<-climate_vars$ID[i]
  tempset<-subset(climate_data_final,ID==climate_vars$ID2[i])
  #Jump timeseries that doesnt fit the condition of two years complete
  #if(length(ts(tempset$TMax))< 2* 12) next		#Check if there is at least two years complete. This line make sure we are working with two years that might not hold in case when there is missing data in the last December - This is probably a temporary issue
  tempout<-My_tsvars(tempset$TMin,freq=12)
  climate_vars$Mean_trend_TMin[i]<-tempout$Mean_trend			#Removed given the high collinearity with same information in TMax
  climate_vars$Stoch_noisesize_TMin[i]<-tempout$Stoch_noisesize	#Removed given the high collinearity with same information in TMax
  climate_vars$Ampli_season_TMin[i]<-tempout$Ampli_season
  climate_vars$Ampli_trend_TMin[i]<-tempout$Ampli_trend		#Removed given high collinearity given the amplitude of seasons
  #add verbose 
  if (i == 1 || i%%25 == 0) {
    message("Calculating mean Matrices", 
            i)
  }
}

#PRECIPTATION
for(i in 1:length(climate_vars$ID)){
  climate_vars$ID2[i]<-climate_vars$ID[i]
  tempset<-subset(climate_data_final,ID==climate_vars$ID2[i])
  #Jump timeseries that doesnt fit the condition of two years complete
  #if(length(ts(tempset$TMax))< 2* 12) next		#Check if there is at least two years complete. This line make sure we are working with two years that might not hold in case when there is missing data in the last December - This is probably a temporary issue
  tempout<-My_tsvars(tempset$Prec,freq=12)
  climate_vars$Mean_trend_Prec[i]<-tempout$Mean_trend
  climate_vars$Stoch_noisesize_Prec[i]<-tempout$Stoch_noisesize
  climate_vars$Ampli_season_Prec[i]<-tempout$Ampli_season
  climate_vars$Ampli_trend_Prec[i]<-tempout$Ampli_trend
  #add verbose 
  if (i == 1 || i%%25 == 0) {
    message("Calculating mean Matrices", 
            i)
  }
}

#========================================================================================
#		MERGING CLIMATIC VARIABLES AND EXPORT
#========================================================================================

climate_df<-do.call(data.frame,climate_vars)%>%select(-ID2)

climate_df%>%head()
climate_df%>%glimpse()
climate_df[,-1]%>%cor()%>%corrplot::corrplot()

#Save climatic data 
#saveRDS(climate_df, "Data/climate_df.RDS")

