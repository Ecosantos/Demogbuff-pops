library(tidyverse)

# Load cleaned data
CleanData<-readRDS("Data/CleanData.RDS")
supertree<-readRDS("Data/supertree.RDS")

Metadata<-CleanData$Metadata
MetadataClean<-CleanData$MetadataClean

#Reduce data to improve redability
MedatadaFinal<-MetadataClean%>%select(-c(lambda,Ecoregion,Binomial))
MedatadaFinal<-MedatadaFinal%>%
  left_join(.,
            Metadata%>%select(ID,StudyStart, StudyDuration, StudyEnd)%>%distinct(),
            by="ID")

rm(CleanData)	#Remove non-used data to improve memory usage



TmaxChelsa<-readRDS(file="./Data/ChelsacrutsData/MaxTemperatureChelsa.rds")
TminChelsa<-readRDS(file="./Data/ChelsacrutsData/MinTemperatureChelsa.rds")
PrecipChelsa<-readRDS("./Data/ChelsacrutsData/PrecChelsa.rds")

# test
TmaxChelsa[1:10,1:10]%>%glimpse()
TminChelsa[1:10,1:10]%>%glimpse()
PrecipChelsa[1:10,1:10]%>%glimpse()

TmaxChelsa%>%dim()
TminChelsa%>%dim()
PrecipChelsa%>%dim()


TESTE<-readRDS("C:\\Artigos e resumos publicados submetidos ideias\\3 - Em desenvolvimento\\Demographic buffering continuum - Plants and animals\\Data and script\\Data\\ChelsacrutsData\\PreciptationChelsa_01Aug23.rds")


TmaxChelsa<-TmaxChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="TMax")%>%as_tibble()
TminChelsa<-TminChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="TMin")%>%as_tibble()
PrecipChelsa<-PrecipChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="Prec")%>%as_tibble()


summary(PrecipChelsa)
summary(TmaxChelsa)

dataclim<-left_join(TmaxChelsa,TminChelsa,by=c("param","Month","Year","ID"))%>%
  left_join(.,PrecipChelsa,by=c("param","Month","Year","ID"))%>%relocate(.,"ID",1)%>%
  mutate(
      Month=as.numeric(Month),
      Year=as.numeric(Year))

dataclim%>%glimpse()



ClimaticData<-Metadata%>%
  dplyr::select(ID,Lat,Lon,StudyStart,StudyEnd)%>%
  distinct(.)

dataclim%>%left_join(.,ClimaticData,by="ID")%>%glimpse()

climate_data_final<-ClimaticData%>%left_join(.,dataclim,by="ID")%>%
  filter(Year>=StudyStart & Year<=StudyEnd)%>%
  filter(complete.cases(.))%>%
  filter(TMax>-3000)%>%
  filter(Prec>-3000)%>%
  filter(TMin>-3000)%>%
  arrange(ID,Year,Month)


#============================================================================
# 				CREATING AN AUTOMATIC FUNCTION
#============================================================================
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
browser()
  }

climate_vars$ID2
tempset

tempset$TMax

plot(tempset$TMax,type="b")


unique(climate_vars$ID2)
climate_vars%>%do.call(data.frame,.)%>%head()

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
  tempout<-My_tsvars(tempset$TMax,freq=12)
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
  tempout<-My_tsvars(tempset$TMax,freq=12)
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

climate_df<-do.call(data.frame,climate_vars)%>%select(-ID2)

climate_df%>%head()
