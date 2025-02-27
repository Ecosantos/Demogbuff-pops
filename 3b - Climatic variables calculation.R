###########################################################
#	  	 CLIMATIC VARIABLES CALCULATION 
# 			a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 		contact by ssantos.gabriel@gmail.com
#			    25 Aug 2023
###########################################################

# 		CLIMATIC DATA 

#============================================================================
#	PREPARE CLIMATE dataset
#============================================================================
ClimaticData<-Metadata%>%
dplyr::select(ID,Lat,Lon,StudyStart,StudyEnd)%>%
distinct(.)

ordem_m<-c(10:12,1:9)

month_ord<-paste0("m",ordem_m)
year_ord<-paste0("y",1901:2016)

nyear<-year_ord%>%length()

yearMx<-data.frame(matrix(1,12,nyear))

dim(yearMx)[1]*dim(yearMx)[2]

colnames(yearMx)<-year_ord

climatic_ord_date<-data.frame(month_ord,yearMx)%>%
pivot_longer(!month_ord,names_to="year",values_to="values")

climatic_ord_date<-climatic_ord_date%>%
mutate(values=paste0(month_ord,year))


#============================================================================
#		Opening climatic data
#============================================================================
TmaxChelsa<-readRDS(file="./Data/ChelsacrutsData/MaxTemperatureChelsa.rds")
TminChelsa<-readRDS(file="./Data/ChelsacrutsData/MinTemperatureChelsa.rds")
PrecipChelsa<-readRDS("./Data/ChelsacrutsData/PrecChelsa.rds")

PrecipChelsa%>%glimpse()
TmaxChelsa

climate_data_final<-TmaxChelsa%>%
mutate(TMin=TminChelsa$TMin)%>%
mutate(Prec=PrecipChelsa$Prec)%>%
mutate(Time= str_replace(Time,"m",""))%>%
separate(.,Time, c("Month", "Year"),"y")%>%
filter(Year>=StudyStart & Year<=StudyEnd)%>%
filter(complete.cases(.))%>%
filter(TMax>-3000)%>%
filter(Prec>-3000)%>%
filter(TMin>-3000)%>%
mutate(Month=as.numeric(Month))%>%
arrange(ID,Year,Month)

climate_data_final%>%summary()

#========================================================================================
#		CALCULATING CLIMATIC VARIABLES
#========================================================================================
#Mean_trend = Mean value of trend
#Stoch_noisesize= STOCHASTICITY in the Time series; Magnitude of Stochasticity/Noise ; var(decompose(ts)$random)	#Here, I'll standardize the value using the Coefficient of variation!
#Stoch_rel = STOCHASTICITY in the Time series; (sd/mean)*100	#Coef. of variation of the timeseries
#Var_trend = STOCHASTICITY of the trend; (sd/mean)*100		#Coef. of variation of the trend	
#Ampli_season = (max(decompose(ts)$season) - min(decompose(ts)$season)/mean(trend)	#Proportional Amplitude 
#Ampli_trend = (max(decompose(ts)$trend) - min(decompose(ts)$trend))/mean(trend)	#Proportional Amplitude 

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
#climate_vars$Ampli_season_TMax[i]<-tempout$Ampli_season
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
#climate_vars$Mean_trend_TMin[i]<-tempout$Mean_trend			#Removed given the high collinearity with same information in TMax
#climate_vars$Stoch_noisesize_TMin[i]<-tempout$Stoch_noisesize	#Removed given the high collinearity with same information in TMax
#climate_vars$Ampli_season_TMin[i]<-tempout$Ampli_season
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
#climate_vars$Stoch_noisesize_Prec[i]<-tempout$Stoch_noisesize
#climate_vars$Ampli_season_Prec[i]<-tempout$Ampli_season
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

#Save climatic data 
#saveRDS(climate_df, paste0(DataDir,"/climate_df.RDS"))

