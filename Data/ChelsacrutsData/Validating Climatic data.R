# METADATA -------
# A temporary script to validate climatic data 
# extracted from Chelsa with Google Colab notebook


# Load data

TmaxChelsa<-readRDS(file="./Data/ChelsacrutsData/MaxTemperatureChelsa.rds")
TminChelsa<-readRDS(file="./Data/ChelsacrutsData/MinTemperatureChelsa.rds")
PrecipChelsa<-readRDS("./Data/ChelsacrutsData/PrecChelsa.rds")

# Check param, Month, Year is fine
TmaxChelsa[1:10,1:10]%>%glimpse()
TminChelsa[1:10,1:10]%>%glimpse()
PrecipChelsa[1:10,1:10]%>%glimpse()

# Check if dimensios are equal to all climatic dataset
TmaxChelsa%>%dim()
TminChelsa%>%dim()
PrecipChelsa%>%dim()


# Hist

TmaxChelsa<-TmaxChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="TMax")%>%as_tibble()
TminChelsa<-TminChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="TMin")%>%as_tibble()
PrecipChelsa<-PrecipChelsa%>%pivot_longer(!c(param,Month,Year),names_to="ID",values_to="Prec")%>%as_tibble()

range(PrecipChelsa$Prec)
range(TmaxChelsa$TMax)
range(TminChelsa$TMin)

hist(PrecipChelsa$Prec,xlim=c(0,3000),breaks=2000)
hist(TmaxChelsa$TMax,xlim=c(-400,500),breaks=2000,col=alpha("red",.3))
hist(TminChelsa$TMin,xlim=c(-400,500),breaks=2000,col=alpha("blue",.3),add=T)
