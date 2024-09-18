#=====================================================================================================================================
#		LOWER-LEVEL (vital rate) STOCHASTIC ELASTICITY WITH RESPECT OF MEAN AND VARIANCE
#=====================================================================================================================================
#Function created in agreement with 
	#vitalRatePerturbation function available in Rage package
		#https://github.com/jonesor/Rage
	#stoch.sens - previous in this code

#--------------------------------------------------------------------------------------------------------------------------------------
# Function: array_to_matrix
# ancillary functio  to transform matrices in arrays
#--------------------------------------------------------------------------------------------------------------------------------------
array_to_matrix<-function(A){
lapply(seq(dim(A)[3]), function(x) A[ , , x])}

print("An ancillary function - 'array_to_matrix' - has been created, be pround!")

#--------------------------------------------------------------------------------------------------------------------------------------
# Function: my.vitalRatePerturbation
# Core function to calculate sensitivity, elasticity and sum of stochasticy elasticities
#--------------------------------------------------------------------------------------------------------------------------------------
my.vitalRatePerturbation<-function(matUs,matFs,matCs){
# Initial Inputs 
#Calculate interatively the stochastic sensitivity for matA, matF, matU and matC 
   k <- ncol(matUs[[1]])
   tlimit=100  
   wvec <- rep(1/k, k)
   w <- cbind(wvec)
   vvec <- rep(1/k, k)
   v <- cbind(vvec)
   matDim <- nrow(matUs[[1]])
 if (is.null(matCs)) {
        matCs <- array_to_matrix(replicate(length(matUs),(matrix(0, matDim, matDim))))
    }

matUs<-rep(matUs,tlimit,length.out=tlimit)
matFs<-rep(matFs,tlimit,length.out=tlimit)
matCs<-rep(matCs,tlimit,length.out=tlimit)

IDS<-rep(paste0("t",(1:length(matUs))),tlimit,length.out=tlimit)
IDS<-data.frame(IDS=IDS,order=1:tlimit)
IDS<-IDS[sample(IDS$order),]

matUs<-matUs[IDS$order]
matFs<-matFs[IDS$order]
matCs<-matCs[IDS$order]

matAs<-NULL
for(i in 1: length(matUs)){
matAs[[i]]<-matUs[[i]] + matFs[[i]] + matCs[[i]]
     }

r <- rep(0, tlimit)
    for (i in 1:tlimit) {
        a <- matAs[[i]]
        wvec <- a %*% wvec
        r[i] <- sum(wvec)
        wvec <- wvec/r[i]
        w <- cbind(w, wvec)
    }

 for (i in rev(1:tlimit)) {
        a <- matAs[[i]]
        vvec <- vvec %*% a
        v <- cbind(t(vvec), v)
    }

elasmean<-elasmat <-sensmat <- matrix(0, nrow = matDim, ncol = matDim)
Meansensmat<-Sigsensmat<-NULL
for (i in 1:tlimit) {
      sensmat <- sensmat + ((v[, i + 1] %*% t(w[, i]))/as.numeric(r[i] * t(v[, i + 1]) %*% w[, i + 1]))
Sigsensmat[[i]]<- (sensmat-sensitivity(mean(matAs)))/tlimit
Meansensmat[[i]]<- (sensmat+Sigsensmat[[i]])/tlimit

# Stochastic elasticity from Matrix elements 
elasmat <- elasmat + ((v[, i + 1] %*% t(w[, i]) * a)/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
elasmean <- elasmean + ((v[, i + 1] %*% t(w[, i]) * mean(matAs))/as.numeric((r[i] * 
            t(v[, i + 1]) %*% w[, i + 1])))
}

sensA<-sensmat <- sensmat/tlimit
Sigsensmat<-mean(Sigsensmat)
Meansensmat<-mean(Meansensmat)
elasmat <- elasmat/tlimit
elasmean <- elasmean/tlimit

sigmas <- lapply(matUs,colSums)
sigMed <- apply(do.call(rbind,sigmas),2,mean)
#==================================================================================
#Vital rate segragation and stochastic contribution to lambda
# Each vital rate is separeted by a block code
#==================================================================================
sigsensSurv<-sigelasSurv<-matsigSensSurv<-elasSurv<-sensSurv <-matSensSurv <-noSurvA <- noSurvA_F <- NULL
meansensSurv<-meanelasSurv<-matmeanSensSurv<-NULL
#---------------------------------------------------------------------------------
#SURVIVAL
#---------------------------------------------------------------------------------
for(i in 1:tlimit){
 noSurvA[[i]] <- noSurvA_F[[i]] <- t(t(matAs[[i]])/sigmas[[i]])
 noSurvA[[i]][, sigmas[[i]] == 0] <- matAs[[i]][, sigmas[[i]] == 0]
 matSensSurv[[i]] <- noSurvA[[i]] * sensA
 matsigSensSurv[[i]] <- noSurvA[[i]] * Sigsensmat
  matmeanSensSurv[[i]] <- noSurvA[[i]] * Meansensmat
sensSurv[[i]]<-colSums(matSensSurv[[i]])
 sigsensSurv[[i]]<-colSums(matsigSensSurv[[i]])
 meansensSurv[[i]]<-colSums(matmeanSensSurv[[i]])
 sensSurv[[i]][sigmas[[i]] == 0 ]<-0
 elasSurv[[i]]<-sigmas[[i]] * sensSurv[[i]]/r[[i]]
 meanelasSurv[[i]]<-sigMed * sensSurv[[i]]/r[[i]]
 sigelasSurv[[i]]<-(sigmas[[i]]-sigMed) * sensSurv[[i]]/r[[i]]
}

#---------------------------------------------------------------------------------
# GROWTH & SHRINKING
#---------------------------------------------------------------------------------
   stage_ageFUN<- function(x){ apply(x, 2, function(x) length(which(x > 0))) == 1}
   stage_agedef<- lapply(matUs, stage_ageFUN)

   lwr <- upr <- matrix(0, nrow = matDim, ncol = matDim)
   lwr[lower.tri(lwr, diag = TRUE)] <- 1
lwr<-array_to_matrix(replicate( tlimit,lwr))
    upr[upper.tri(upr, diag = TRUE)] <- 1;
upr<-array_to_matrix(replicate( tlimit,upr))
matmeanSensGrowShri<-matsigSensGrowShri<- matSensGrowShri <- array_to_matrix(replicate( tlimit,matrix(NA_real_, matDim, matDim)))

matsigElasGrow <- matElasGrow <- NULL
matmeanElasShri<-matmeanElasGrow<-matmeanElasGrowShri <- matsigElasShri<-matsigElasGrowShri <- matElasShri<-matElasGrowShri <- NULL
meanelasGrow<-meanelasShri<-NULL
sigelasShri<- sigelasGrow <- elasShri<- elasGrow <- NULL
meansensShri <- sigsensShri<- meansensGrow<-sigsensGrow<- sensShri<- sensGrow<- NULL
matelasSensShri<-matelasSensGrow <-NULL
matsigSensShri<-matsigSensGrow <-matmeanSensGrow<-matmeanSensShri<-matSensShri<-matSensGrow <-NULL

for(i in 1:tlimit){
#middle step
    for (j in 1:matDim) {
        matSensGrowShri[[i]][, j] <- sensA[j, j] * (-sigmas[[i]][j]) + sensA[,j] * (sigmas[[i]][j])
  matmeanSensGrowShri[[i]][, j] <- sensA[j, j] * (-sigMed[j]) + sensA[,j] * (sigMed[j])
  matsigSensGrowShri[[i]][, j] <- sensA[j, j] * (-(sigmas[[i]]-sigMed)[j]) + sensA[,j] * ((sigmas[[i]]-sigMed)[j])
 }
 # ~~~~ GROWTH ~~~~ 
matSensGrowShri[[i]][which(matUs[[i]] == 0)] <- 0
matsigSensGrowShri[[i]][which(matUs[[i]] == 0)] <- 0
matmeanSensGrowShri[[i]][which(matUs[[i]] == 0)] <- 0
 matSensGrowShri[[i]][, stage_agedef[[i]] == TRUE] <- 0
 matsigSensGrowShri[[i]][, stage_agedef[[i]] == TRUE] <- 0
 matmeanSensGrowShri[[i]][, stage_agedef[[i]] == TRUE] <- 0
  matSensGrow[[i]] <- lwr[[i]] * matSensGrowShri[[i]]
   matsigSensGrow[[i]] <- lwr[[i]] * matsigSensGrowShri[[i]]
   matmeanSensGrow[[i]] <- lwr[[i]] * matmeanSensGrowShri[[i]]
   sensGrow[[i]] <- colSums(matSensGrow[[i]], na.rm = TRUE)
    sigsensGrow[[i]] <- colSums(matsigSensGrow[[i]], na.rm = TRUE)
    meansensGrow[[i]] <- colSums(matmeanSensGrow[[i]], na.rm = TRUE)
matSensGrow[[i]] <- lwr[[i]] * matSensGrowShri[[i]]
matsigSensGrow[[i]] <- lwr[[i]] * matsigSensGrowShri[[i]]
matmeanSensGrow[[i]] <- lwr[[i]] * matmeanSensGrowShri[[i]]
    sensGrow[[i]] <- colSums(matSensGrow[[i]], na.rm = TRUE)
     sigsensGrow[[i]] <- colSums(matsigSensGrow[[i]], na.rm = TRUE)
     meansensGrow[[i]] <- colSums(matmeanSensGrow[[i]], na.rm = TRUE)
 matElasGrowShri[[i]] <- noSurvA [[i]] * matSensGrowShri[[i]]/r[[i]]
 matsigElasGrowShri[[i]] <- noSurvA [[i]] * matsigSensGrowShri[[i]]/r[[i]]
 matmeanElasGrowShri[[i]] <- noSurvA [[i]] * matmeanSensGrowShri[[i]]/r[[i]]
    matElasGrow[[i]] <- lwr[[i]] * matElasGrowShri[[i]]
    matsigElasGrow[[i]] <- lwr[[i]] * matsigElasGrowShri[[i]]
    matmeanElasGrow[[i]] <- lwr[[i]] * matmeanElasGrowShri[[i]]

elasGrow[[i]] <- colSums(matElasGrow[[i]], na.rm = TRUE)
sigelasGrow[[i]] <- colSums(matsigElasGrow[[i]], na.rm = TRUE)
 meanelasGrow[[i]] <- colSums(matmeanElasGrow[[i]], na.rm = TRUE)
    elasGrow[[i]][sigmas[[i]] == 0] <- 0
 sigelasGrow[[i]][sigmas[[i]] == 0] <- 0
 meanelasGrow[[i]][sigmas[[i]] == 0] <- 0

# ~~~~ SHRINKING ~~~~ 
 matSensShri[[i]] <- upr[[i]]  * matSensGrowShri[[i]]
 matsigSensShri[[i]] <- upr[[i]]  * matsigSensGrowShri[[i]]
 matmeanSensShri[[i]] <- upr[[i]]  * matmeanSensGrowShri[[i]]
  sensShri[[i]]  <- colSums(matSensShri[[i]] , na.rm = TRUE)
   sigsensShri[[i]]  <- colSums(matsigSensShri[[i]] , na.rm = TRUE)
   meansensShri[[i]]  <- colSums(matmeanSensShri[[i]] , na.rm = TRUE)
matElasShri[[i]] <- upr[[i]] * matElasGrowShri[[i]]
 matsigElasShri[[i]] <- upr[[i]] * matsigElasGrowShri[[i]]
 matmeanElasShri[[i]] <- upr[[i]] * matmeanElasGrowShri[[i]]
   elasShri[[i]] <- colSums(matElasShri[[i]], na.rm = TRUE)
   sigelasShri[[i]] <- colSums(matsigElasShri[[i]], na.rm = TRUE)
   meanelasShri[[i]] <- colSums(matmeanElasShri[[i]], na.rm = TRUE)
     elasShri[[i]][sigmas[[i]] == 0] <- 0
     sigelasShri[[i]][sigmas[[i]] == 0] <- 0
     meanelasShri[[i]][sigmas[[i]] == 0] <- 0
	 }  #Growth and Shrinking end!

#---------------------------------------------------------------------------------
### FECUNDITY & CLONALITY 
#---------------------------------------------------------------------------------
   sigma_rep <- sigmas
   sigmaMed_rep <- apply(do.call(rbind,sigma_rep),2,mean)
   
for(i in 1:tlimit){ sigma_rep[[i]][sigma_rep[[i]] == 0] <- 1}
matSensClo <- matSensFec <- array_to_matrix(replicate( tlimit,sensA))
matsigSensClo <- matsigSensFec <- array_to_matrix(replicate( tlimit,Sigsensmat))
matmeanSensClo <- matmeanSensFec <- array_to_matrix(replicate( tlimit,Meansensmat))

elasClo <- matElasClo <- elasFec <- matElasFec <- sensClo<-sensFec<-NULL
sigelasClo <- matsigElasClo <- sigelasFec <- matsigElasFec <- sigsensClo<-sigsensFec<-NULL
meanelasClo <- matmeanElasClo <- meanelasFec <- matmeanElasFec <- meansensClo<-meansensFec<-NULL
for(i in 1:tlimit){ 

 # ~~~ FECUNDITY ~~~~
 matSensFec[[i]][which(matFs[[i]] == 0)] <- 0
 matsigSensFec[[i]][which(matFs[[i]] == 0)] <- 0
 matmeanSensFec[[i]][which(matFs[[i]] == 0)] <- 0
   matSensFec[[i]] <- t(t(matSensFec[[i]]) * sigma_rep[[i]])
   matmeanSensFec[[i]] <- t(t(matmeanSensFec[[i]]) * sigmaMed_rep)
   matsigSensFec[[i]] <- t(t(matsigSensFec[[i]]) * (sigma_rep[[i]]-sigmaMed_rep))
     sensFec[[i]] <- colSums(matSensFec[[i]])
     sigsensFec[[i]] <- colSums(matsigSensFec[[i]])
     meansensFec[[i]] <- colSums(matmeanSensFec[[i]])
 matElasFec[[i]] <- noSurvA[[i]] * matSensFec[[i]]/r[[i]]
 matsigElasFec[[i]] <- noSurvA[[i]] * matsigSensFec[[i]]/r[[i]]
 matmeanElasFec[[i]] <- noSurvA[[i]] * matmeanSensFec[[i]]/r[[i]]
    elasFec[[i]] <- colSums(matElasFec[[i]])
    sigelasFec[[i]] <- colSums(matsigElasFec[[i]])
    meanelasFec[[i]] <- colSums(matmeanElasFec[[i]])
 
 # ~~~~ CLONALITY ~~~~
matSensClo[[i]][which(matCs[[i]] == 0)] <- 0
matsigSensClo[[i]][which(matCs[[i]] == 0)] <- 0
matmeanSensClo[[i]][which(matCs[[i]] == 0)] <- 0
    matSensClo[[i]] <- t(t(matSensClo[[i]]) * sigma_rep[[i]])
    matmeanSensClo[[i]] <- t(t(matSensClo[[i]]) * sigmaMed_rep)
    matsigSensClo[[i]] <- t(t(matSensClo[[i]]) * (sigma_rep[[i]]-sigmaMed_rep))
    sensClo[[i]] <- colSums(matSensClo[[i]])
    sigsensClo[[i]] <- colSums(matsigSensClo[[i]])
    meansensClo[[i]] <- colSums(matmeanSensClo[[i]])
matElasClo[[i]] <- noSurvA[[i]] * matSensClo[[i]]/r[[i]]
matsigElasClo[[i]] <- noSurvA[[i]] * matsigSensClo[[i]]/r[[i]]
matmeanElasClo[[i]] <- noSurvA[[i]] * matmeanSensClo[[i]]/r[[i]]
    elasClo[[i]] <- colSums(matElasClo[[i]])
    sigelasClo[[i]] <- colSums(matsigElasClo[[i]])
    meanelasClo[[i]] <- colSums(matmeanElasClo[[i]])
}#Reproduction end!

#---------------------------------------------------------------------------------
#Lower-level Stochastic Sensitivities
#----------------------------------------------------------------------
## Mean
Sens_mean=c(mean(unlist(lapply(sensFec,sum))),
mean(unlist(lapply(sensGrow,sum))),
mean(unlist(lapply(sensShri,sum))),
mean(unlist(lapply(sensClo,sum))),
mean(unlist(lapply(sensSurv,sum))))
## SD
Sens_sd=c(sd(unlist(lapply(sensFec,sum))),
sd(unlist(lapply(sensGrow,sum))),
sd(unlist(lapply(sensShri,sum))),
sd(unlist(lapply(sensClo,sum))),
sd(unlist(lapply(sensSurv,sum))))
#----------------------------------------------------------------------
#Lower-level Stochastic Elasticities
#----------------------------------------------------------------------
## Mean
Elas_mean=c(mean(unlist(lapply(elasFec,sum))),
mean(unlist(lapply(elasGrow,sum))),
mean(unlist(lapply(elasShri,sum))),
mean(unlist(lapply(elasClo,sum))),
mean(unlist(lapply(elasSurv,sum))))
## SD
Elas_sd=c(sd(unlist(lapply(elasFec,sum))),
sd(unlist(lapply(elasGrow,sum))),
sd(unlist(lapply(elasShri,sum))),
sd(unlist(lapply(elasClo,sum))),
sd(unlist(lapply(elasSurv,sum))))

#----------------------------------------------------------------------
#Lower-level Stochastic in respect to the mean
#----------------------------------------------------------------------
## Mean
Elas_stochmean_mean=c(mean(unlist(lapply(meanelasFec,sum))),
mean(unlist(lapply(meanelasGrow,sum))),
mean(unlist(lapply(meanelasShri,sum))),
mean(unlist(lapply(meanelasClo,sum))),
mean(unlist(lapply(meanelasSurv,sum))))

## SD
Elas_stochmean_sd=c(sd(unlist(lapply(meanelasFec,sum))),
sd(unlist(lapply(meanelasGrow,sum))),
sd(unlist(lapply(meanelasShri,sum))),
sd(unlist(lapply(meanelasClo,sum))),
sd(unlist(lapply(meanelasSurv,sum))))

#----------------------------------------------------------------------
# Lower-level Elasticities in respect to the variance
#----------------------------------------------------------------------
##Mean
Elas_sig_mean=c(mean(unlist(lapply(sigelasFec,sum))),
mean(unlist(lapply(sigelasGrow,sum))),
mean(unlist(lapply(sigelasShri,sum))),
mean(unlist(lapply(sigelasClo,sum))),
mean(unlist(lapply(sigelasSurv,sum))))

##SD
Elas_sig_sd=c(sd(unlist(lapply(sigelasFec,sum))),
sd(unlist(lapply(sigelasGrow,sum))),
sd(unlist(lapply(sigelasShri,sum))),
sd(unlist(lapply(sigelasClo,sum))),
sd(unlist(lapply(sigelasSurv,sum))))

#----------------------------------------------------------------------
# ADD CUMULATIVE mean and sd values for Lower-level Elasticities in respect to the mean and variance
#----------------------------------------------------------------------
## Mean - in respect to the MEAN
  Elas_stochmean_mean_SUM=mean(
	unlist(lapply(meanelasFec,sum))+
	unlist(lapply(meanelasGrow,sum))+
	unlist(lapply(meanelasShri,sum))+
	unlist(lapply(meanelasClo,sum))+
	unlist(lapply(meanelasSurv,sum)))
Elas_stochmean_mean=c(Elas_stochmean_mean,Elas_stochmean_mean_SUM)

#-----

## sd - in respect to the MEAN
Elas_stochmean_sd_SUM=sd(
	unlist(lapply(meanelasFec,sum))+
	unlist(lapply(meanelasGrow,sum))+
	unlist(lapply(meanelasShri,sum))+
	unlist(lapply(meanelasClo,sum))+
	unlist(lapply(meanelasSurv,sum)))
Elas_stochmean_sd=c(Elas_stochmean_sd,Elas_stochmean_sd_SUM)

#-----

## Mean - in respect to the VARIANCE
  Elas_sig_mean_SUM=mean(
	unlist(lapply(sigelasFec,sum))+
	unlist(lapply(sigelasGrow,sum))+
	unlist(lapply(sigelasShri,sum))+
	unlist(lapply(sigelasClo,sum))+
	unlist(lapply(sigelasSurv,sum)))

Elas_sig_mean=c(Elas_sig_mean,Elas_sig_mean_SUM)


## sd - in respect to the VARIANCE
  Elas_sig_sd_SUM=sd(
	unlist(lapply(sigelasFec,sum))+
	unlist(lapply(sigelasGrow,sum))+
	unlist(lapply(sigelasShri,sum))+
	unlist(lapply(sigelasClo,sum))+
	unlist(lapply(sigelasSurv,sum)))

Elas_sig_sd=c(Elas_sig_sd,Elas_sig_sd_SUM)


# FUNCTION OUTPUT

 label<-c("Reproduction","Growth","Shrinking","Clonality","Survival")
 Sens=data.frame(Mean= Sens_mean, SD=Sens_sd,row.names=label)
 Elas=data.frame(Mean= Elas_mean, SD=Elas_sd,row.names=label)
 Elas_stochmean=data.frame(Mean= Elas_stochmean_mean, SD=Elas_stochmean_sd,row.names=c(label,"Cumulative"))
 Elas_sigma=data.frame(Mean= Elas_sig_mean, SD=Elas_sig_sd,row.names=c(label,"Cumulative"))
out<-list(Sens,Elas,Elas_stochmean,Elas_sigma,elasmean)
names(out)<-c(
"Lower-level Stochastic Sensitivities",
 "Lower-level Stochastic Elasticities",
 "Lower-level Stochastic Elasticity with respect to the mean", 
 "Lower-level Stochastic Elasticity with respect to the Variance",
 "ElasMean")
return(out)
}

#===============================================================================================

print("Core function 'my.vitalRatePerturbation' - has been created, You nailed it!")
print(cat(paste0("'my.vitalRatePerturbation' is based on: \n",
"- VitalRatePerturbation function in the Rage package: https://github.com/jonesor/Rage. \n",
"- Stochastic elasticticity in Caswell 2001 p.406-407. \n",
"Outputs are: \n",
"[[1]] Lower-level Stochastic Sensitivities \n",
"[[2]] Lower-level Stochastic Elasticities \n",
"[[3]] Lower-level Stochastic Elasticities \n",
"[[4]] Lower-level Stochastic Elasticities with respect  to variance\n",
"[[5]] Matrix level stochastic Elasticitity with respect to the mean")))

