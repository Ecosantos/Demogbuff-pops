###########################################################
#	  	 LIFE HISTORY TRAITS CALCULATION 
# 			a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 		contact by ssantos.gabriel@gmail.com
#			    25 Aug 2023
###########################################################


#==================================================================
# 			Create mean MPMs 
#==================================================================
unique_ID_all<-unique_ID<-unique(Metadata$ID)
 MeanMPM<-Metadata%>%
		group_by(ID)%>%
		  group_map( ~ mpm_mean(.$mat))


#==================================================================
#!!!!!!!!  WARNING MENSAGE   !!!!!!!!!!!!!!!!!!!!!!
#"In mpm_mean(.$mat) :
# CompadreMat objects in given list do not all have the same MatrixClassAuthor. Returning MatrixClassAuthor from first list element
#------------------------------------------------------------------------
#Incongruences seems small. For instance check "MatrixClassAuthor"
#TESTE<-filter(Metadata,SpeciesAuthor=="Cimicifuga_rubifolia")%>%
#		  group_map( ~ mpm_mean(.$mat))
#TESTE
#lapply(filter(Metadata,SpeciesAuthor=="Cimicifuga_rubifolia")$mat,matrixClass)
#==================================================================
#------------------------------------------------------------------------
# Split MPMs in A, U, F, and C
#------------------------------------------------------------------------
meanA<-meanF<-meanU<-meanC<-NULL

   meanA<-matA(MeanMPM)
   meanU<-matU(MeanMPM) 
   meanF<-matF(MeanMPM)
   meanC<-matC(MeanMPM)

#==================================================================
# 	Check reproductive stages for each MPM
#==================================================================
L_repro_stages<-lapply(meanF,repro_stages)

#==================================================================
# 		Calculate life history traits
#==================================================================
#La = mean age at first reproduction
#LaProp = probability of achieving reproductive maturity
#e = mean life expectancy
#shape_surv(H)= Shape of survirvorship - Equivalent to Keyfitz' entropy (see Capdevila 2020 Func. Ecology https://doi.org/10.1111/1365-2435.13604
#shape_rep(D) = Shape of reproduction - Equivalent to Demetrius' entropy (see Baudisch & Stott 2019 -  https://doi.org/10.1111/2041-210X.13289
#growth = 
#------------------------------------------------------------------------

# first step - life history traits from MPMs
MPM_LHtraits<-n1<-NULL

for(i in 1:length(unique_ID)){
MPM_LHtraits$La[i]<-mature_age(matU = meanU[[i]], matR = meanF[[i]],start=1)
MPM_LHtraits$LaProb[i]<-mature_prob(matU = meanU[[i]], matR = meanF[[i]],start=1) 
MPM_LHtraits$e[i]<-life_expect_mean(matU = meanU[[i]], start = 1L)
n1[[i]]<-mature_distrib(meanU[[i]],L_repro_stages[[i]],start=1)
MPM_LHtraits$ID[i]<-unique_ID[i]
 #add verbose 
	if (i == 1 || i%%50 == 0) {
                message("Calculating Life history traits from mean MPM ", 
                  i)
	}
}

# second step - life history traits from life table
LifeTable_LHtraits<-NULL

for(i in 1:length(unique_ID)){
#j<-i
LifeTable_LHtraits$ID[i]<-MPM_LHtraits$ID[i]
n1[[i]]<-mature_distrib(meanU[[i]],L_repro_stages[[i]],start=1)
tryCatch({
LifeTable_LHtraits$D[i]<-shape_surv(meanU[[i]],trunc=TRUE)
LifeTable_LHtraits$Growth[i]<-vr_growth(meanU[[i]])
LifeTable_LHtraits$H[i]<-shape_rep(rep=meanF[[i]], surv=meanU[[i]])},
	error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
 #add verbose 
	if (i == 1 || i%%100 == 0) {
                message("Calculating Life history traits from Life table ", 
                  i)
	}
}

# Third step - Merge life history traits from MPMs and life tables
LHtraits<-left_join(by="ID",
data.frame(MPM_LHtraits),
data.frame(LifeTable_LHtraits))%>%
relocate(ID)%>%
as_tibble()

#----------------------------------------------------------------------------
#		DETECTING OUTLIERS USING PSYCH 
#----------------------------------------------------------------------------
PsychMD<-LHtraits%>%
as_tibble()%>%
filter(complete.cases(.))%>%
select(-c(ID))%>%
psych::outlier()

#Plot
LHtraits%>%
as_tibble()%>%
select(-c(ID))%>%
filter(complete.cases(.))%>%
mutate_if(is.numeric,scale)%>%
psych::pairs.panels(.,
#bg=c("blue","red")[(PsychMD>26)+1],pch=21,stars=T)
bg=c("blue","red")[(PsychMD>10)+1],pch=21,stars=T)


#----------------------------------------------------------------------------
#		REMOVING MULTIVARIATE OUTLIERS
#	Check outliers using Mahalanobis distance
# check here: https://stats.stackexchange.com/questions/203968/threshold-for-mahalanobis-distance#:~:text=What%20is%20a%20reasonable%20threshold,standard%20deviations%20would%20be%20reasonable.
#----------------------------------------------------------------------------
#Store outliers apart
LHtraits<-LHtraits%>%
filter(complete.cases(.))%>%
mutate_if(is.numeric,scale)%>%
column_to_rownames(var = "ID")%>%
mahalanobis_distance(.)%>%
rownames_to_column(var = "ID")

#==================================================================
#				EXPORT
#==================================================================
#Save LIFE HISTORY TRAITS
#saveRDS(LHtraits, "Data/LHtraits.RDS")



