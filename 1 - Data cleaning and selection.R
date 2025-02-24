###########################################################
#		Data and Metadata selection and cleaning	
# 			a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 		contact by ssantos.gabriel@gmail.com
#			    25 Aug 2023
###########################################################

#==================================================================
#		LOAD & FILTERING
#==================================================================

mosaic <- Rmosaic::mos_fetch("v1.0.0")	

#mos_fetch for some reason is bugged so, it clean everything else, so, run again the directory
DataDir<-"Data"


load(paste0(DataDir,"/COMADRE_v.4.23.3.1.RData"))
load(paste0(DataDir,"/COMPADRE_v.6.23.5.0.RData"))

#---------------------------------------------
#	Check red flags
# definition in: https://jonesor.github.io/Rcompadre/reference/cdb_flag.html
#---------------------------------------------
compadre<-compadre <- as_cdb(compadre)
comadre <- as_cdb(comadre)

compadre_flags <- cdb_flag(compadre)
comadre_flags <- cdb_flag(comadre)

#---------------------------------------------
#Add ID for each study
#---------------------------------------------
compadre_flags<-compadre_flags%>%
mutate(StudyID=cdb_id_studies(compadre_flags))%>%
mutate(SpeciesAuthor2=str_replace_all(SpeciesAuthor, "_", ""))%>%
group_by(SpeciesAuthor2,StudyID,Lat,Lon)%>%
mutate(Pop=cur_group_id())%>%
mutate(ID=paste0(abbreviate(SpeciesAuthor2,dot = TRUE),StudyID,"_",Pop))%>%
ungroup()%>%select(-SpeciesAuthor2)

comadre_flags<-comadre_flags%>%
  mutate(StudyID=cdb_id_studies(comadre_flags))%>%
  mutate(SpeciesAuthor2=str_replace_all(SpeciesAuthor, "_", ""))%>%
  group_by(SpeciesAuthor2,StudyID,Lat,Lon)%>%
  mutate(Pop=cur_group_id())%>%
  mutate(ID=paste0(abbreviate(SpeciesAuthor2,dot = TRUE),StudyID,"_",Pop))%>%
  ungroup()%>%select(-SpeciesAuthor2)


#---------------------------------------------
# Subset relevant (adequated) data 
#---------------------------------------------
compadre_sub <- subset(
  compadre_flags,
 	check_NA_A == FALSE & check_NA_U==FALSE & check_NA_F==FALSE
	& check_zero_F == FALSE & check_zero_U==FALSE 
	& check_singular_U == FALSE 
	& check_component_sum == TRUE 
	& check_ergodic == TRUE 
 	& MatrixComposite == "Individual" 
& StudyDuration >= 3
& MatrixSplit == "Divided"
& MatrixFec == "Yes"
& MatrixCaptivity != "C"
& SurvivalIssue<=1)

comadre_sub <- subset(
  comadre_flags,
 	check_NA_A == FALSE & check_NA_U==FALSE & check_NA_F==FALSE
	& check_zero_F == FALSE & check_zero_U==FALSE 
	& check_singular_U == FALSE 
	& check_component_sum == TRUE 
	& check_ergodic == TRUE 
 	& MatrixComposite == "Individual"  # REMOVED in 20/01/2025
& StudyDuration >= 3
& MatrixSplit == "Divided"
& MatrixFec == "Yes"
& MatrixCaptivity != "C"
& SurvivalIssue<=1)

#---------------------------------------------
# Subset columns to improve human readability
#---------------------------------------------
compadre_sub%>%names()

#Columns to keep in metadata
VarKeep<-c("mat","ID","Pop","SpeciesAuthor","StudyStart","StudyDuration","StudyEnd",
"Species","Genus","Family","CommonName","Phylum","Order",
"OrganismType","AngioGymno","DicotMonoc","ProjectionInterval","MatrixCriteriaOntogeny",
"MatrixPopulation","Lat","Altitude","Continent",
"StudiedSex","MatrixSeasonal","MatrixCaptivity","MatrixStartSeason","MatrixEndYear",
"MatrixEndMonth","MatrixSplit","Observations","StudyID","MatrixID",
"SpeciesAccepted","Kingdom","Class","MatrixCriteriaSize","MatrixCriteriaAge",
"NumberPopulations","Lon","Country","Ecoregion","MatrixComposite","MatrixTreatment",
"MatrixStartYear","MatrixStartMonth","MatrixEndSeason","CensusType",
"MatrixFec","MatrixDimension")


compadre_sub%>%select(all_of(VarKeep))
compadre_traits<-compadre_sub%>%select(all_of(VarKeep))
comadre_traits<-comadre_sub%>%select(all_of(VarKeep))

#---------------------------------------------
#	Filtering unrealistic lambdas
#---------------------------------------------
Metadata<-bind_rows(compadre_traits%>%as_tibble(),comadre_traits%>%as_tibble())%>%
mutate(lambda=sapply(matA(mat), popdemo::eigs, what = "lambda"))%>%
filter(lambda < 2)%>%
filter(lambda > 0.5)

#Check
Metadata%>%ggplot(.,aes(x=lambda))+geom_histogram()

#============================================================================
#	REMOVE DATA WITHOUT PHYLOGENETIC CORRESPONDENCE
#============================================================================
Metadata<-Metadata%>%mutate(Binomial=paste0(Genus,"_",Species))

Metadata<-Metadata%>%
filter(Binomial%in%c(mosaic@phylogeny$animalPhylogeny$tip,mosaic@phylogeny$plantPhylogeny$tip))

MetadataClean<-Metadata%>%select(ID,SpeciesAccepted,Kingdom,Phylum,Class,Order,OrganismType,AngioGymno,lambda,Ecoregion,Binomial)

MetadataClean


#============================================================================
#	MERGE PHYLOGENIES IN A SINGLE SUPER TREE
#============================================================================

supertree<-ape::bind.tree(mosaic@phylogeny$animalPhylogeny,
	mosaic@phylogeny$plantPhylogeny,
		where = "root", position = -.10^10, interactive = FALSE)



#============================================================================
#				EXPORT
#============================================================================

#Save Metadata and Metadata Clean
Metadatas = list(Metadata= Metadata, MetadataClean=MetadataClean)
#saveRDS(Metadatas , "Data/CleanData.RDS")

#Save mosaic rethrived data 
#saveRDS(mosaic, "Data/mosaicdata.RDS")

#Save mosaic phylogenies merged
#saveRDS(supertree, "Data/supertree.RDS")




