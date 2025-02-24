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
comadre_flags<-comadre_flags%>%
  mutate(StudyID=cdb_id_studies(comadre_flags))%>%
  mutate(SpeciesAuthor2=str_replace_all(SpeciesAuthor, "_", ""))%>%
  group_by(SpeciesAuthor2,StudyID,Lat,Lon)%>%
  mutate(Pop=cur_group_id())%>%
  mutate(ID=paste0(abbreviate(SpeciesAuthor2,dot = TRUE),StudyID,"_",Pop))%>%
  ungroup()%>%select(-SpeciesAuthor2)

compadre_flags<-compadre_flags%>%
  mutate(StudyID=cdb_id_studies(compadre_flags))%>%
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
  # 	& MatrixComposite == "Individual" # REMOVED in 20/01/2025
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
  # 	& MatrixComposite == "Individual"  # REMOVED in 20/01/2025
  & StudyDuration >= 3
  & MatrixSplit == "Divided"
  & MatrixFec == "Yes"
  & MatrixCaptivity != "C"
  & SurvivalIssue<=1)