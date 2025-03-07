#'##########################################################
#	DEMOGRAPHIC BUFFERING CONTINUUM 
#		 plants and animals
#		 by Gabriel Santos
#	contact by ssantos.gabriel@gmail.com
#		 	16 September 2024
#'##########################################################

set.seed(1)

rm(list=ls())

#Packages
library(tidyverse)	#v2.0.0
library(popbio)		#v2.7		
library(popdemo)		#v1.3.1
library(ggplot2)		#v3.5.0
library(scales)		#v1.3.0
library(tidyr)		#v1.3.0
library(viridis)		#v0.6.3
library("FactoMineR")	#v2.8
library("factoextra")	#v1.0.7
library(Rcompadre)	#v1.2.1
library(Rage)		#v1.4.0
library(vegan)		#v2.6.4
library(rstatix)		#0.7.2
#library(Rmosaic) ; remotes::install_github("mosaicdatabase/Rmosaic")
#library(imputeTS)	
library(tidybayes)	#v3.0.6
library("treeio")		#v1.24.1
library("ggtree")		
library(MCMCglmm)		#v2.35
#phytools			#v1.9.16 - In use but not loaded


rm(list=ls())


#'==================================================================
#	0. SessionInfo and package versions ------
#'==================================================================
#Check package versions
cbind(unlist(loadedNamespaces()),
  unlist(lapply(
	lapply(loadedNamespaces(),packageVersion),
		as.character)))%>%
		data.frame()%>%arrange(-desc(X1))

#'==================================================================
#		1.	README -----
#'==================================================================
# The following script provide the general framework used to analyse
# demographic buffering continuum in Santos et al. in review:
# Population responses to environmental stochasticity are primarily driven by survival-reproduction trade-offs and mediated by aridity
#
# Intermediary code and data steps are sourced along the framework
# Intermediary code and data include:
#	1 - Data cleaning and selection and its intermediary data
#	2 - Life history traits calculation and its intermediary data
#	3 - Climatic variables calculation
#	4 - Analyses with MCMCglmm 
#	5 - Core function: Stochastic elasticities of variance


#'==================================================================
#		COMPADRE, COMADRE and MOSAIC 
#	  2. DATA SELECTION AND CLEANING ------
#'------------------------------------------------------------------
# Script avaliable in "1 - Data cleaning and selection.R"
#	file.edit("1 - Data cleaning and selection.R")
# Produce two datasets:
# 	- CleanData: 
#		Filter matrix singularity, presence of fecundity, individual matrices only
#	- supertree:
#		Extract supertree from MOSAIC database (Bernard et al. 2023 Scientific Data).
#'==================================================================

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

#'==================================================================
#   3. LIFE HISTORY TRAITS     -------
#		Calculate life-history traits
#'------------------------------------------------------------------
# Script avaliable in "2 - Life history traits calculation.R"
#	file.edit("2 - Life history traits calculation.R")
# Produce dataset LHtraits.RDS:
# 	- Calculate life history traits and detect outliers
#'==================================================================
LHtraits<-readRDS("Data/LHtraits.RDS")
LHtraits
## .. Filter outliers ----
spLHmat<-LHtraits%>%
	filter(is.outlier=="FALSE")%>%
		select(-c(is.outlier,mahal.dist))	%>%
column_to_rownames(var = "ID")

#Correlation plot - check colinearity
#spLHmat%>%cor()%>%corrplot::corrplot(.,title="Check Colinerarity on LH traits")

#'============================================================================
#   3.2	FAST-SLOW CONTINUUM - PCA -----
#'============================================================================
LHpca<-spLHmat%>%filter(complete.cases(.))%>%PCA(.,graph=FALSE,scale=TRUE) # Using raw data without imputatation

# LHpca - Map variables
LHpca$eig%>%t()	#Explained variables
LHpca$ind$coord	#Eigenvalues
dimdesc(LHpca, axes = 1:3, proba = 0.05)

# 3.2.2. Check PCA -----
cowplot::plot_grid(nrow=1,
cowplot::plot_grid(ncol=1,
fviz_eig(LHpca)+theme_bw(base_size=14),
fviz_pca(LHpca, geom=c("point"))+theme_bw(base_size=14)),
LHpca$var$cor%>%ggcorrplot::ggcorrplot(method = "circle")+theme_bw(base_size=14)+
  theme(axis.text.x=element_text(angle=75,hjust=1)))+
coord_fixed(ratio=.50)

#ggsave(file="Figures/PCALH.svg")


#Use Principal Components as variables - add to LHtraits
LHtraits<-left_join(by="ID",
	LHtraits,
	LHpca$ind$coord%>%data.frame()%>%
	rownames_to_column(var="ID")%>%as_tibble()%>%
rename_with(., ~ gsub("Dim", "LHAxis", .x, fixed = TRUE))%>%
	select(ID,LHAxis.1:LHAxis.2))%>%
filter(complete.cases(.))

#'============================================================================
#   4. CLIMATIC VARIABLES and ENVIRONMENTAL PCA ----
#		Extract and summarise climatic information
#'----------------------------------------------------------------------------
# Extracted climatic data is available in three different files:
#  1. MinTemperaturaChelsa.rds
#  2. MaxTemperaturaChelsa.rds
#  3. PrecChelsa.rds
## Because data download and information extraction is quite time consuming:
# A dedicated Jupyter notebook to run on Google Colab is provided in:
#     - "ChelsaData Download and extraction - Google Colab.ipynb"
#     - Run this Google Colab by typing in your brownser:
#     - https://githubtocolab.com/Ecosantos/Demogbuff-pops
# 
# 3 - Climatic variables calculation.R use the already extracted climatic data to:
#   - Data summary and metrics extraction described in methods
#	  - Produce final dataset to analyse: "climate_df.RDS"
#	  - Environmental trend, amplitude, and stochasticisticy
#'============================================================================
climate_df<-readRDS("Data/climate_df.RDS")


## 4.1. Checking collinearity -------
cor_matrix <- climate_df %>% select(-ID) %>% cor()
diag(cor_matrix) <- NA  

# Find values higher than 0.65
high_corr <- which(abs(cor_matrix) > 0.65, arr.ind = TRUE)

# Output to selection
data.frame(
  Var1 = rownames(cor_matrix)[high_corr[, 1]],
  Var2 = colnames(cor_matrix)[high_corr[, 2]],
  Correlation = cor_matrix[high_corr]
) %>% arrange(desc(abs(Correlation)))  

# Plot all variables to visualize collinearity
#climate_df%>%select(-ID)%>%
#  select(-c(Ampli_season_TMax,
#            Ampli_season_TMin,
#            Ampli_season_Prec))%>%
#  cor()%>%corrplot::corrplot()

## 4.2 Removing collinar climatic variables -----
climate_df<-climate_df%>%select(
  -c(Mean_trend_TMin,
  Stoch_noisesize_TMin,Stoch_noisesize_Prec,
  Ampli_season_TMax,Ampli_season_TMin,Ampli_season_Prec))%>%
  as_tibble()


# Produce the climatic/environmental PCA
ClimPCA<-climate_df%>%
mutate_at(vars(-c(ID)),scale)%>%
column_to_rownames("ID")%>%
#cor(.)%>%corrplot::corrplot(.)
PCA(.,graph=F,scale=FALSE)


cowplot::plot_grid(nrow=1,
cowplot::plot_grid(ncol=1,
fviz_eig(ClimPCA)+theme_bw(base_size=14),
fviz_pca(ClimPCA, geom=c("point"))+theme_bw(base_size=14)),
ClimPCA$var$cor%>%ggcorrplot::ggcorrplot(method = "circle")+theme_bw(base_size=14)+
  theme(axis.text.x=element_text(angle=75,hjust=1)))+
coord_fixed(ratio=.60)


#ggsave(file="Figures/PCAEnv.svg")

dimdesc(ClimPCA, axes = 1:3, proba = 0.05)


#'------------------------------------------------------------------------------
##	4.3	Mapping climatic axes ----
#'------------------------------------------------------------------------------
Plot_clim<-Metadata%>%dplyr::select(.,c(ID,Lat,Lon))%>%
distinct(.,.keep_all=T)%>%
left_join(.,
data.frame(ClimPCA$ind$coord)%>%rownames_to_column(.,var="ID"),
by="ID")%>%
filter(complete.cases(.))%>%
mutate(across(Dim.1:Dim.3,scales::rescale,to=c(0,1)))

#Plot map
world <- map_data("world")

gmap<-ggplot() +
geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "grey50", fill = "grey95", linewidth= 0.05)

gmap +
 geom_point(
    data = left_join(Plot_clim,distinct(MetadataClean,ID,.keep_all=T)%>%select(ID,Kingdom)),
    aes(x=Lon, y=Lat, fill = Dim.1,stroke=.3,size=Dim.2,shape=Kingdom))+
scale_shape_manual(values=c(23,21))+
scale_fill_viridis(option="magma")+
scale_size(range = c(1.5, 5))+
 theme_void()

#ggsave(file="Figures/ClimMap_all_taxa2.svg")

#'============================================================================
# 	5. CALCULATE DEMOGRAPHIC BUFFERING FROM LOWER LEVEL VITAL RATES -----
#
#  Source mainfunction used in the calculation of 
#	stochastic elasticities within respect to variance of lower level vital rates
#	MainFunction - Stochastic elasticities of variance lower level.R
# 	  - This script creates the following functions:
#		1. "my.vitalRatePerturbation": Produce the stochasticity elasticity with respect to variance.
#		2. "array_to_matrix": Ancilliary function
#'-----------------------------------------------------------------------------

source("MainFunction - Stochastic elasticities of variance lower level.R")

#'---------------------------------------------------------------------------------------
# EXAMPLE: Case where shrinking is possible in Animal
# Alternation between reproductive non-reproductive stage
# Reproductive estage <--> non-reproductive stage
#'---------------------------------------------------------------------------------------
filter(Metadata,ID=="Urs2.183_726")$mat
lapply(filter(Metadata,ID=="Urs2.183_726")$mat,matA)

my.vitalRatePerturbation(
lapply(filter(Metadata,ID=="Urs2.183_726")$mat,matU),
lapply(filter(Metadata,ID=="Urs2.183_726")$mat,matF),
lapply(filter(Metadata,ID=="Urs2.183_726")$mat,matC))
#'---------------------------------------------------------------------------------------

# Organize unique ID to rethrive matrices to metadata after analyses
uniqueID<-unique(Metadata$ID)

# Create vectors to accomated estimated buffering and ancilliary statistics
Buffmx<-temp<-MatRep<-StochElasVR<-ElasSigVR<-SigRatioVR<-NULL


for(i in 1:length(uniqueID)){
temp<-my.vitalRatePerturbation(	
lapply(filter(Metadata,ID==uniqueID[i])$mat,matU),	#Determines matU
lapply(filter(Metadata,ID==uniqueID[i])$mat,matF),	#Determines matF
lapply(filter(Metadata,ID==uniqueID[i])$mat,matC))	#Determines matC
Buffmx[[i]]<-sum(temp[[5]])
StochElasVR[[i]]<-temp[[3]]
ElasSigVR[[i]]<-temp[[4]]
names(StochElasVR)[[i]]<-uniqueID[i]
names(ElasSigVR)[[i]]<-uniqueID[i]
# quantify timeseries length 
MatRep[[i]]<-length(lapply(filter(Metadata,ID==uniqueID[i])$mat,matU))
 #add verbose 
	if (i == 1 || i%%40 == 0) {
                message("Calculating mean Matrices", 
                  i)
	}
rm(temp)
}

# Quantify the number of populations with timeseries longer than three
unlist(MatRep)[unlist(MatRep)>2]%>%length()

#'---------------------------------------------------------------------------
#		Merge Demographic buffering calculations in a single Data frame
#'---------------------------------------------------------------------------
ElasSigVR_full<-lapply(ElasSigVR,rownames_to_column, var = "VR")%>%
Map(cbind, ID = names(.), .)%>% do.call(rbind,.)%>% as_tibble()%>%
pivot_wider(names_from = VR,values_from=c(Mean,SD))%>%
filter(complete.cases(.))

#'----------
#A simpler version is created without standard desviation
# as standard desviation was not included in the model
# The full version ("ElasSigVR_full") is still required to summary statistics below
#	see RANKING BUFFERING
#'----------

ElasSigVR_data<-lapply(ElasSigVR,rownames_to_column, var = "VR")%>%
Map(cbind, ID = names(.), .)%>% do.call(rbind,.)%>% as_tibble()%>%
select(-SD)%>%				#Remove standard desviation
pivot_wider(names_from = VR,values_from=Mean)%>%
filter(complete.cases(.))


StochElasVR_data<-lapply(StochElasVR,rownames_to_column, var = "VR")%>%
Map(cbind, ID = names(.), .)%>%
do.call(rbind,.)%>% as_tibble()%>% select(-SD)%>%
pivot_wider(names_from = VR,values_from=Mean)%>%
filter(complete.cases(.))



#'---------------------------------------------------------------------------
#		Final demographic buffering data preparation. 
# Then, it will be ready for merging with climatic data and life history traits
#'---------------------------------------------------------------------------
databuff_data<-left_join(
ElasSigVR_data%>%setNames(paste0(names(.),'_SigElas')),
StochElasVR_data%>%setNames(paste0(names(.),'_Base')),
by=c("ID_SigElas"="ID_Base"))

MPMinfo<-data.frame(
 Buffmx=unlist(Buffmx),
  MatRep=unlist(MatRep),
   ID=uniqueID)%>%as_tibble()


databuff_data%>%glimpse()


databuff_data<-left_join(
 databuff_data,MPMinfo,
	by=c("ID_SigElas"="ID"))%>%
		filter(MatRep>2)%>%
		  rename(ID="ID_SigElas")

#Comparing buffering in vital rates and life history traits
databuff_data%>%
left_join(.,LHtraits,by="ID")%>%
filter(complete.cases(.))%>%
select(-"ID")%>%cor()%>%corrplot::corrplot()

#'========================================================================================
#	---- 6. MERGING CLIMATIC + LIFE HISTORY + BUFFERING DATA -----
#'========================================================================================
colnames(ClimPCA$ind$coord)<-gsub("Dim.", "ClimPC.", colnames(ClimPCA$ind$coord))
colnames(LHpca$ind$coord)<-gsub("Dim.", "LHPC.", colnames(LHpca$ind$coord))

LHpca$ind$coord%>%head()
LHpca_axes12<-LHpca$ind$coord[,c(1,2)]%>%data.frame()%>%rownames_to_column("ID")
ClimPCA_axes123<-ClimPCA$ind$coord[,c(1,2,3)]%>%data.frame()%>%rownames_to_column("ID")

ClimPCA_axes123%>%glimpse()

merged_data<-databuff_data%>%
left_join(.,LHpca_axes12,by="ID")%>%
left_join(.,ClimPCA_axes123,by="ID")%>%
filter(complete.cases(.))%>%
left_join(.,MetadataClean,by="ID")%>%
distinct(ID,.keep_all=TRUE)


merged_data%>%glimpse()

#Correlation plot
merged_data%>%select_if(is.numeric)%>%cor()%>%corrplot::corrplot()


#'========================================================================================
#	  ---- 7. THE SUPER TREE ---- 
# BUILT IT AND MAKE SURE IT WORKS ON PHYLOGENETIC ANALYSES
#'========================================================================================
#Make subtree
sppINphylo<-unique(merged_data$Binomial)[unique(merged_data$Binomial)%in%supertree$tip]

subtree<-keep.tip(supertree, sppINphylo)

#Check and Avoid duplicates
any(duplicated(subtree$node.label))
subtree<-makeNodeLabel(subtree); any(duplicated(subtree$node.label))

#Avoid polytomyes
subtree_backup<-subtree<-multi2di(subtree)


#Make sure branches are comparable, non-negative and non-zero
subtree$edge.length<-scales::rescale(subtree$edge.length,to=c(0.00001,.99999))

#Force ultrametric
subtree<-phytools::force.ultrametric(subtree, method="extend")

#Check if structure are kept
identical(subtree,subtree_backup)
all.equal.phylo(subtree,subtree_backup,use.edge.length=FALSE)

is.rooted(subtree)
is.binary(subtree)
is.ultrametric(subtree)
any(subtree$edge.length==0)
subtree$edge.length[subtree$edge.length==0]

#'================================================================================
#		----- 8. FINAL DATASETS AND PHYLOGENETIC ANALYSES -----
#'================================================================================

final_data<-merged_data%>%
mutate(phylo=Binomial)%>%
mutate(inPhylo=Binomial%in%sppINphylo)%>%
filter(inPhylo!=FALSE)%>%
select(-inPhylo)%>%
arrange(.,match(Binomial,subtree$tip))%>%
data.frame()

#Total populations
final_data%>%dim()

#Total species
final_data$SpeciesAccepted%>%unique()%>%length()

#Populations by kingdom
final_data%>%select(Kingdom)%>%table()

#Species by kingdom
final_data%>%select(SpeciesAccepted,Kingdom)%>%
distinct()%>%select(Kingdom)%>%table()

# USED MPMs
filter(Metadata, ID %in% final_data$ID)$mat%>%length()

#'=======================================================================================
## -------	8.1 DATASET PREPARATION AND GLMM ANALYSES -----
#'=======================================================================================
#Prepare dataset for phylogenetic analyses
final_data%>%
  as_tibble()%>%
  distinct(final_data$Binomial,.keep_all=TRUE)$Binomial%>%
  column_to_rownames("Binomial")%>%
  select(is.numeric,Kingdom)%>%
  select(-c(Clonality_Base,Clonality_SigElas))

#Prepare phylogeny with respective populations 
subtree_Animals<-keep.tip(subtree, filter(final_data,Kingdom=="Animalia")$phylo)
subtree_Plants<-keep.tip(subtree, filter(final_data,Kingdom=="Plantae")$phylo)

# Create a new data only to run MCMCglmm
data_model<-final_data%>%select(-c(Reproduction_Base:Cumulative_Base))

#Transform CUMULATIVE IN ABSOLUTE VALUE - VERY IMPORTANT! 
data_model$Cumulative_SigElas<-abs(data_model$Cumulative_SigElas)


#save(data_model,subtree_Animals,subtree_Plants,
#     file = "Data/GLMMdata.Rdata")

#'--------------------------------------------------------------------------------------
###	---------- 8.2.1 GLMM parameteres --------------
#'--------------------------------------------------------------------------------------
# Determines the fixed effect component
fixEffect<-fixEffect<-"~LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.2 * ClimPC.3"

# Determines all variables of interest to make multiple models
InterestingVars<-c("Survival","Growth","Shrinking","Reproduction","Clonality","Buffmx","Cumulative")


traits<-traits_glmm<- unique (grep(paste(InterestingVars,collapse="|"), 
                                   colnames(data_model), value=TRUE))

#'--------------------------------------------------------------------------------------
##    ---- 8.2.2 Export GLMM DATA & accessory info ---- 
##	Export GLMM data and accessory information to run externally if necessary
# I prefer to run with Google Colab Notebook to improve efficiency
#'--------------------------------------------------------------------------------------
# For model selection use
# file.edit("5b - MCMCglmm.R")

# Once best model is reached use
# file.edit("5b - MCMCglmm.R")


#'=======================================================================================
#		---- 9. BUFFERING PATTERNS ---- 
#'=======================================================================================
final_data%>%glimpse()

final_data_meta<-final_data%>%
left_join(.,MetadataClean%>%select(-lambda)%>%distinct(),by="ID")

final_data_meta%>%glimpse()

#'-----------------------------------------------------------------------------------
##	---- 9.1 RANKING BUFFERING ---- 
#'-----------------------------------------------------------------------------------
ToRank<-final_data%>%
mutate(Cumulative_SigElas=abs(Cumulative_SigElas))%>%
group_by(Kingdom)%>%
select(ID,Cumulative_SigElas,SpeciesAccepted,AngioGymno,Class)%>%
  mutate(rank  = rank(Cumulative_SigElas, ties.method = "random"))

ToRank%>%glimpse()

ToRank<-ToRank%>%
group_by(Kingdom)%>%
mutate(min=min(Cumulative_SigElas),
		max=max(Cumulative_SigElas))%>%
filter(rank==min(rank)|rank==max(rank))%>%
select(-c(min,max))%>%
left_join(.,Metadata%>%select(ID, CommonName,Family)%>%distinct(),by="ID")%>%
select(2,3,7,4,5,6,8,9)


# Merging
ToRank%>%
left_join(.,ElasSigVR_full,by="ID")%>%mutate_if(is.numeric,round,3)%>%as.data.frame()

#'-----------------------------------------------------------------------------------
##	----- 9.2. Check overlapping Animals x plants 	----- 
#'-----------------------------------------------------------------------------------

t.test(abs(final_data$Cumulative_SigElas)~final_data$Kingdom)

#'-----------------------------------------------------------------------------------
##		----- 9.3. PROPORTIONAL CONTRIBUTION 	----- 
#'-----------------------------------------------------------------------------------
RelContrib<-final_data%>%as_tibble()%>%
select(ID,Reproduction_SigElas:Survival_SigElas)%>%
mutate_if(is.numeric,abs)%>%
mutate(NewCumulative=rowSums(across(where(is.numeric))))%>%
#mutate_if(is.numeric,round,4)%>%glimpse()
pivot_longer(!c(ID,NewCumulative))%>%
mutate(RelativeProp=(value/NewCumulative)*100)%>%
left_join(.,MetadataClean%>%select(-lambda)%>%distinct(),by="ID")


### 9.3.1. Relative contribution by organismType -------
RelContrib%>%
group_by(Kingdom,name,OrganismType)%>%
summarise(Relativemean=mean(RelativeProp),n=n())%>%
pivot_wider(names_from=name,values_from = Relativemean)%>%
group_by(Kingdom)%>%group_split()

### 9.3.2. Relative contribution by kingdon -------
RelContribtax<-RelContrib%>%
group_by(name,Kingdom)%>%
summarise(Relativemean=mean(RelativeProp),
		SD=sd(RelativeProp),
		n=n())


RelContribtax%>%
  group_by(Kingdom)%>%group_split()


#'======================================================================================================
#	---- 10. LOADING GLMM OUTPUTS AND PROCESS RESULTS   -----
#'======================================================================================================

MCMCglmm_output<-readRDS("Data/MCMCglmm_output.rds")

lapply(MCMCglmm_output,names)

## Phylogenetic models
MCMCglmm_phylo_plants<-MCMCglmm_output$Phylogenetic_models[[1]]
MCMCglmm_phylo_animals<-MCMCglmm_output$Phylogenetic_models[[2]]

## Non-Phylogenetic (simple) models
MCMCglmm_simple_plants<-MCMCglmm_output$Simple_models[[1]]
MCMCglmm_simple_animals<-MCMCglmm_output$Simple_models[[2]]

rm(MCMCglmm_output)

MCMCglmm_simple_plants


#'------------------------------------------------------------------------------------------------------
##	---- 10.1 DATA HARMONIZATION 	---- 
#'------------------------------------------------------------------------------------------------------
#'------------------------------------------------------------------------------------------------------
### Acessory function
#'------------------------------------------------------------------------------------------------------
process_models<-function(X,Taxa=Taxa,Model=Model){
  # Verify if arguments match existing models
  if (!Taxa %in% c("Plants", "Animals") || !Model %in% c("Phylo", "Simple")) {
    stop("Please use Taxa = 'Plants' or 'Animals' and Model = 'Phylo' or 'Simple'") }
  summ<-lapply(X,function(inner_list) summary(inner_list))
  coff<-lapply(summ,function(inner_list) inner_list[[5]])
  coff_df<-lapply(coff,as.data.frame)
  coff_df2<-lapply(Phylo_models_plants_coefs,rownames_to_column, var = "Statistics")
  do.call(rbind,Map(cbind, Trait = names(coff_df2),Taxa=Taxa, Model=Model, coff_df2))
}
#'------------------------------------------------------------------------------------------------------

# Phylogenetic corrected models
Phylo_models_df_plants<-process_models(MCMCglmm_phylo_plants,Taxa="Plants", Model="Phylo")
Phylo_models_df_animals<-process_models(MCMCglmm_phylo_animals,Taxa="Animals", Model="Phylo")

# Simple models
Simple_models_df_plants<-process_models(MCMCglmm_phylo_plants,Taxa="Plants", Model="Simple")
Simple_models_df_animals<-process_models(MCMCglmm_phylo_animals,Taxa="Animals", Model="Simple")

#'------------------------------------------------------------------------------------------------------

GLMMs_df<-rbind(Phylo_models_df_plants,Phylo_models_df_animals,
                Simple_models_df_plants,Simple_models_df_animals)


colnames(GLMMs_df)<-c("Trait","Taxa","Model","Statistics","post.mean","low95","high95","eff.samp","pMCMC")

#'---------------------------------------------------------------------------------------------------
## -----  10.1.1. SUMMARY SINTHESIS  -----  
#'---------------------------------------------------------------------------------------------------

GLMMs_df_summary<-GLMMs_df%>%
  mutate(sig=ifelse(pMCMC<=0.05,"Sig","Non-Sig"))%>%
  filter(pMCMC<=0.1)%>%
  filter(Statistics!="(Intercept)")


# PLANTs x Cumulative only
GLMMs_df_summary%>%
  filter(Taxa=="Plants" & Trait == "Cumulative_SigElas")%>%
  group_by(Trait)


# By vital rates - Plants & Animals
GLMMs_df_summary%>%
  filter(Trait != "Cumulative_SigElas")%>%
  group_by(Trait)%>%group_split()


#'======================================================================================================
#			----- 10.1.2. Posterior distributions -------
#'======================================================================================================
Phylo_posterior_animals<-lapply(MCMCglmm_phylo_animals,function(inner_list) data.frame(inner_list$Sol,Taxa="Animals",Model="Phylo"))
Phylo_posterior_plants<-lapply(MCMCglmm_phylo_plants,function(inner_list) data.frame(inner_list$Sol,Taxa="Plants",Model="Phylo"))

Simple_posterior_animals<-lapply(MCMCglmm_simple_animals,function(inner_list) data.frame(inner_list$Sol,Taxa="Animals",Model="Simple"))
Simple_posterior_plants<-lapply(MCMCglmm_simple_plants,function(inner_list) data.frame(inner_list$Sol,Taxa="Plants",Model="Simple"))

Posterior_data<-rbind(
  do.call(rbind,Phylo_posterior_animals),
  do.call(rbind,Phylo_posterior_plants),
  do.call(rbind,Simple_posterior_animals),
  do.call(rbind,Simple_posterior_plants))%>%
  rownames_to_column(., var = "VAR")%>%
  separate(VAR,c("Trait"))%>%
  pivot_longer(!c(Trait,Taxa,Model),values_to="Values",names_to="Variables")%>%
  mutate(Variables=ifelse(Variables=="LHPC.1.LHPC.2","LHPC.1:LHPC.2",Variables))%>%
  mutate(Variables=ifelse(Variables=="ClimPC.1.ClimPC.2","ClimPC.1:ClimPC.2",Variables))%>%
  as_tibble()%>%
  filter(Variables!="X.Intercept.")

Posterior_data<-left_join(Posterior_data,
                          GLMMs_df%>%separate(Trait,"Trait"),
                          by=c("Trait","Taxa","Model","Variables"="Statistics"))%>%
  mutate(sig=ifelse(pMCMC<=0.05,"Sig","Non-Sig"))


#==========================================================================
#	---- 11. Estimating PHYLOGENETIC SIGNAL ---- 
#==========================================================================

# Create an axilliary function
my.fake.lamb<-function(model){
  out =  model$VCV[,"phylo"]/
    (model$VCV[,"phylo"]+
       #            model$VCV[,"species"]+
       model$VCV[,"units"])
  mean.Lambda=mean(out)
  SE.lambda=(sd(out)/sqrt(length(out)))
  return(out)
}



#InterestingVars<-c("Survival","Growth","Shrinking","Reproduction","Clonality","Buffmx","Cumulative")
traits<-c("Reproduction_SigElas", "Growth_SigElas", "Shrinking_SigElas",  "Clonality_SigElas", "Survival_SigElas",
          "Cumulative_SigElas","Buffmx")  

lapply(MCMCglmm_phylo_animals,function(inner_list) inner_list$VCV)

Phylo_signal_animals<-lapply(MCMCglmm_phylo_animals,function(inner_list) as.data.frame(my.fake.lamb(inner_list)))%>%
  do.call(cbind,.)

Phylo_signal_plants<-lapply(MCMCglmm_phylo_plants,function(inner_list) as.data.frame(my.fake.lamb(inner_list)))%>%
  do.call(cbind,.)

colnames(Phylo_signal_animals)<-str_split_i(traits, "_", 1)
colnames(Phylo_signal_plants)<-str_split_i(traits, "_", 1)

Phylo_signal_df<-rbind(
  data.frame(Phylo_signal_plants,Taxa="Plants"),
  data.frame(Phylo_signal_animals,Taxa="Animals"))%>%
  pivot_longer(!Taxa,names_to="Trait",values_to="Values")




#'====================================================================
# -----   12. H2: Phylogenetic signal x E_sig -----
#		TEST IF PHYLOGENETIC SIGNAL IS GREATER 
#   IN VITAL RATES THAT CONTRIBUTE MOST
#'====================================================================

RelContribtax

Phylo_summary # Originally from 5c - Tidy MCMCglmm output

# TESTING Relative contribution
cor.test_df<-RelContribtax %>%
  filter(Kingdom=="Plantae")%>%
  mutate(name = gsub("_SigElas", "", name, fixed = TRUE))%>%
  rename(Trait = name, E_sig = Relativemean)%>%
  left_join(., Phylo_summary, by = "Trait")


  cor.test(cor.test_df$E_sig,cor.test_df$MEDIAN,method="kendall")
  
plot(cor.test_df$E_sig~cor.test_df$MEDIAN)
  
# TESTING raw values - Marginally significant! 
cor.test_df<-final_data%>%as_tibble()%>%
  select(ID,Reproduction_SigElas:Survival_SigElas)%>%
#  mutate_if(is.numeric,abs)%>%
pivot_longer(!c(ID))%>%
  left_join(.,MetadataClean%>%select(ID,Kingdom)%>%distinct(),by="ID")%>%
  group_by(name,Kingdom)%>%
  summarise(E_sig=mean(value),
            SD=sd(value),
            n=n())%>%
  mutate(name = gsub("_SigElas", "", name, fixed = TRUE))%>%
  rename(Trait = name)%>%
  left_join(., Phylo_summary, by = "Trait")%>%
  filter(Kingdom=="Plantae")

cor.test_df
cor.test(cor.test_df$E_sig,cor.test_df$MEDIAN,method="spearman")

plot(cor.test_df$E_sig~cor.test_df$MEDIAN)


#'====================================================================
#		----- 11. FIGURES ------
#'====================================================================

#'====================================================================
##		11.1. FIGURE 1 -----
#'====================================================================
#'------------------------------------------------------------------------
#	ACCESSORY FUNCTION
#'------------------------------------------------------------------------

scientific_10 <- function(x) { 
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#'------------------------------------------------------------------------

# Calculate stochastic Growth Rate (Stoch Lambs)
temp<-Stochslambs<-NULL
for(i in 1:length(unique(final_data$ID))){
temp<-filter(
	Metadata,ID==unique(final_data$ID)[i])$mat%>%
		matA()%>%
		stoch.growth.rate(.,maxt = 200,verbose = FALSE)
Stochslambs$approx[[i]]<-temp$approx
Stochslambs$sim[[i]]<-temp$sim
Stochslambs$CIlow[[i]]<-temp$sim.CI[[1]]
Stochslambs$CIhigh[[i]]<-temp$sim.CI[[2]]
#add verbose 
	if (i == 1 || i%%40 == 0) {
              message("Calculating mean Matrices", 
                 i)
	}
rm(temp)
	}

Stochslambs_df<-data.frame(
	unique(final_data$ID),
		lapply(Stochslambs,cbind))%>%
	as_tibble()%>% unnest()%>%
		mutate_if(is.numeric,exp)%>%
	rename(., ID = unique.final_data.ID.)

# Plot
final_data%>%
 select(ID:Cumulative_SigElas,Buffmx:MatRep,OrganismType,Kingdom,lambda)%>%
 pivot_longer(!c(ID,MatRep:lambda))%>%
 mutate(values=ifelse(value==0,rnorm(1,1e-7,1e-7),value))%>%
 left_join(.,Stochslambs_df,by='ID')%>%
 filter(name=="Cumulative_SigElas")%>%
 #filter(!(value>quantile(value,.95)))%>%
 #filter(!(value<quantile(value,.025)))%>%
 ggplot(.,aes(y=approx,x=abs(value)))+
# geom_point(aes(size=sqrt(MatRep),fill=Kingdom),shape=21)
#geom_point(aes(size=sqrt(MatRep),fill=Kingdom),shape=21)+
geom_pointrange(aes(ymin=CIlow,ymax=CIhigh,fill=Kingdom,shape=Kingdom),size=.9)+
scale_fill_manual(values=c("#144375","#1f9a59"))+
scale_shape_manual(values=c(23,21))+
geom_hline(yintercept=1,linetype=2,color="grey50",size=.9)+
scale_x_continuous(trans ='log',
breaks=c(0,0.000001,0.0001,0.001,0.001,0.01,0.1,1),
label = scientific_10)+
#annotation_logticks() + 
#xlab(expression("(More Buffered)            -       " %<-% "        "~Sigma~"E"^"s"^~mu~"        " %->%  "       +            (More Labile)")) +
xlab(bquote("More Buffered     -         " 
	%<-% "Relatative effect of environmental variation ("~ abs(sum(E[v]^sigma))~" )" %->%  
									"       +            More Labile"))+
ylab(bquote("Stochastic growth rate ("~lambda[s]~")"))+
#facet_grid(.~name,scales="free")+
#theme(legend.position="top")+
theme_minimal(base_size=21)+
theme(aspect.ratio=6/18)

#ggsave(file="Figures/Buffering continuum.svg")


#'---------------------------------------------------------------------------------------
##	----- 11.2. FIGURE 2. BUFFERING AND VITAL RATES ----- 
#'---------------------------------------------------------------------------------------
taxaLev<-c("Chordata","Arthropoda","Cnidaria",
	"Annual","Herbaceous perennial","Succulent",
		"Epiphyte","Shrub","Tree","Palm")

traitLev<-c("Cumulative","Survival","Reproduction","Growth","Shrinking","Clonality")

# Create the database
plotbuff_data<-final_data%>%
mutate(Taxa=ifelse(Kingdom=="Plantae", OrganismType,Phylum))%>%
select(c(Reproduction_SigElas:Survival_Base),MatRep,Kingdom,Taxa)%>%
pivot_longer(!c(Kingdom,MatRep,Taxa),names_to="Variable",values_to="Value")%>%
group_by(Kingdom,Taxa,Variable)%>%
summarise(Mean=mean(Value),
		SD=sd(Value),
		Populations=n(),
		MatRep=sum(MatRep),
		SE=SD/sqrt(Populations),
		CI=1.96*(SD/sqrt(Populations)))%>%
separate(Variable, c("Variable", "Form"))%>%
mutate(Taxa = factor(Taxa , levels = taxaLev))%>%
mutate(Variable = factor(Variable, levels = traitLev))

plotbuff_data

plotbuff<-plotbuff_data%>%
mutate(Variable=factor(Variable, levels = c("Cumulative","Survival", "Growth", "Shrinking", "Reproduction","Clonality")))%>%
ggplot(.,aes(x=Variable,y=Mean,group=Taxa))+
geom_pointrange(position=position_dodge(0.4),
 aes(ymin = Mean-SE, ymax = Mean+SE,
 	fill=Taxa,shape=Kingdom),color="black",linewidth = 1.3,size= 1.1,alpha=.9)+
ylab(NULL)+xlab(NULL)+
scale_shape_manual(values=c(23,21))+
scale_fill_viridis_d()+
#scale_color_viridis_d()+
#guides(fill = guide_legend(override.aes=list(shape=c(23,21))))+
guides(
alpha= FALSE,
fill = guide_legend(override.aes=list(shape=c(23,21),alpha=1)))+
theme_light(base_size=18)+
facet_grid(Form~Kingdom,scales="free")+
coord_flip()

plotbuff


p1<-p2<-plotbuff

p1$data<-plotbuff$data%>%
#filter(Kingdom=="Animalia")%>%
filter(Form=="Base")

p1<-p1+ylab(expression(
	paste("Vital rate elasticity ( ",E[v]," )")))

p2$data<-plotbuff$data%>%
filter(Form=="SigElas")

p2<-p2+ylab(bquote("Effect of environmental variation ("~ sum(E[v]^sigma)~" )"))


x11(height=862, width=1106);cowplot::plot_grid(p1,p2,ncol=1,labels="AUTO")

#ggsave(file="Figures/Vital rate contribution_new.svg")


#======================================================================================================
#			----- 11.3 FIGURE 3. Posterior distribution -------
#======================================================================================================
ggplot_posteriors<-Posterior_data%>%
  filter(Trait!="Buffmx")%>%
  ggplot(.,aes(x=Variables,y=Values,group=Model))+
  geom_hline(yintercept=0,linetype=2,color="grey50",linewidth=1.4)+
  stat_pointinterval(position=position_dodge(.5),.width = c(.66, .95),
                     aes(x = Variables,color=sig,shape=Model,fill=Model,alpha=sig))+
  scale_alpha_manual(values=c(.2,1))+
  scale_fill_manual(values=c("#da1438","#a783ce"))+
  scale_color_manual(values=c("grey30","black"))+
  scale_shape_manual(values=c(21,21))+
  xlab("Variables")+ylab("Posterior distribution")+
  ggh4x::facet_grid2(Trait~Taxa,scales="free_x",independent = "x")+
  theme_minimal(base_size=16)+coord_flip()+
  theme(
    #axis.text.x = element_text(angle = 45, hjust=1),
    legend.position="top")+
  guides(shape = guide_legend(override.aes = list(size = 5)))



ggplot_posteriors$data%>%glimpse()

ggplot_posteriors_vr<-ggplot_posteriors_cumu<-ggplot_posteriors

ggplot_posteriors_cumu$data<-filter(ggplot_posteriors$data,Taxa=="Plants" & Trait == "Cumulative")
ggplot_posteriors_vr$data<-filter(ggplot_posteriors$data,Taxa=="Plants" & Trait != "Cumulative")

# All vital rates + cumulative + All taxa
ggplot_posteriors

#  cumulative  & Plants only
ggplot_posteriors_cumu

#  cumulative vr & Plants only
ggplot_posteriors_vr









#'-------------------------------------------------------------------------
#	----- 11.3. BUFFERING TREE OF LIFE ----- 
#           NOT INCLUDED!!!
#'-------------------------------------------------------------------------
library(ggtree)
library(ggnewscale)
groupTaxa<-NULL
groupTaxa$Animals<-subtree_Animals$tip
groupTaxa$Plants<-subtree_Plants$tip

subtree_grouped<- groupOTU(subtree, groupTaxa)

phylo_circ <- phylo_circ <- ggtree(subtree_grouped,layout = "circular", 
	aes(color=group))+
scale_color_manual(values=c("black","#144375","#1f9a59"))+
geom_tiplab(size=4,offset = 2)


data_tree_df<-final_data%>%
group_by(phylo)%>%
summarise(
	Rep=mean(Reproduction_SigElas),
	Sur = mean(Survival_SigElas),
	Clo = mean(Clonality_SigElas),
	Shr = mean(Shrinking_SigElas),
	Gro = mean(Growth_SigElas),
	Sum = mean(Cumulative_SigElas))%>%
#mutate(across(R:C,scales::rescale,to=c(-1,1)))%>%
#mutate(across(R:C,log10))%>%
data.frame()%>%column_to_rownames(., var = "phylo") 

data_tree_df%>%glimpse()

gheatmap(phylo_circ, data_tree_df, offset=0.1, width=.2,
               colnames_angle=95, colnames_offset_y = .5) +
    scale_fill_viridis(option="magma")


# ALTERNATIVAMENTE

phylo1<-gheatmap(phylo_circ, dplyr::select(data_tree_df,Rep), offset=-.5, width=.15,
               colnames_angle=95, colnames_offset_y = -0.5) +
 scale_fill_gradientn(colours = colorRampPalette(c('#009593', '#E9E6C5', '#D35C79'))(3))

phylo1 ; phylo2<-phylo1+new_scale_fill()

phylo2<-gheatmap(phylo2, dplyr::select(data_tree_df,Sur ), offset=0.1, width=.15,
               colnames_angle=95, colnames_offset_y = .5) +
 scale_fill_gradientn(colours = colorRampPalette(c('#089392', '#E9E6C5', '#CF597E'))(5))

phylo2; phylo3<-phylo2+new_scale_fill()

phylo3<-gheatmap(phylo3, dplyr::select(data_tree_df,Sum), offset=.75, width=.15,
               colnames_angle=95, colnames_offset_y = .5) +
 scale_fill_gradientn(colours = colorRampPalette(c('#008585', '#E9E6C5', '#C7522B'))(5))

phylo3



#'===================================================================
#	---- 12. FINAL METADATA ---- 
#'===================================================================

#'-----------------------------------------------------------------------
# Include DOI to retrive full reference list of studies used
#'-----------------------------------------------------------------------
load("Data/COMADRE_v.4.23.3.1.RData")
load("Data/COMPADRE_v.6.23.5.0.RData")

compadre <- as_cdb(compadre)
compadre$StudyID<-cdb_id_studies(compadre)

comadre <- as_cdb(comadre)
comadre$StudyID<-cdb_id_studies(comadre)

compadre_db<-rbind(comadre@data,compadre@data)%>%
distinct(SpeciesAccepted,StudyID,StudyStart,StudyEnd,Lat,Lon,DOI_ISBN,Authors,Journal,YearPublication,SourceType)
#'-----------------------------------------------------------------------
library(rcrossref)

Final_metadata<-final_data%>%
mutate(lambda=Stochslambs_df$approx)%>%
left_join(by="ID",.,
Metadata%>%select(ID,StudyID,StudyStart,StudyEnd,Lat,Lon)%>%distinct(ID,.keep_all=T))%>%
select(Kingdom,SpeciesAccepted,OrganismType,StudyID,StudyStart,StudyEnd,Ecoregion,Lat,Lon,lambda,Cumulative_SigElas)%>%
left_join(.,compadre_db,by=c("SpeciesAccepted","StudyID","StudyStart","StudyEnd","Lat","Lon"))

Final_metadata%>%glimpse()

Final_metadata$StudyStart%>%min()
Final_metadata$StudyEnd%>%max()

mean(as.numeric(Final_metadata$StudyEnd)-as.numeric(Final_metadata$StudyStart))
sd(as.numeric(Final_metadata$StudyEnd)-as.numeric(Final_metadata$StudyStart))

DOIs<-unique(Final_metadata$DOI_ISBN)
REFS<-NULL

for (i in 1:length(DOIs)){
tryCatch({
REFS[[i]]<-cr_cn(DOIs[[i]], format = "text", style = "nature-neuroscience-brief-communications")},
	error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
print(i)
}

refere_list<-cbind(DOIs,REFS=unlist(REFs))%>%data.frame()

Final_metadata<-Final_metadata%>%
left_join(.,refere_list,by=c("DOI_ISBN"="DOIs"))%>%
mutate(REFS=ifelse(REFS!="NA",REFS,paste0(Authors," ",SourceType," ","(",YearPublication,").")))


Final_metadata
#openxlsx::write.xlsx(Final_metadata, file = "Supplementary material Table S2.xlsx",append = FALSE)


final_data%>%
left_join(by="ID",.,
Metadata%>%select(ID,StudyID,StudyStart,StudyEnd)%>%distinct(ID,.keep_all=T))%>%glimpse()

Metadata%>%glimpse()