#'##########################################################
#	  	      MCMC GLMM
# 	a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 	contact by ssantos.gabriel@gmail.com
#			    05 March 2025
#' ---------------------------------------------------------
# MCMCglmms run separately so we can run it
# with intermediary data produced
#'##########################################################

set.seed(1)
rm(list=ls())


# Define necessary packages
need_pkgs <- c("tidyverse", "plotMCMC", "mcmcr", "MuMIn", "MCMCglmm","MuMIn")

exist_pckgs<-exist_pckgs <- installed.packages()[, "Package"]

if (any(!need_pkgs %in% exist_pckgs)) {   # Check for inexisting packages and install them
  install.packages(need_pkgs[!need_pkgs %in% exist_pckgs])
}



rm(list=ls())

#'===========================================================
# ----- LOAD data ready for GLMM analyses ----
# GLMM data contains:
#  - Data model
#  - subtree_Animals - Phylogenetic tree for Animals ready for analysis
#  - subtree_Animals - Phylogenetic tree for Plants ready for analysis
#'===========================================================
# Download data In Google colab
#  GLMMdata_link<-"https://github.com/Ecosantos/Demogbuff-pops/raw/refs/heads/incorporating-MCMCGlmm/Data/GLMMdata.Rdata"
# download.file(GLMMdata_link, "Data/GLMMdata.Rdata", mode = "wb")

load("Data/GLMMdata.Rdata")

#'===========================================================
# ---- Define GLMM SETTINGS ----
#'===========================================================

fixEffect<-fixEffect<-"~LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.2"
InterestingVars<-c("Survival","Growth","Shrinking","Reproduction","Clonality","Buffmx","Cumulative")

traits<-traits_glmm<- unique (grep(paste(InterestingVars,collapse="|"), 
                                   colnames(data_model), value=TRUE))

#'-----------------------------------------------------------
# ---- Define priors ----
#'-----------------------------------------------------------

prior_phylo<-list(G=list(G1=list(V=1,nu=0.02)),
                  R=list(V=1,nu=0.02))

nitt=100000; #nitt=1000
burnin=1000; #burnin=100 
thin=10	

glmmScale<-"FALSE"


#'===========================================================
# ---- RUN MCMCglmm MODELS  ----
#'===========================================================

#'-----------------------------------------------------------
##	----- PHYLOGENETIC MCMC GLMM -----
#'-----------------------------------------------------------
MCMCglmm_phylo_animals<-MCMCglmm_phylo_plants<-NULL

#Animals	phylo
for(i in 1:length(traits)){
print(paste("Running model:",traits[i],"~",fixEffect))
MCMCglmm_phylo_animals[[i]]<-MCMCglmm(formula(paste0(traits[i], fixEffect)),
    random=~phylo,family="gaussian",
		ginverse=list(phylo=inverseA(subtree_Animals,nodes="TIPS",scale=TRUE)$Ainv),
				prior=prior_phylo,data=subset(data_model,Kingdom=="Animalia"),
   						nitt=nitt,burnin=burnin,thin=thin,singular.ok=TRUE, scale = glmmScale)

names(MCMCglmm_phylo_animals)[i]<-traits[[i]]
}


#Plants
for(i in 1:length(traits)){
print(paste("Running model:",traits[i],"~",fixEffect))
MCMCglmm_phylo_plants[[i]]<-MCMCglmm(formula(paste0(traits[i], fixEffect)),
    random=~phylo,family="gaussian",
		ginverse=list(phylo=inverseA(subtree_Plants,nodes="TIPS",scale=TRUE)$Ainv),
				prior=prior_phylo,data=subset(data_model,Kingdom=="Plantae"),
   						nitt=nitt,burnin=burnin,thin=thin,singular.ok=TRUE, scale = glmmScale)

names(MCMCglmm_phylo_plants)[i]<-traits[[i]]
}

#'===========================================================
##	----- SIMPLE MCMC GLMM -----
#'===========================================================
MCMCglmm_simple_animals<-MCMCglmm_simple_plants<-NULL

#Animals
for(i in 1:length(traits)){
print(paste("Running model:",traits[i],"~",fixEffect))
MCMCglmm_simple_animals[[i]]<-MCMCglmm(formula(paste0(traits[i], fixEffect)),
    family="gaussian",data=subset(data_model,Kingdom=="Animalia"),
   			nitt=nitt,burnin=burnin,thin=thin,singular.ok=TRUE, scale = glmmScale)

names(MCMCglmm_simple_animals)[i]<-traits[[i]]
}

#Plants
for(i in 1:length(traits)){
print(paste("Running model:",traits[i],"~",fixEffect))
MCMCglmm_simple_plants[[i]]<-MCMCglmm(formula(paste0(traits[i], fixEffect)),
   family="gaussian",data=subset(data_model,Kingdom=="Plantae"),
   			nitt=nitt,burnin=burnin,thin=thin,singular.ok=TRUE, scale = glmmScale)

names(MCMCglmm_simple_plants)[i]<-traits[[i]]
}


#'===========================================================
#	----- EXPORT GLMM OUTPUTs -----
#'===========================================================

glmmOUT<-list(Simple_models=
       list(Plants = MCMCglmm_simple_plants,Animals = MCMCglmm_simple_animals),
     Phylogenetic_models= 
       list(Plants = MCMCglmm_phylo_plants,Animals = MCMCglmm_phylo_animals)
     )

#saveRDS(glmmOUT,"Data/MCMCglmm_output.rds")

