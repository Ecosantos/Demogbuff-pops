#'##########################################################
#	  	 MODEL SELECTION
# 	a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 	contact by ssantos.gabriel@gmail.com
#			    05 March 2025
#' ---------------------------------------------------------
# Rationale: We must decide between better model composition
# as Environmental PCA included until 3 axes 
# revealing these three potential relevant source of information #
#'##########################################################

# Define necessary packages
need_pkgs <- c("tidyverse", "plotMCMC", "mcmcr", "MuMIn", "MCMCglmm","MuMIn")

exist_pckgs<-exist_pckgs <- installed.packages()[, "Package"]

if (any(!need_pkgs %in% exist_pckgs)) {   # Check for inexisting packages and install them
  install.packages(need_pkgs[!need_pkgs %in% exist_pckgs])
}


# load non-existing packages
lapply(need_pkgs, require, character.only = TRUE)


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

 

data_model%>%glimpse()

hist(data_model$Cumulative_SigElas) # must be absolute values!

#'===========================================================
# ---- Define priors ----
#'===========================================================


prior_phylo<-list(G=list(G1=list(V=1,nu=0.02)),
                  R=list(V=1,nu=0.02))

# Prior for simple model will be removed. 
# It will run with prior_phylo for standardization!
#prior_simple<-list(G=list(R=list(V=1,nu=0.02))) 

nitt=100000; #nitt=1000
burnin=1000; #burnin=100 
thin=10	

glmmScale<-"FALSE"

#'===========================================================
# ---- Create models manually ----
#'===========================================================
# A function was created to reduce code
fit_mcmcglmm <- function(formula) {
  MCMCglmm(
    formula,  # Garantir que a fórmula seja interpretada corretamente
    random = ~phylo,
    family = "gaussian",
    ginverse = list(phylo = inverseA(subtree_Plants, nodes = "TIPS", scale = TRUE)$Ainv),
    prior = prior_phylo,
    data = subset(data_model, Kingdom == "Plantae"),
    nitt = nitt,
    burnin = burnin,
    thin = thin,
    singular.ok = TRUE,
    scale = glmmScale
  )
}


# Models with interactions
Clim123 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.2 * ClimPC.3)
Clim12 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.2)
Clim13 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.3)
Clim23 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.2 * ClimPC.3)

#Aditive models
Clim123_plus<-fit_mcmcglmm(Cumulative_SigElas~LHPC.1 * LHPC.2 + ClimPC.1 + ClimPC.2 + ClimPC.3)
Clim12_plus<-fit_mcmcglmm(Cumulative_SigElas~LHPC.1 * LHPC.2 + ClimPC.1 + ClimPC.2)

# Aditive for life history and climatic data
Clim12_plus_only<-fit_mcmcglmm(Cumulative_SigElas~LHPC.1 + LHPC.2 + ClimPC.1 + ClimPC.2) 

# Models with single climatic variable
Clim1 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.1)
Clim2 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.2)
Clim3 <- fit_mcmcglmm(Cumulative_SigElas ~ LHPC.1 * LHPC.2 + ClimPC.3)


forms<-lapply(
  lapply(
    list(Clim1,Clim2,Clim3,
         Clim123,
         Clim12,Clim13,Clim23,
         Clim123_plus,Clim12_plus,
         Clim12_plus_only),formula),
          function(f) paste(deparse(formula(f)), collapse = ""))%>%do.call(rbind,.)

mod_out<-cbind(formulas=forms,
  BIC(Clim1,Clim2,Clim3,Clim123,Clim12,Clim13,Clim23,Clim123_plus,Clim12_plus,Clim12_plus_only),
  AIC=AIC(Clim1,Clim2,Clim3,Clim123,Clim12,Clim13,Clim23,Clim123_plus,Clim12_plus,Clim12_plus_only)[,2],
  DIC=DIC(Clim1,Clim2,Clim3,Clim123,Clim12,Clim13,Clim2,Clim123_plus,Clim12_plus,Clim12_plus_only)[,2])%>%
  arrange(DIC)

mod_out


Clim12%>%summary()

allChains <- mcmcr::as.mcmc(cbind(Clim12$Sol,Clim12$VCV))
plotMCMC::plotTrace(allChains)


# Compare some models if necessary
#summary(Clim123)
#summary(Clim123_plus)
#summary(Clim12_plus)
#summary(Clim1)
#summary(Clim12)
#summary(Clim13)







