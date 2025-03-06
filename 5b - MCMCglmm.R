set.seed(1)

library(MCMCglmm)

GLMMdata_link<-"https://github.com/Ecosantos/Demogbuff-pops/raw/refs/heads/incorporating-MCMCGlmm/Data/GLMMdata.Rdata"
download.file(GLMMdata_link, "Data/GLMMdata.Rdata", mode = "wb")

load("https://github.com/Ecosantos/Demogbuff-pops/raw/refs/heads/incorporating-MCMCGlmm/Data/GLMMdata.Rdata")

data_model<-final_data%>%select(-c(Reproduction_Base:Cumulative_Base))


MCMCglmm_phylo_plants<-MCMCglmm_phylo_animals<-NULL



fixEffect<-fixEffect<-"~LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.2"
InterestingVars<-c("Survival","Growth","Shrinking","Reproduction","Clonality","Buffmx","Cumulative")

traits<-traits_glmm<- unique (grep(paste(InterestingVars,collapse="|"), 
                                   colnames(data_model), value=TRUE))

#------------------------------------------------------------------------------------------------------
#	Define priors and GLMM iteractions
#------------------------------------------------------------------------------------------------------

prior_phylo<-list(G=list(G1=list(V=1,nu=0.02)),
   R=list(V=1,nu=0.02))

#prior_simple<-list(G=list(R=list(V=1,nu=0.02)))

nitt=51000; #nitt=1000
burnin=1000; #burnin=100 
thin=3	

glmmScale<-"FALSE"


#------------------------------------------------------------------------------------------------------
#	PHYLOGENETIC MCMC GLMM
#------------------------------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------------------------------
#	SIMPLE MCMC GLMM
#------------------------------------------------------------------------------------------------------

MCMCglmm_simple_animals<-MCMCglmm_simple_plants<-NULL

# ATUALMENTE USING RANDOM AS DEFAULT! CHANGE IT IN THE FUTURE!

#Animals	phylo
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



glmmOUT<-list(Simple_models=
       list(Plants = MCMCglmm_simple_plants,Animals = MCMCglmm_simple_animals),
     Phylogenetic_models= 
       list(Plants = MCMCglmm_phylo_plants,Animals = MCMCglmm_phylo_animals)
     )

saveRDS(glmmOUT,"Data/MCMCglmm_output.rds")

