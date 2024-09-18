set.seed(1)
data_model<-final_data%>%select(-c(Reproduction_Base:Cumulative_Base))

data_model%>%glimpse()
#Transform CUMULATIVE IN ABSOLUTE VALUE

data_model$Cumulative_SigElas<-abs(data_model$Cumulative_SigElas)


data_model%>%glimpse()

fixEffect<-fixEffect<-"~LHPC.1 * LHPC.2 + ClimPC.1 * ClimPC.2"
InterestingVars<-c("Survival","Growth","Shrinking","Reproduction","Clonality","Buffmx","Cumulative")

traits_glmm<- unique (grep(paste(InterestingVars,collapse="|"), 
                        colnames(data_model), value=TRUE))


traits<-traits_glmm

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

MCMCglmm_phylo_plants<-MCMCglmm_phylo_animals<-NULL

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

#======================================================================================================
#	MERGING RESULTS
#======================================================================================================

#------------------------------------------------------------------------------------------------------
#			As summary
#------------------------------------------------------------------------------------------------------
Phylo_models_animals<-lapply(MCMCglmm_phylo_animals,function(inner_list) summary(inner_list))
Phylo_models_plants<-lapply(MCMCglmm_phylo_plants,function(inner_list) summary(inner_list))

Simple_models_animals<-lapply(MCMCglmm_simple_animals,function(inner_list) summary(inner_list))
Simple_models_plants<-lapply(MCMCglmm_simple_plants,function(inner_list) summary(inner_list))


Phylo_models_animals_coefs<-lapply(Phylo_models_animals,function(inner_list) inner_list[[5]])
Phylo_models_plants_coefs<-lapply(Phylo_models_plants,function(inner_list) inner_list[[5]])

Simple_models_animals_coefs<-lapply(Simple_models_animals,function(inner_list) inner_list[[5]])
Simple_models_plants_coefs<-lapply(Simple_models_plants,function(inner_list) inner_list[[5]])


#---------------------------------------------------------------------------------------------------
#			PHYLO 1
#---------------------------------------------------------------------------------------------------
Phylo_df<-lapply(Phylo_models_plants,function(inner_list) inner_list[[6]])%>%
Map(cbind, Trait = names(.),.)%>%do.call(rbind,.)%>%
as_tibble()%>%
mutate(across(c("post.mean":"eff.samp"),as.numeric))


#---------------------------------------------------------------------------------------------------
#Transform to data.frame
#---------------------------------------------------------------------------------------------------
Phylo_models_animals_coefs<-lapply(Phylo_models_animals_coefs,as.data.frame)
Phylo_models_plants_coefs<-lapply(Phylo_models_plants_coefs,as.data.frame)

Simple_models_animals_coefs<-lapply(Simple_models_animals_coefs,as.data.frame)
Simple_models_plants_coefs<-lapply(Simple_models_plants_coefs,as.data.frame)


Phylo_models_df_plants<-lapply(Phylo_models_plants_coefs,rownames_to_column, var = "Statistics")%>%
Map(cbind, Trait = names(.),Taxa="Plants", Model="Phylo", .)%>%do.call(rbind,.)

Phylo_models_df_animals<-lapply(Phylo_models_animals_coefs,rownames_to_column, var = "Statistics")%>%
Map(cbind, Trait = names(.),Taxa="Animals",Model="Phylo", .)%>%do.call(rbind,.)


Simple_models_df_plants<-lapply(Simple_models_plants_coefs,rownames_to_column, var = "Statistics")%>%
Map(cbind, Trait = names(.),Taxa="Plants", Model="Simple", .)%>%do.call(rbind,.)

Simple_models_df_animals<-lapply(Simple_models_animals_coefs,rownames_to_column, var = "Statistics")%>%
Map(cbind, Trait = names(.),Taxa="Animals",Model="Simple", .)%>%do.call(rbind,.)



GLMMs_df<-rbind(Phylo_models_df_plants,Phylo_models_df_animals,
		Simple_models_df_plants,Simple_models_df_animals)


colnames(GLMMs_df)<-c("Trait","Taxa","Model","Statistics","post.mean","low95","high95","eff.samp","pMCMC")

GLMMs_df%>%
mutate(sig=ifelse(pMCMC<=0.05,"Sig","Non-Sig"))%>%
filter(Statistics!="(Intercept)")%>%
ggplot(.,aes(x=Statistics,y=post.mean,group=Model))+
geom_hline(yintercept=0,linetype=2,color="grey50",size=1.4)+
geom_pointrange(position=position_dodge(0.5),
	aes(ymin = low95, ymax = high95,fill=Taxa,alpha=sig,shape=Model))+
scale_alpha_manual(values=c(.4,1))+
scale_shape_manual(values=c(21,22))+
facet_grid(Taxa~Trait,scales="free_y")+
theme_bw(base_size=16)+
theme(axis.text.x = element_text(angle = 45, hjust=1))


GLMMs_df%>%
mutate(sig=ifelse(pMCMC<=0.05,"Sig","Non-Sig"))%>%
filter(pMCMC<=0.1)%>%
filter(Statistics!="(Intercept)")%>%
group_by(Trait)%>%group_split()



GLMMs_df%>%
#filter(Taxa=="Plants")%>%
separate(Trait, c("Trait", "Form"))%>%
filter(Form!="Base")%>%
mutate(sig=ifelse(pMCMC<=0.05,"Sig","Non-Sig"))%>%
filter(Statistics!="(Intercept)")%>%
ggplot(.,aes(x=Statistics,y=post.mean,group=Model))+
geom_hline(yintercept=0,linetype=2,color="grey50",size=1.4)+
geom_pointrange(position=position_dodge(0.75),
	aes(ymin = low95, ymax = high95,fill=Model,alpha=sig,shape=Model))+
scale_alpha_manual(values=c(.4,1))+
scale_shape_manual(values=c(21,22))+
scale_fill_brewer(palette = "Dark2")+
ggh4x::facet_grid2(Taxa~Trait,scales="free_x",independent = "x")+
theme_minimal(base_size=16)+
theme(axis.text.x = element_text(angle = 45, hjust=1))+coord_flip()




#------------------------------------------------------------------------------------------------------
#			As Raw posterior distribution
#------------------------------------------------------------------------------------------------------
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

Posterior_data%>%glimpse()

ggplot_posteriors<-Posterior_data%>%
filter(Trait!="Buffmx")%>%
ggplot(.,aes(x=Variables,y=Values,group=Model))+
geom_hline(yintercept=0,linetype=2,color="grey50",size=1.4)+
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

ggplot_posteriors

#==========================================================================
#	PHYLOGENETIC SIGNAL
#==========================================================================

my.fake.lamb<-function(model){
 out =  model$VCV[,"phylo"]/
          (model$VCV[,"phylo"]+
#            model$VCV[,"species"]+
             model$VCV[,"units"])
mean.Lambda=mean(out)
SE.lambda=(sd(out)/sqrt(length(out)))
return(out)
}



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

Phylo_signal_df%>%
group_by(Taxa,Trait)%>%
summarise(Mean=mean(Values),
		SD=sd(Values),
		N=n(),
		gaussianCI=1.96*(SD/sqrt(N)),
lower95=quantile(Values,.025),
higher95=quantile(Values,.975))



library(ggridges)

Phylo_signal_df%>%glimpse()

Phylo_signal_df%>%
group_by(Taxa,Trait)%>%
filter(Trait!="Buffmx")%>%
filter(!(Values>quantile(Values,.975)))%>%
filter(!(Values<quantile(Values,.025)))%>%
ggplot(.,aes(x=Values, y =  fct_rev(Trait),fill=Taxa))+
#geom_density_ridges(aes(fill = Taxa), rel_min_height = 0.01)+
geom_density_ridges(aes(fill=Taxa))+
xlab("Pagel's lambda )")+
scale_fill_manual(values=c("#144375","#1f9a59"))+
theme_minimal(base_size=16)+
ylab(NULL)+
theme(axis.text.y=element_blank(),legend.position="top",
panel.spacing = unit(2, "lines"))+
facet_grid(Trait~.,scale="free_y")



cowplot::plot_grid(
rel_widths=c(1,.3),
ggplot_posteriors,
Phylo_signal_df%>%
group_by(Taxa,Trait)%>%
filter(Trait!="Buffmx")%>%
filter(!(Values>quantile(Values,.975)))%>%
filter(!(Values<quantile(Values,.025)))%>%
ggplot(.,aes(x=Values, y =  fct_rev(Trait),fill=Taxa))+
#geom_density_ridges(aes(fill = Taxa), rel_min_height = 0.01)+
geom_density_ridges(aes(fill=Taxa))+
xlab("Pagel's lambda )")+
scale_fill_manual(values=c("#144375","#1f9a59"))+
theme_minimal(base_size=16)+
ylab(NULL)+
theme(axis.text.y=element_blank(),legend.position="top",
panel.spacing = unit(2, "lines"))+
facet_grid(Trait~.,scale="free_y"))

#ggsave(file="Figures/MCMCglmm result.svg")

#==========================================================================
#	TRACEPLOT
#==========================================================================

#require(plotMCMC)

Check_traces<-function(X){
  windows(record=TRUE) # opens a window and starts recording
    op <- par(ask=TRUE)
allChains <- NULL
for(i in 1:length(X)){
	allChains <-as.mcmc(cbind(X[[i]]$Sol,X[[i]]$VCV))
plotMCMC::plotTrace(allChains,
	main=X[[i]]$Fixed$formula,cex.main = .6)
print ("Click on plot to continue")
#readline(prompt="Press [enter] to continue")
	}
    windows.options(record=FALSE) #stops recording.
}

#Check_traces(MCMCglmm_phylo_plants)
#Check_traces(MCMCglmm_phylo_animals)
#Check_traces(MCMCglmm_simple_plants)
#Check_traces(MCMCglmm_simple_animals)


#==========================================================================
#	PREDICTED x OBSERVED
#==========================================================================
#---------------------------------------------------------------------------
#		PLOT PLANTS
#---------------------------------------------------------------------------
par(mfrow=c(2,7))
TEMP_df<-TEMP_pred<-TEMP_obs<-NULL
for ( i in 1:length(traits)){
TEMP_pred<-MCMCglmm_phylo_plants[[i]]%>%predict()%>%as.vector()
TEMP_obs<-subset(data_model,Kingdom=="Plantae")%>%select(traits[[i]])%>%as.matrix()
TEMP_df<-filter(data.frame(TEMP_pred,TEMP_obs=TEMP_obs[,1]),TEMP_obs<=quantile(TEMP_obs,0.975))
TEMP_df<-filter(TEMP_df,TEMP_obs >= quantile(TEMP_obs,0.25))
plot(TEMP_df$TEMP_pred~TEMP_df$TEMP_obs,main=traits[[i]])
abline(a=0,b=1)
}

for ( i in 1:length(traits)){
TEMP_pred<-MCMCglmm_simple_plants[[i]]%>%predict()%>%as.vector()
TEMP_obs<-subset(data_model,Kingdom=="Plantae")%>%select(traits[[i]])%>%as.matrix()
TEMP_df<-filter(data.frame(TEMP_pred,TEMP_obs=TEMP_obs[,1]),TEMP_obs<=quantile(TEMP_obs,0.975))
TEMP_df<-filter(TEMP_df,TEMP_obs >= quantile(TEMP_obs,0.25))
plot(TEMP_df$TEMP_pred~TEMP_df$TEMP_obs,main=traits[[i]])
abline(a=0,b=1)
}



#---------------------------------------------------------------------------
#		PLOT PLANTS
#---------------------------------------------------------------------------
par(mfrow=c(2,7))
TEMP_df<-TEMP_pred<-TEMP_obs<-NULL
for ( i in 1:length(traits)){
TEMP_pred<-MCMCglmm_phylo_animals[[i]]%>%predict()%>%as.vector()
TEMP_obs<-subset(data_model,Kingdom=="Animalia")%>%select(traits[[i]])%>%as.matrix()
TEMP_df<-filter(data.frame(TEMP_pred,TEMP_obs=TEMP_obs[,1]),TEMP_obs<=quantile(TEMP_obs,0.975))
TEMP_df<-filter(TEMP_df,TEMP_obs >= quantile(TEMP_obs,0.25))
plot(TEMP_df$TEMP_pred~TEMP_df$TEMP_obs,main=traits[[i]])
abline(a=0,b=1)
}

for ( i in 1:length(traits)){
TEMP_pred<-MCMCglmm_simple_animals[[i]]%>%predict()%>%as.vector()
TEMP_obs<-subset(data_model,Kingdom=="Animalia")%>%select(traits[[i]])%>%as.matrix()
TEMP_df<-filter(data.frame(TEMP_pred,TEMP_obs=TEMP_obs[,1]),TEMP_obs<=quantile(TEMP_obs,0.975))
TEMP_df<-filter(TEMP_df,TEMP_obs >= quantile(TEMP_obs,0.25))
plot(TEMP_df$TEMP_pred~TEMP_df$TEMP_obs,main=traits[[i]])
abline(a=0,b=1)
}

