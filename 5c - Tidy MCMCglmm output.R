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

rm(list)

library(tidyverse)
library(tidybayes)	#v3.0.6
library(ggridges)
library(MuMIn);options(na.action = "na.fail")




#======================================================================================================
#	---- LOADING GLMM OUTPUTS -----
#======================================================================================================

MCMCglmm_output<-readRDS("Data/MCMCglmm_output.rds")

lapply(MCMCglmm_output,names)

# Phylogenetic models
MCMCglmm_phylo_plants<-MCMCglmm_output$Phylogenetic_models[[1]]
MCMCglmm_phylo_animals<-MCMCglmm_output$Phylogenetic_models[[2]]

# Non-Phylogenetic (simple) models
MCMCglmm_simple_plants<-MCMCglmm_output$Simple_models[[1]]
MCMCglmm_simple_animals<-MCMCglmm_output$Simple_models[[2]]

rm(MCMCglmm_output)


#======================================================================================================
#	---- DATA HARMONIZATION 	---- 
#======================================================================================================

#------------------------------------------------------------------------------------------------------
##			---- 	As summary 	---- 
#------------------------------------------------------------------------------------------------------
Phylo_models_animals<-lapply(MCMCglmm_phylo_animals,function(inner_list) summary(inner_list))
Phylo_models_plants<-lapply(MCMCglmm_phylo_plants,function(inner_list) summary(inner_list))

Simple_models_animals<-lapply(MCMCglmm_simple_animals,function(inner_list) summary(inner_list))
Simple_models_plants<-lapply(MCMCglmm_simple_plants,function(inner_list) summary(inner_list))


Phylo_models_animals_coefs<-lapply(Phylo_models_animals,function(inner_list) inner_list[[5]])
Phylo_models_plants_coefs<-lapply(Phylo_models_plants,function(inner_list) inner_list[[5]])

Simple_models_animals_coefs<-lapply(Simple_models_animals,function(inner_list) inner_list[[5]])
Simple_models_plants_coefs<-lapply(Simple_models_plants,function(inner_list) inner_list[[5]])


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


#---------------------------------------------------------------------------------------------------
# -----  SUMMARY SINTHESIS  -----  
#---------------------------------------------------------------------------------------------------

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


#======================================================================================================
#			----- Posterior distributions -------
#======================================================================================================
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


#======================================================================================================
#			----- GGPLOT - Posterior distribution -------
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


#---------------------------------------------------------------------------------------------------
#			PHYLO 1
# Don't remember why it is necessary
# Exclude in next versions
#---------------------------------------------------------------------------------------------------
#Phylo_df<-lapply(Phylo_models_plants,function(inner_list) inner_list[[6]])%>%
#  Map(cbind, Trait = names(.),.)%>%do.call(rbind,.)%>%
#  as_tibble()%>%
#  mutate(across(c("post.mean":"eff.samp"),as.numeric))

#==========================================================================
#	---- Estimating PHYLOGENETIC SIGNAL ---- 
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

Phylo_signal_df%>%
  group_by(Taxa,Trait)%>%
  summarise(Mean=mean(Values),
            SD=sd(Values),
            N=n(),
            gaussianCI=1.96*(SD/sqrt(N)),
            lower95=quantile(Values,.025),
            higher95=quantile(Values,.975))




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
#	----- COMPARING PHYLOGENETIC SIGNAL ---------
#==========================================================================

Phylo_signal_df%>%glimpse()

filter(Phylo_signal_df,Trait=="Reproduction" & Taxa=="Plants")$Values%>%range()


Phylo_signal_df%>%
  ggplot(.,aes(x=Trait))+
  geom_bar()


Phylo_summary <- Phylo_signal_df %>%
  filter(Taxa=="Plants" & !(Trait %in% c("Cumulative","Buffmx")))%>%
  group_by(Taxa,Trait) %>%
  summarise(
    MEDIAN = median(Values),
    SD = sd(Values),
    SE = SD / sqrt(n()) )%>%
  mutate(Trait=factor(Trait,levels = c("Survival", "Growth", "Shrinking", "Reproduction", "Clonality")))


ggplot(Phylo_summary, aes(x = Trait, y = MEDIAN, fill = Taxa)) +
  geom_bar(stat = "identity", position = position_dodge()) +  # Barras com transparência leve
  geom_pointrange(aes(ymin = MEDIAN - SD, ymax=pmin(MEDIAN + SD, 1)), 
                  position = position_dodge(width = 0.9), color = "black", size = 0.8) +
  scale_fill_manual(values=c("#1f9a59"))+
  labs(x = NULL,
       y = "Phylogenetic signal \n (Pagel's lambda)") +
  theme_minimal(base_size=18)+
  theme(legend.position="none",
        panel.spacing = unit(2, "lines"))



data_model


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

