###########################################################
#	  		PHYLOGENETIC SIGNAL  
#		Blomberg's K and Pagels' Lambda
# 			a sript for  the project
#	DEMOGRAPHIC BUFFERING CONTINUUM in PLANTS AND ANIMALS
#			 by Gabriel Santos
# 		contact by ssantos.gabriel@gmail.com
#			    25 Aug 2023
###########################################################


library(phytools)

#================================================================================
#				PHYLOGENETIC SIGNAL
# 				!!!!!IMPORTANT!!!!!
# From now on, all analyses use absolute value of demographic buffering
#================================================================================


#================================================================================
#	CALCULATE BLOMBERG'S K FOR MULTIPLE VARIABLES
#================================================================================
#	 Acessory funtion
#--------------------------------------------------------------------------------

multiBlombergK<-function(X,tree=tree,nsim=nsim,test=TRUE){
MeanSignal<-temp<-out<-NULL
SimSignal<-matrix(NA,nsim,dim(X)[2])
	for(i in 1:dim(X)[2]){	
print(paste0("Calculing phylogenetic signal of: ",names(X)[i]))
	temp<- phylosig(tree,
			setNames(X[,i],rownames(X)),
				method="K",test=test,nsim=nsim)
MeanSignal[[1]]$Variables[i]<-names(X)[i]
MeanSignal[[1]]$Signal[i]<-temp$K
MeanSignal[[1]]$p_values[i]<-round(temp$P,4)
SimSignal[,i]<-temp[[3]]
			}
 colnames(SimSignal)<-names(X)
 rownames(SimSignal)<-paste0("Sim",rep(1:nsim, each = 1))
return(list(as.data.frame(MeanSignal),SimSignal))
}
#------------------------------------------------------------------------

PhyloSig_data_ALL%>%head()
Ks_All<-multiBlombergK(PhyloSig_data_ALL,tree=subtree,nsim=1000,test=TRUE)
Ks_Plants<-multiBlombergK(PhyloSig_data_PLANTS,tree=subtree_Plants,nsim=1000,test=TRUE)
Ks_Animals<-multiBlombergK(PhyloSig_data_ANIMALS,tree=subtree_Animals,nsim=1000,test=TRUE)


#================================================================================
#	Pagel's lambda
#================================================================================

multiPagelslamb<-function(X,tree=tree,test=TRUE){
MeanSignal<-temp<-out<-NULL
	for(i in 1:dim(X)[2]){	
print(paste0("Calculing phylogenetic signal of: ",names(X)[i]))
	temp<- phylosig(tree,
			setNames(X[,i],rownames(X)),
				method="lambda",test=test)
MeanSignal[[1]]$Variables[i]<-names(X)[i]
MeanSignal[[1]]$Signal[i]<-temp$lambda
MeanSignal[[1]]$p_values[i]<-round(temp$P,4)
			}
return(list(as.data.frame(MeanSignal)))
}


Lambs_All<-multiPagelslamb(PhyloSig_data_ALL,tree=subtree,test=TRUE)
Lambs_Plants<-multiPagelslamb(PhyloSig_data_PLANTS,tree=subtree_Plants,test=TRUE)
Lambs_Animals<-multiPagelslamb(PhyloSig_data_ANIMALS,tree=subtree_Animals,test=TRUE)



#================================================================================
#				RESULTS
#================================================================================

#BLOMBERG'S K
Ks_All[[1]]%>%filter(p_values<.1)
Ks_Plants[[1]]%>%filter(p_values<.1)
Ks_Animals[[1]]%>%filter(p_values<.1)

rbind(
data.frame(Ks_All[[2]],Taxa="ALL"),
data.frame(Ks_Plants[[2]],Taxa="Plants"),
data.frame(Ks_Animals[[2]],Taxa="Animals"))%>%
pivot_longer(.,!Taxa,names_to="Variables",values_to="Signal")%>%
ggplot(.,aes(x=Signal,group=Taxa))+
geom_histogram(aes(fill=Taxa),position="identity",alpha=0.7)+
geom_vline(xintercept=1,color="orange",linetype=2)+
facet_wrap(.~Variables,scale='free')+
theme_bw()


#Pagel's Lambda
Lambs_All[[1]]%>%filter(p_values<.1)
Lambs_Plants[[1]]%>%filter(p_values<.1)
Lambs_Animals[[1]]%>%filter(p_values<.1)

