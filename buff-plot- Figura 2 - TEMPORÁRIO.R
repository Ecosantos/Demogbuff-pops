# PLOT BUFF - Figura 2 - 1a tentativa. Acho que não será necessário. Todas as mudanças serão feitas no Inkscape!
# Se não der certo usar aqui, por exemplo para alargar as figuras. 
# Esse script não está completo mas alguns bugs que apareciam já foram consertados

plotbuff<-plotbuff_data%>%
  mutate(Variable=factor(Variable, levels = c("Cumulative","Survival", "Growth", "Shrinking", "Reproduction","Clonality")))%>%
  ggplot(.,aes(x=Variable,y=Mean,group=Taxa))+
  geom_pointrange(position=position_dodge(0.4),
                  aes(ymin = Mean-SE, ymax = Mean+SE,
                      fill=Taxa,shape=Kingdom),color="black",linewidth = 1.3,size= 1.1,alpha=.9)+
  ylab(NULL)+xlab(NULL)+
  scale_fill_viridis_d()+
  #scale_color_viridis_d()+
  #guides(fill = guide_legend(override.aes=list(shape=c(23,21))))+
  guides(alpha= FALSE)+
  theme_light(base_size=18)+
  facet_grid(Form~Kingdom,scales="free")+
  coord_flip()

p1_plants<-p2_plants<-p1<-p2<-plotbuff

p1_plants<-p1<-p1+ylab(expression(paste("Vital rate elasticity ( ",E[v]," )")))

p1$data<-plotbuff$data%>%filter(Form=="Base")
p1<-p1+ scale_shape_manual(values=c(23,21)) +
        guides(alpha= FALSE,
           fill = guide_legend(override.aes=list(shape=c(23,21),alpha=1)))

p1