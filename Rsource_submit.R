
### Perform MCMCglmm

All_family <- read.table("/path/to/All_family_data.csv",sep=",",header=T)

install.packages("MCMCglmm")
library(MCMCglmm)
install.packages("ape")
library(ape)
install.packages("phytools")
library(phytools)

# Brassicaceae
Bra_phylo <- read.tree("/path/to/Bra_newick_genus_sp.txt")
Bra_phylo$node.label <- seq(1, length(Bra_phylo$node.label))
Bra_phylo_forced <- force.ultrametric(Bra_phylo, method = "extend")
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
inv.phylo<-inverseA(Bra_phylo_forced,nodes="TIPS",scale=TRUE)
Bra_sub <- All_family[All_family$Family=="Brassicaceae",]
Bra_model <-MCMCglmm(log10(Seed_weight)~as.factor(Mating_system)+as.factor(Growth_form_herb),random=~Species_name,
                     family="gaussian",ginverse=list(Species_name=inv.phylo$Ainv),prior=prior,
                     data=Bra_sub,nitt=1000000,burnin=1000,thin=1)

# Asteraceae
Ast_phylo <- read.tree("/path/to/Ast_newick_genus_sp.txt")
Ast_phylo$node.label <- seq(1, length(Ast_phylo$node.label))
Ast_phylo_forced <- force.ultrametric(Ast_phylo, method = "extend")
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
inv.phylo<-inverseA(Ast_phylo_forced,nodes="TIPS",scale=TRUE)
Ast_sub <- All_family[All_family$Family=="Asteraceae",]
Ast_model <- MCMCglmm(log10(Seed_weight)~as.factor(Mating_system)+as.factor(Growth_form_herb),random=~Species_name,
                      family="gaussian",ginverse=list(Species_name=inv.phylo$Ainv),prior=prior,
                      data=Ast_sub,nitt=1000000,burnin=1000,thin=1)

# Solanaceae
Sol_phylo <- read.tree("/path/to/Sol_newick_genus_sp.txt")
Sol_phylo$node.label <- seq(1, length(Sol_phylo$node.label))
Sol_phylo_forced <- force.ultrametric(Sol_phylo, method = "extend")
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
inv.phylo<-inverseA(Sol_phylo_forced,nodes="TIPS",scale=TRUE)
Sol_sub <- All_family[All_family$Family=="Solanaceae",]
Sol_model <- MCMCglmm(log10(Seed_weight)~as.factor(Mating_system)+as.factor(Growth_form_herb),random=~Species_name,
                      family="gaussian",ginverse=list(Species_name=inv.phylo$Ainv),prior=prior,
                      data=Sol_sub,nitt=1000000,burnin=1000,thin=1)

### Perform linear regression analysis

anova(lm(log10(Seed_weight)~as.factor(Mating_system)+as.factor(Growth_form_herb),data=Ast_sub))
anova(lm(log10(Seed_weight)~as.factor(Mating_system)+as.factor(Growth_form_herb),data=Bra_sub))
anova(lm(log10(Seed_weight)~as.factor(Mating_system)+as.factor(Growth_form_herb),data=Sol_sub))

anova(lm(log10(Seed_weight)~as.factor(Genus),data=Ast_sub))
anova(lm(log10(Seed_weight)~as.factor(Genus),data=Bra_sub))
anova(lm(log10(Seed_weight)~as.factor(Genus),data=Sol_sub))

### Draw Fig.1 of the paper

install.packages("ggplot2")
install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

g1<- ggplot(All_family, aes(Mating_system, log10(Seed_weight), color = Mating_system))+  
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = 0.1, alpha = 0.7)+
  facet_wrap(~Family, nrow = 1)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")
g2<- ggplot(All_family, aes(Growth_form_herb, log10(Seed_weight), color = Growth_form_herb))+  
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = 0.1, alpha = 0.7)+
  facet_wrap(~Family, nrow = 1)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

grid.arrange(g1, g2, ncol=1)

