######################## Load directory and Packages ######################## 
library(igraph)
library(plotly)
library(data.table)
library(betapart)
library(vegan)
library(reshape2)
library(ggplot2)
library(betalink)
library(imager)
library(knitr)
library(plyr)
library(bipartite)
library(igraph)
library(tidyverse)
library(janitor)
library(stringr)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(lmPerm)

######################## Data ######################## 
#setwd("C:/Users/Windows 7/OneDrive - ufpr.br/Doutorado/GUARAGUA?U/TESE")
int=read.csv2("Fish.parasite.csv")
host.traits=read.csv2("traits.csv")
parasite.traits=read.csv2("traits.p.csv")
host.sp.traits=read.csv2("traits.inter.csv")
env.var = read.csv("Env.form.csv",  sep = ";")
kn= read.csv2("kn.host.csv")
#parasite.traits = parasite.traits[,-c(12:14)]


########################  Group species by site ######################## 
int.pform= full_join(int, parasite.traits, by="Species")
int.ecto= int.pform %>% filter(Parasitism.form == "Ectoparasite")
int.endo= int.pform %>% filter(Parasitism.form == "Endoparasite")

nets.ect = dlply(int.ecto, .(LOC), acast, Cod.ind~Species, sum, value.var="ABD.x")
nets.endo = dlply(int.endo, .(LOC), acast, Cod.ind~Species, sum, value.var="ABD.x")
nets = dlply(int.pform, .(LOC), acast, Cod.ind~Species, sum, value.var="ABD.x")

nets.graph.ct = lapply(nets.ect, function(x) graph_from_incidence_matrix(x, weighted=TRUE))
nets.graph.nd = lapply(nets.endo, function(x) graph_from_incidence_matrix(x, weighted=TRUE))

########################  Modularity ######################## 

output.nd <- data.frame(raw=1, null=2)
CI.nd= data.frame(NULL)
for(len in 1:length(nets.endo)){
  mod = metaComputeModules(nets.endo[[len]], N=5, method="Beckett",steps = 999)
  output.nd[len,1] = mod@likelihood
  nulls = nullmodel(nets.endo[[len]], N=100, method="vaznull")
  mod.null = sapply(nulls, metaComputeModules, step=999)
  like.nulls = sapply(mod.null[len], function(x) x@likelihood)
  ci= data.frame(quantile(unlist(like.nulls), c(0.025, 0.975)))
  output.nd[len,2] = (mean(like.nulls))
  CI.nd= rbind(CI.nd, ci)
}


#modularity.endo = output %>% transform(Qs = raw-null) %>% transform(Z.value= Qs/null)

output.ct <- data.frame(raw=1, null=2)
CI.ct= data.frame(NULL)
for(len in 1:length(nets.ect))
{
  mod = metaComputeModules(nets.ect[[len]], N=5, method="Beckett",steps = 999)
  output.ct[len,1] = mod@likelihood
  nulls = nullmodel(nets.ect[[len]], N=100, method="vaznull")
  mod.null = sapply(nulls, metaComputeModules, step=999)
  like.nulls = sapply(mod.null, function(x) x@likelihood)
  ci= data.frame(quantile(unlist(like.nulls), c(0.025, 0.975)))
  output.ct[len,2] = (mean(like.nulls))
  CI.ct= rbind(CI.ct, ci)
}

modularity.ecto = output.ct %>% transform(Qs = raw-null) %>% transform(Z.value= Qs/null, Sector= c(1, 2, 3, 4), Form = "Ecto", Sample = c("1.ct", "2.ct", "3.ct","4.ct"))

modularity.endo = output.nd %>% transform(Qs = raw-null) %>% transform(Z.value= Qs/null, Sector= c(1, 2, 3, 4), Form = "Endo", Sample = c("1.nd", "2.nd", "3.nd","4.nd"))

mod.tot<- full_join(modularity.ecto, modularity.endo)
mod.env<- full_join(env.var, mod.tot)
write.csv(mod.env, "mod.Form.env.csv")

########################  Relation modularity ######################## 
mod.env.form = mod.env %>% group_by(Sector, Form, Sample) %>% summarise_all(mean) %>% ungroup()

#with(mod.env.nd, leveneTest(Rich.prst, Form))
mod.env.ct = mod.env.form%>% filter(Form == "Ecto")
mod.env.nd = mod.env.form %>% filter(Form == "Endo")
#vvariance not different from each other
str(mod.env.ct)
#Shapiro test para normalidade de dados

with(mod.env.ct, shapiro.test(Human_Structures))
with(mod.env.nd, shapiro.test(log(Conductivity)))
with(mod.env.ct, shapiro.test(Rich.prst))
with(mod.env.ct, shapiro.test(Int.prst))

#Rich and abund vs env variables
#Huma, structures
lm.h.ct = lm(Human_Structures~Int.prst, data=mod.env.ct)
summary(lm.h.ct)#p=0.03646***
lm.h.ct = lm(scale(Human_Structures)~Rich.prst, data=mod.env.ct)
summary(lm.h.ct)#p=0.5838

lm.h.nd = lm(Human_Structures~Int.prst, data=mod.env.nd)
summary(lm.h.nd)#p=0.7103
lm.h.nd = lm(Human_Structures~Rich.prst, data=mod.env.nd)
summary(lm.h.nd)#p=0.8598

#Transparency
lm.t.ct = lm(Secchi.cm.~Int.prst, data=mod.env.ct)
summary(lm.t.ct)#p=0.01535***
lm.t.ct = lm(Secchi.cm.~Rich.prst, data=mod.env.ct)
summary(lm.t.ct) #p=0.5781

lm.t.nd = lm(Secchi.cm.~Int.prst, data=mod.env.nd)
summary(lm.t.nd)#p=0.9751
lm.t.nd = lm(Secchi.cm.~Rich.prst, data=mod.env.nd)
summary(lm.t.nd) #p=0.5442

#Conductivity
lm.t.nd = lmp(Int.prst~log(Conductivity), data=mod.env.nd)
summary(lm.t.nd)#p=0.1711
lm.t.nd = lm(Rich.prst~log(Conductivity), data=mod.env.nd)
summary(lm.t.nd) #p=0.06217**


lm.t.ct = lm(Int.prst~log(Conductivity), data=mod.env.ct)
summary(lm.t.ct)#p=0.6363
lm.t.ct = lm(Rich.prst~log(Conductivity), data=mod.env.ct)
summary(lm.t.ct) #p=0.1524

#Modularity vs env. rich and abnd 
lm.h.ct = lm(scale(Human_Structures)~Z.value, data=mod.env.ct)
summary(lm.h.ct)#p=0.062***
lm.h.nd = lm(scale(Human_Structures)~Z.value, data=mod.env.nd)
summary(lm.h.nd) #p=0.85
lm.t.ct = lm(Z.value~ scale(Secchi.cm.)*scale(Human_Structures) + (1|Form), data=mod.env)
summary(lm.t.ct)#p=0.0562***
lm.t.nd = lm(scale(Secchi.cm.)~Z.value, data=mod.env.nd)
summary(lm.t.nd) #p=0.904
lm.i.ct = lm(Int.prst~Z.value, data=mod.env.ct)
summary(lm.i.ct) #p=0.068*
lm.r.ct = lm(Rich.prst~Z.value, data=mod.env.ct)
summary(lm.r.ct) #p=0.32
lm.i.nd = lm(Int.prst~Z.value, data=mod.env.nd)
summary(lm.i.nd) #p=0.01****
lm.r.nd = lm(Rich.prst~Z.value, data=mod.env.nd)
summary(lm.r.nd) #p=0.10

######################## Modules Identity#################
######Test if species are in the same module in each run####
randomised_ec = function(nets.graph.ct, limit){
  # create empty dataframe
  df = NULL
  # initialise counter
  randomisation_id = 100
  # randomise until limit
  while (randomisation_id - 1 < limit) {
    
    # create empty list
    modules3 <- igraph::cluster_leading_eigen(nets.graph.ct$`3`)
    MD3 <- igraph::modularity(modules3)
    mod3<-data.frame(Species= modules3$names, Modules= modules3$membership)
    betadiv_output = cbind(mod3, randomisation_id)
    df = rbind(df, betadiv_output)
    # update counter
    randomisation_id = randomisation_id + 1
  }
  return(df)
}

test4.ec= randomised_ec(nets.graph.ct,150)
str(test4.ec)
test4.ec$Modules=as.factor(test4.ec$Modules)


randomised_en = function(nets.graph.nd, limit){
  # create empty dataframe
  df = NULL
  # initialise counter
  randomisation_id = 50
  # randomise until limit
  while (randomisation_id - 1 < limit) {
    
    # create empty list
    modules3 <- igraph::cluster_leading_eigen(nets.graph.nd$`4`)
    MD3 <- igraph::modularity(modules3)
    mod3<-data.frame(Species= modules3$names, Modules= modules3$membership)
    betadiv_output = cbind(mod3, randomisation_id)
    df = rbind(df, betadiv_output)
    # update counter
    randomisation_id = randomisation_id + 1
  }
  return(df)
}

test4.en= randomised_en(nets.graph.nd,100)
test4.ec.u=test4.ec %>% select(-randomisation_id) %>% distinct()
test1.en.u=test1.en %>% select(-randomisation_id) %>% distinct()

######################## Make dataframe with module id for ecto ######################## 
#S1
modules1.ct <- igraph::cluster_optimal(nets.graph.ct$`1`)
mod1.ct<-data.frame(Species= modules1.ct$names, Modules= modules1.ct$membership)
m1.ct.h = mod1.ct[c(1:28),]
m1.ct.p =mod1.ct[c(29:39),]
#S2
modules2.ct <- igraph::cluster_optimal(nets.graph.ct$`2`)
mod2.ct<-data.frame(Species= modules2.ct$names, Modules= modules2.ct$membership)
m2.ct.h = mod2.ct[-c(24:32),]
m2.ct.p =mod2.ct[-c(1:23),]
#S3
modules3.ct <- igraph::cluster_optimal(nets.graph.ct$`3`)
mod3.ct<-data.frame(Species= modules3.ct$names, Modules= modules3.ct$membership)
m3.ct.h = mod3.ct[c(1:32),]
m3.ct.p =mod3.ct[c(33:52),]
#S4
modules4.ct <- igraph::cluster_optimal(nets.graph.ct$`4`)
mod4.ct<-data.frame(Species= modules4.ct$names, Modules= modules4.ct$membership)
m4.ct.h = mod4.ct[c(1:39),]
m4.ct.p =mod4.ct[c(40:66),]

####Include modules id and czvalues in traits dataframe
m1.ct.h<- m1.ct.h %>% rename("COD.IND" = Species) 
m1.ct.h=left_join(m1.ct.h, host.traits, by= "COD.IND")
m1.ct.h=left_join(m1.ct.h, kn, by= "COD.IND")

mod1.p.ct<- mod1.p.ct %>% rename("Species" = Col_names) 
m1.ct.p=left_join(m1.ct.p, parasite.traits, by= "Species")

m2.ct.h<- m2.ct.h %>% rename("COD.IND" = Species) 
m2.ct.h=left_join(m2.ct.h, host.traits, by= "COD.IND")
m2.ct.h=left_join(m2.ct.h, kn, by= "COD.IND")

mod2.p.ct<- mod2.p.ct %>% rename("Species" = Col_names) 
m2.ct.p=left_join(m2.ct.p, parasite.traits, by= "Species")

m3.ct.h<- m3.ct.h %>% rename("COD.IND" = Species) 
m3.ct.h=left_join(m3.ct.h, host.traits, by= "COD.IND")
m3.ct.h=left_join(m3.ct.h, kn, by= "COD.IND")

mod3.p.ct<- mod3.p.ct %>% rename("Species" = Col_names) 
m3.ct.p=left_join(m3.ct.p, parasite.traits, by= "Species")

m4.ct.h<- m4.ct.h %>% rename("COD.IND" = Species) 
m4.ct.h=left_join(m4.ct.h, host.traits, by= "COD.IND")
m4.ct.h=left_join(m4.ct.h, kn, by= "COD.IND")

mod4.p.ct<- mod4.p.ct %>% rename("Species" = Col_names) 
m4.ct.p=left_join(m4.ct.p, parasite.traits, by= "Species")


##### Make dataframe with module id for endo##################
#S1
modules1.nd <- igraph::cluster_optimal(nets.graph.nd$`1`)
mod1.nd<-data.frame(Species= modules1.nd$names, Modules= modules1.nd$membership)
m1.nd.h = mod1.nd[c(1:21),]
m1.nd.p =mod1.nd[c(22:31),]
#S2
modules2.nd <- igraph::cluster_optimal(nets.graph.nd$`2`)
mod2.nd<-data.frame(Species= modules2.nd$names, Modules= modules2.nd$membership)
m2.nd.h = mod2.nd[c(1:27),]
m2.nd.p =mod2.nd[c(28:40),]
#S3
modules3.nd <- igraph::cluster_optimal(nets.graph.nd$`3`)
mod3.nd<-data.frame(Species= modules3.nd$names, Modules= modules3.nd$membership)
m3.nd.h = mod3.nd[c(1:18),]
m3.nd.p =mod3.nd[c(19:27),]
#S4
modules4.nd <- igraph::cluster_optimal(nets.graph.nd$`4`)
mod4.nd<-data.frame(Species= modules4.nd$names, Modules= modules4.nd$membership)
m4.nd.h = mod4.nd[c(1:27),]
m4.nd.p =mod4.nd[c(28:46),]

####Include modules id in traits dataframe
m1.nd.h<- m1.nd.h %>% rename("COD.IND" = Species) 
m1.nd.h=left_join(m1.nd.h, host.traits, by= "COD.IND")
m1.nd.h=left_join(m1.nd.h, kn, by= "COD.IND")

mod1.p.nd<- mod1.p.nd %>% rename("Species" = Col_names) 
m1.nd.p=left_join(m1.nd.p, parasite.traits, by= "Species")

m2.nd.h<- m2.nd.h %>% rename("COD.IND" = Species) 
m2.nd.h=left_join(m2.nd.h, host.traits, by= "COD.IND")
m2.nd.h=left_join(m2.nd.h, kn, by= "COD.IND")

mod2.p.nd<- mod2.p.nd %>% rename("Species" = Col_names) 
m2.nd.p=left_join(m2.nd.p, parasite.traits, by= "Species")

m3.nd.h<- m3.nd.h %>% rename("COD.IND" = Species) 
m3.nd.h=left_join(m3.nd.h, host.traits, by= "COD.IND")
m3.nd.h=left_join(m3.nd.h, kn, by= "COD.IND")

mod3.p.nd<- mod3.p.nd %>% rename("Species" = Col_names) 
m3.nd.p=left_join(m3.nd.p, parasite.traits, by= "Species")

m4.nd.h<- m4.nd.h %>% rename("COD.IND" = Species) 
m4.nd.h=left_join(m4.nd.h, host.traits, by= "COD.IND")
m4.nd.h=left_join(m4.nd.h, kn, by= "COD.IND")

mod4.p.nd<- mod4.p.nd %>% rename("Species" = Col_names) 
m4.nd.p=left_join(m4.nd.p, parasite.traits, by= "Species")

#####Write results#####################

write.csv(m1.nd.h, "T1.hsp.end.csv")
write.csv(m2.nd.h, "T2.hsp.end.csv")
write.csv(m3.nd.h, "T3.hsp.end.csv")
write.csv(m4.nd.h, "T4.hsp.end.csv")

write.csv(m1.nd.p , "T.p.mod1nd.csv")
write.csv(m2.nd.p , "T.p.mod2nd.csv")
write.csv(m3.nd.p , "T.p.mod3nd.csv")
write.csv(m4.nd.p , "T.p.mod4nd.csv")

#################### Mean values for host traits###############
T1.hsp.end$Modules = as.factor(T1.hsp.end$Modules)
meant1nd = T1.hsp.end %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.endo) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant2nd = T2.hsp.end %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.endo) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant3nd = T3.hsp.end %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.endo) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant4nd = T4.hsp.end %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.endo) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant1ct = T1.hsp.ect %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.inf) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant2ct = T2.hsp.ect %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.inf) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant3ct = T3.hsp.ect %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.inf) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

meant4ct = T4.hsp.ect %>% dplyr::select(Modules, COMP.P.cm., PESO.g., Kn, Int.inf) %>% group_by(Modules) %>% summarise_all(.funs= funs(mean, sd), na.rm=T)

write.csv(meant1nd, "meant1nd.csv")
write.csv(meant2nd, "meant2nd.csv")
write.csv(meant3nd, "meant3nd.csv")
write.csv(meant4nd, "meant4nd.csv")

write.csv(meant1ct, "meant1ct.csv")
write.csv(meant2ct, "meant2ct.csv")
write.csv(meant3ct, "meant3ct.csv")
write.csv(meant4ct, "meant4ct.csv")
