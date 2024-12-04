library(dplyr)
library(tidyr)
library(ggplot2)
library(pls)
library(doParallel)
library(caret)
library(reshape2)

# setwd("/work/xu2/HyperspectralML/")
setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")

### Result Section 1 ###

# traits=list()
# years<-c(2021,2022,2023)
# for (y in 1:3){
#   year=years[y]
#   ini<-T
#   for (filen in c("SLA","CN_isotopes","fluorescence","ACi")){
#     Xdata<-load_traits("maize",year,filen,"acc",FALSE)
#     
#     ind<-which(is.na(Xdata),arr.ind = T)
#     if (dim(ind)[1]!=0){
#       Xdata<-Xdata[-ind[,2],]
#     }
#     ind<-which(colnames(Xdata)%in%"Repeat")
#     if (length(ind)>0){
#       Xdata<-Xdata[,-ind]
#     }
#     if (ini){
#       Xdatatot<-Xdata
#       ini<-F
#     }else{
#       Xdatatot<-merge(Xdatatot,Xdata)
#     }
#   }
#   traits[[y]]<-Xdatatot
# }



vars<-c("SLA","C",  "N",  "CN", "d13C","d15N","NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
"NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
"maxNPQ","endNPQ","endFvFm","initialFvFm","Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE" )

year_comb=matrix(c(2021,2022,2021,2023,2022,2023),3,2,byrow = T)
cormat<-matrix(0,26,3)
ncount<-matrix(0,26,3)
rownames(cormat)<-vars

for (c in 1:nrow(year_comb)){
  for (filen in c("SLA","CN_isotopes","fluorescence","ACi")){
    trait1=load_traits("maize",year_comb[c,1],filen,"acc",FALSE)
    trait2=load_traits("maize",year_comb[c,2],filen,"acc",FALSE)
    
    if (year_comb[c,1]==2021){
      for (i in 1:length(trait1$accession)){
        temp<-strsplit(trait1$accession[i],'_')
        trait1$accession[i]<-paste0(temp[[1]][1],'_00',temp[[1]][2])
      }
    }
    
    if (year_comb[c,2]==2021){
      for (i in 1:length(trait2$accession)){
        temp<-strsplit(trait2$accession[i],'_')
        trait2$accession[i]<-paste0(temp[[1]][1],'_00',temp[[1]][2])
      }
    }
    commonacc<-intersect(trait1$accession,trait2$accession)
    temp1<-trait1[trait1$accession%in%commonacc,]
    temp2<-trait2[trait2$accession%in%commonacc,]
    
    ind<-which(colnames(temp1)%in%c("accession","Repeat"))
    if (length(ind)>0){
      temp1<-temp1[,-ind]
    }
    ind<-which(colnames(temp2)%in%c("accession","Repeat"))
    if (length(ind)>0){
      temp2<-temp2[,-ind]
    }
    
    var1<-colnames(temp1)
    temp2<-temp2[,match(var1,colnames(temp2))]
    
    for (v in 1:length(var1)){
      indvar<-match(var1[v],vars)
      cormat[indvar,c]<-cor(temp1[,v],temp2[,v],method="spearman")
      ncount[indvar,c]<-nrow(temp1)
    }
  }
}

load(file="correlation_traits_spearman.RData")


############# Count the trait numbers ###############

types=c("none","plot","acc") 
years<-c(2021,2022,2023)
files<-c("SLA","CN_isotopes","fluorescence","ACi")
ncount<-array(0,c(4,3,3))
rownames(ncount)<-files

for (g in 1:3){
  for (y in 1:3){
    for (f in 1:4){
      filen<-files[f]
      trait=load_traits("maize",years[y],filen,types[g],FALSE)
      ncount[f,y,g]<-nrow(trait)
    }
  }
}

