library(dplyr)
library(tidyr)
library(ggplot2)
library(pls)
library(doParallel)
library(caret)
library(reshape2)
library(vegan)
library(geosphere)
# setwd("/work/xu2/HyperspectralML/")
setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")

ref=list()
years<-c(2021,2022,2023)
for (y in 1:3){
  year<-years[y]
  filen <- paste("data/reflectance/", year, "_hyperspectral_reflectance_maize.csv", sep = "")
  Xdata <- read.csv(filen)
  ind<-which(is.na(Xdata),arr.ind = T)
  if (dim(ind)[1]!=0){
    Xdata<-Xdata[-ind[,2],]
  }
  acc_ref <- read.csv(paste0("data/accessions/",year,"_maize_reference.csv"))
  ind<-match(Xdata$Plot,acc_ref$Plot)
  Xdata$Plot<-acc_ref$Accession[ind]
  colnames(Xdata)[1]<-"Accession"
  data <- Xdata %>%
    group_by(Accession) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  ref[[y]]<-data
}


for (i in 1:length(ref[[1]]$Accession)){
  temp<-strsplit(ref[[1]]$Accession[i],'_')
  ref[[1]]$Accession[i]<-paste0(temp[[1]][1],'_00',temp[[1]][2])
}


############## Mantel correlation ##################
year_comb=matrix(c(1,2,1,3,2,3),3,2,byrow = T)
mantel_cor<-list()
for (c in 1:nrow(year_comb)){
  ind1=year_comb[c,1]
  ind2=year_comb[c,2]
  
  commonacc<-intersect(ref[[ind1]]$Accession,ref[[ind2]]$Accession)
  
  temp1<-ref[[ind1]][ref[[ind1]]$Accession%in%commonacc,]
  temp2<-ref[[ind2]][ref[[ind2]]$Accession%in%commonacc,]
  
  covar.mat1<-cov(temp1[,3:ncol(temp1)])
  covar.mat2<-cov(temp2[,3:ncol(temp1)])
  
  mantel_cor[[c]]<-mantel(covar.mat1, covar.mat2, method = "spearman", permutations = 50, na.rm = TRUE)
}



save(mantel_cor,file="./results/tables/HSR_mantel_corr.RData")




###################################################

year_comb=matrix(c(1,2,1,3,2,3),3,2,byrow = T)
cormat<-matrix(0,3,(ncol(ref[[1]])-2))
for (c in 1:nrow(year_comb)){
  ind1=year_comb[c,1]
  ind2=year_comb[c,2]
  
  commonacc<-intersect(ref[[ind1]]$Accession,ref[[ind2]]$Accession)
  temp1<-ref[[ind1]][ref[[ind1]]$Accession%in%commonacc,]
  temp2<-ref[[ind2]][ref[[ind2]]$Accession%in%commonacc,]
  

  mean1<-colMeans(temp1[,3:ncol(temp1)],na.rm = T)
  mean2<-colMeans(temp2[,3:ncol(temp2)],na.rm = T)
  
  for (w in 1:(ncol(ref[[1]])-2)){
    cormat[c,w]<-cor(temp1[,2+w],temp2[,2+w])
  }
}
colnames(cormat)<-colnames(temp1)[3:ncol(temp1)]

## Correlation between 2022&2021 and 2022&2023 for diff wv

hstart<-which(colnames(cormat)%in%paste0("mean_",350))
hend<-which(colnames(cormat)%in%paste0("mean_",700))

mean(cormat[,hstart:hend])
mean(cormat[1,hstart:hend])
mean(cormat[2,hstart:hend])
mean(cormat[3,hstart:hend])


hstart<-which(colnames(cormat)%in%paste0("mean_",700))
hend<-which(colnames(cormat)%in%paste0("mean_",1300))

mean(cormat[1,hstart:hend])
mean(cormat[2,hstart:hend])
mean(cormat[3,hstart:hend])

hstart<-which(colnames(cormat)%in%paste0("mean_",1300))
hend<-which(colnames(cormat)%in%paste0("mean_",1900))

mean(cormat[,hstart:hend])
mean(cormat[2,hstart:hend])
mean(cormat[3,hstart:hend])

hstart<-which(colnames(cormat)%in%paste0("mean_",1300))
hend<-which(colnames(cormat)%in%paste0("mean_",2500))

mean(cormat[,hstart:hend])
mean(cormat[2,hstart:hend])
mean(cormat[3,hstart:hend])#


hstart<-which(colnames(cormat)%in%paste0("mean_",350))
hend<-which(colnames(cormat)%in%paste0("mean_",2500))
mean(cormat[1,hstart:hend])
mean(cormat[2,hstart:hend])
mean(cormat[3,hstart:hend])


### Result Section 1 ###

traits=list()
years<-c(2021,2022,2023)
for (y in 1:3){
  year=years[y]
  ini<-T
  for (filen in c("SLA","CN_isotopes","fluorescence","ACi")){
    Xdata<-load_traits("maize",year,filen,"acc",FALSE)
    
    ind<-which(is.na(Xdata),arr.ind = T)
    if (dim(ind)[1]!=0){
      Xdata<-Xdata[-ind[,2],]
    }
    ind<-which(colnames(Xdata)%in%"Repeat")
    if (length(ind)>0){
      Xdata<-Xdata[,-ind]
    }
    if (ini){
      Xdatatot<-Xdata
      ini<-F
    }else{
      Xdatatot<-merge(Xdatatot,Xdata)
    }
  }
  traits[[y]]<-Xdatatot
}


for (i in 1:length(traits[[1]]$accession)){
  temp<-strsplit(traits[[1]]$accession[i],'_')
  traits[[1]]$accession[i]<-paste0(temp[[1]][1],'_00',temp[[1]][2])
}

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

write.csv(cormat,file="./results/tables/trait_corr_acc_averaged.csv")
load(file="correlation_traits_spearman.RData")

# cormatspearman=cormat
# load(file="correlation_traits.RData")


# save(cormat,ncount,file="correlation_traits_spearman.RData")

df<-melt(cormat)
colnames(df)<-c("Variable","Pair","Correlation")
df$Pair[df$Pair==1]="2021 vs 2022"
df$Pair[df$Pair==2]="2021 vs 2023"
df$Pair[df$Pair==3]="2022 vs 2023"

ggplot(df, aes(x=Pair, y=Correlation)) +
  geom_boxplot() + stat_compare_means(label = "p.signif",label.y.npc = "top")

my_comparisons <- list( c("2021", "2022"), c("2021", "2023"), c("2022", "2023"))


diff<-cor1-cor2
