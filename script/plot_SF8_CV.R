library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)
library(rstatix)

rm(list = ls())
setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")

var1<-c("SLA","C","N","CN","Vmax","Vpmax","Vpmax.Vmax","Area","C","d13C","d15N","SL","a400","gsw","iWUE","A_sat","gs_sat","ci_sat","iWUE_sat",
        "gs","sl",
        "NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res",
        "phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear",
        "maxNPQ","endNPQ","endFvFm","initialFvFm")
var2<-c("SLA","%C","%N","C:N","Vmax","Vpmax","Vpmax/Vmax","Leaf area","%C","d13C","d15N","SL","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","iWUE\n sat light\n 400 CO2","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
        "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n recovery\n offset","NPQ\n induction\n slope",
        "Maximum\n NPQ","Final\n NPQ","Final\n phiPSII","Fv/Fm")



yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=c(2021,2022)
yearsused[[4]]=2023
yearsused[[5]]=c(2021,2022,2023)
types=c("none","plot","acc") 
typesp=c("Raw data","Plot-averaged","Genotype-averaged") 
scenarios<-c("a. Seen genotype\n(R²)","b. Unseen genotype\n(R²)","c. Unseen year\n(Cor²)")
files=c("CN_isotopes","ACi","SLA","lightSaturatedGasExchange") #

trait_types<-c("Structural & biochemical traits", "Photosynthetic traits", "Chlorophyll fluorescence")

vars<-list()
vars[[1]]<-c("SLA")
vars[[2]]<-c("N","CN")
vars[[3]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
vars[[4]]<-c("Vpmax","Vmax","a400","gsw")


files<-c("SLA","CN_isotopes","fluorescence","ACi")

variables<-c("SLA","N","CN","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm","Vpmax","Vmax","a400","gsw")

types=c("none","plot","acc") 
type_lab=c("Raw data","Plot-averaged","Genotype-averaged")
ini<-T
for (t in 1:3){
  aggtype=types[t]
  cv_mat<-matrix(0,100,length(variables))
  colnames(cv_mat)<-var2[match(variables,var1)]
  
  for (f in 1:4){
    filen<-files[f]
    data<-get_XYdata("maize",2021, filen,aggtype) 
    ind=which(colnames(data)=="Repeat")
    if (length(ind)>0){
      data<-data[,-ind]
    }
    data0<-data
    for (varname in vars[[f]]){
      data<-data0
      accessions<-unique(data$Accession)
      acctot<-data$Accession
      yind<-which(colnames(data)%in%varname)
      data<-data[,yind]
      data<-as.matrix(data)
      nrep<-20
      
      for (r in 1:nrep){
        accfolds <- mycreateFolds(accessions, k = 5)
        for (j in 1:length(accfolds)){
          test_acc<-accessions[accfolds[[j]]]
          train_acc<-accessions[setdiff(1:length(accessions),accfolds[[j]])]
          test_ind<-which(acctot%in%test_acc)
          train_ind<-setdiff(1:nrow(data),test_ind)
          
          data_train<-data[train_ind]
          cv_mat[(r-1)*5+j,match(varname,variables)]<-sd(data_train)/mean(data_train)*100
        }
      }
    }
  }
  df<-melt(cv_mat)
  df<-df[,-1]
  colnames(df)<-c("Variable","CV")
  df$Aggregation<-type_lab[t]
  if (ini){
    dfCV<-df
    ini<-F
  }else{dfCV<-rbind(dfCV,df)}
}


##

trait_types<-list()
trait_types[[1]]<-c("SLA","%C","%N","C:N")
trait_types[[2]]<-c("Vmax","Vpmax","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2")
trait_types[[3]]<-c("phiPSII\n recovery\n offset","Maximum\n NPQ","Final\n phiPSII","Fv/Fm")

dfCV$Trait_type<-"Structural & biochemical"
dfCV$Trait_type[dfCV$Variable%in%trait_types[[2]]]<-"Gas exchange"
dfCV$Trait_type[dfCV$Variable%in%trait_types[[3]]]<-"Chlorophyll fluorescence"

dfCV$Variable<-factor(dfCV$Variable, levels=unique(dfCV$Variable))
dfCV$Aggregation<-factor(dfCV$Aggregation, levels=type_lab)
dfCV$Trait_type<-factor(dfCV$Trait_type, levels=c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))

##
save(dfCV,file="CV_scenario2_3aggtypes.RData")

load(file="CV_scenario2_3aggtypes.RData")

ggplot(dfCV, aes(x=Variable, y=CV,fill=Aggregation)) +
  geom_boxplot()+
  facet_grid( ~ Trait_type,scales="free",space="free")+
  labs(y="Coefficient of variation (%)",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Aggregation), size=4,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 15,face="bold"),
        strip.text = element_text(size=19,face="bold"),axis.title=element_text(size = 19,face="bold"))+
  # stat_compare_means(label = "p.signif",label.y.npc = "top",vjust=-2) +
  scale_fill_brewer(palette = "Dark2")


ggsave("figures/SF8_CV.png",width = 12, height = 7)

