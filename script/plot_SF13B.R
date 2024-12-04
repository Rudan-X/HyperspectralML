library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)
rm(list = ls())
setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")

var1<-c("SLA","C","N","CN","Vmax","Vpmax","Vpmax.Vmax","Area","C","d13C","d15N","SL","a400","gsw","iWUE","A_sat","gs_sat","ci_sat","iWUE_sat",
        "gs","sl",
        "NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res",
        "phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear",
        "maxNPQ","endNPQ","endFvFm","initialFvFm")
var2<-c("SLA","%C","%N","C/N","Vmax","Vpmax","Vpmax/Vmax","Leaf area","%C","d13C","d15N","SL","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","iWUE\n sat light\n 400 CO2","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
        "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n induction\n residual","NPQ\n induction\n slope",
        "Maximum\n NPQ","Final\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")



types=c("none","plot","acc") 
typesp=c("Raw data","Plot-averaged","Genotype-averaged") 
scenarios<-c("Scenario i\n(R²)","Scenario ii\n(R²)","Scenario iii\n(Cor²)")
files=c("CN_isotopes","ACi","SLA","lightSaturatedGasExchange") #

trait_types<-c("Structural & biochemical", "Photosynthesis-related", "Chlorophyll fluorescence")

vars<-list()
vars[[1]]<-c("SLA","N","CN")
vars[[2]]<-c("Vpmax","Vmax","a400")
vars[[3]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")


methods=c("PLSR")#"PLSR"

yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=2023

genotypes<-c("seen","unseen")

ini<-T
aggtype=types[2]
for (m in 1:length(methods)){
  for (y in 1:1){
    testing_year<-yearsused[[y]]
    
    years=list()
    years[[1]]=2021
    years[[2]]=2022
    years[[3]]=2023
    years[[4]]=c(2021,2022)
    years[[5]]=c(2021,2023)
    years[[6]]=c(2022,2023)
    
    check=c()
    for (j in years){
      if (sum(testing_year==j)==0){
        check=c(check,TRUE)
      }else{check=c(check,FALSE)}
    }
    training_years<-years[check]
    
    yearsp<-c()
    for(k in 1:length(training_years)){
      temp<-training_years[[k]]
      check<-temp[1]
      if (length(temp)>1){
        for (j in 2:length(temp)){
          check<-paste0(check,"&",temp[j])
        }
      }
      yearsp[k]=check
    }
    
    
    for (f in 1:3){
      for (varname in vars[[f]]){
        for (r in 1:3){
          for (g in 1:2){
            filen=paste0("results/scenario3/",genotypes[g],"_maize_",methods[m],"_", varname,"_testY",testing_year,"_trainY",yearsp[r],"_",aggtype,".RData")
            if (file.exists(filen)){
              load(filen)
              dat=melt(testing_corr)
              newvar<-var2[match(varname,var1)]
              df=data.frame(Correlation=dat$value^2,Variable=newvar,Trait_type=trait_types[f],Genotype=genotypes[g],Training_years=paste0(yearsp[r]))
              
              if (ini){
                dfscore<-df
                ini<-F
              }else{
                dfscore=rbind(dfscore,df)
              }
            }
          }
          
        }
        
      }
    }
  }
}


dfscore$Variable<-factor(dfscore$Variable, levels=unique(dfscore$Variable))
dfscore$Training_years<-factor(dfscore$Training_years, levels=yearsp)
dfscore$Trait_type<-factor(dfscore$Trait_type, levels=trait_types)
# dfscore$Year<-factor(dfscore$Year, levels=c("2021","2022","2023","2021_2022","2021_2022_2023"))


ggplot(dfscore, aes(x=Variable, y=Correlation,fill=Training_years)) +
  geom_boxplot()+
  facet_grid(Genotype ~ Trait_type,scales="free")+
  scale_y_continuous(limits = c(-1, 0.9) )+ #
  labs(y="Correlation²",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Training_years), size=4,position=position_dodge(0.9)) +
  coord_cartesian(ylim =c(-1, 1.3))+ #min(dfscore$Score)
  # scale_fill_brewer() + #legend.position = c(.93, .95)
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 15,face="bold"),axis.title.y=element_text(size=18,face="bold"),
        strip.text = element_text(size=19,face="bold"))+
  # stat_compare_means(label = "p.signif",label.y.npc = "top",vjust=-2) +
  scale_fill_brewer(palette = "Dark2")


ggsave("figures/SF13.png",width = 16, height = 10)
