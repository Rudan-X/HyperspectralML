library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pls)
library(caret)
library(readxl)
library(e1071)
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
var2<-c("SLA","%C","%N","C/N","Vmax","Vpmax","Vpmax/Vmax","Leaf area","%C","d13C","d15N","SL","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","iWUE\n sat light\n 400 CO2","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
        "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n induction\n residual","NPQ\n induction\n slope",
        "Maximum\n NPQ","Final\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")



yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=c(2021,2022)
yearsused[[4]]=2023
yearsused[[5]]=c(2021,2022,2023)
types=c("none","plot","acc") 
typesp=c("Raw data","Plot-averaged","Genotype-averaged") 
scenarios<-c("Scenario i\n(R²)","Scenario ii\n(R²)","Scenario iii\n(Cor²)")
files=c("CN_isotopes","ACi","SLA","lightSaturatedGasExchange") #

trait_types<-c("Structural & biochemical", "Photosynthesis-related", "Chlorophyll fluorescence")

vars<-list()
vars[[1]]<-c("d15N","C","d13C")
vars[[2]]<-c("Vpmax.Vmax","SL","iWUE")#,"gs_sat","ci_sat","iWUE_sat"
vars[[3]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","NPQ_ind_linear","endNPQ")

scenas=c(1)

methods=c("PLSR","SVR")#"PLSR"
## Plot

for (scena in scenas){
  for (t in 1:2){
    aggtype=types[t]
    for (m in 1:length(methods)){
      for (y in 1:1){
        years<-yearsused[[y]]
        if (length(years)>1){
          yearsp=c()
          for (Y in years){
            if (Y==years[1]){
              yearsp=Y
            }else{yearsp=paste0(yearsp,"_",Y)}
          }
        }else{
          yearsp=years
        }
        for (f in 1:length(vars)){
          for (varname in vars[[f]]){
            
            filen=paste0("results/scenario",scena,"/","maize","_",methods[m],"_", varname,"_",yearsp,"_",aggtype,".RData")
            if (file.exists(filen)){
              load(filen)
              dat=melt(testing_R2)
              # dat=melt(testing_corr2)
              # df=data.frame(Score=dat$value,Variable=varname,Year=as.character(yearsp),Aggregation=typesp[t])
              newvar<-var2[match(varname,var1)]
              df=data.frame(Score=dat$value,Variable=newvar,Year=as.character(yearsp),Method=methods[m],Aggregation=typesp[t],Trait_type=trait_types[f],Scenario=scenarios[scena])
              
              if (!exists("dfscore")){
                dfscore=df
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


#############

dfscore$Variable<-factor(dfscore$Variable, levels=unique(dfscore$Variable))
dfscore$Trait_type<-factor(dfscore$Trait_type, levels=trait_types)
dfscore$Aggregation<-factor(dfscore$Aggregation, levels=unique(dfscore$Aggregation))
# dfscore$Year<-factor(dfscore$Year, levels=c("2021","2022","2023","2021_2022","2021_2022_2023"))


ggplot(dfscore, aes(x=Variable, y=Score,fill=Method)) +
  geom_boxplot()+
  facet_grid(Aggregation ~ Trait_type,scales="free", switch = "y",space="free")+
  scale_y_continuous(position = "right" ,limits = c(-1, 0.9) )+ #
  labs(y="Coefficient of determination",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Method), size=5,position=position_dodge(0.9)) +
  coord_cartesian(ylim =c(-1, 1.3))+ #min(dfscore$Score)
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 13,face="bold"),
        axis.title.y=element_text(size=18,face="bold"),
        strip.text = element_text(size=17,face="bold"))+
  stat_compare_means(label = "p.signif",label.y.npc = "top",vjust=-2) +
  scale_fill_brewer(palette = "Dark2")


ggsave("figures/SF3.png",width = 15, height = 8)
