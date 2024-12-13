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
var2<-c("SLA","%C","%N","C:N","Vmax","Vpmax","Vpmax/Vmax","Leaf area","%C","d13C","d15N","SL","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","iWUE\n sat light\n 400 CO2","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
        "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n recovery\n offset","NPQ\n induction\n slope",
        "Maximum\n NPQ","Final\n NPQ","Final\n phiPSII","Fv/Fm")


# yearsused=list()
# yearsused[[1]]=2021
# yearsused[[2]]=2022
# yearsused[[3]]=2023
# yearsused[[4]]=c(2021,2022)
# yearsused[[5]]=c(2021,2022,2023)
yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=2023
yearsused[[4]]=c(2021,2022)
yearsused[[5]]=c(2021,2023)
yearsused[[6]]=c(2022,2023)
yearsused[[7]]=c(2021,2022,2023)
types=c("none","plot","acc") #c("none","plot","acc")
typesp=c("Raw data","Averaged by plot","Averaged by genotype") #,"accession"
scenarios<-c("Scenario i","Scenario ii","Scenario iii")
files=c("CN_isotopes","ACi","SLA","lightSaturatedGasExchange") #

trait_types<-c("Structural & biochemical", "Gas exchange", "Chlorophyll fluorescence")

vars<-list()
vars[[1]]<-c("SLA","N","CN")
vars[[2]]<-c("Vpmax","Vmax","a400","gsw")
vars[[3]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

scenas=c(1)

methods=c("SVR")
## Plot
# rm(dfscore)
for (scena in scenas){
  for (t in 2:2){
    aggtype=types[t]
    for (m in 1:length(methods)){
      for (y in 1:7){
        yearsp=c()
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
        yearsp2<-gsub("_","&",yearsp)
        for (f in 1:length(vars)){
          for (varname in vars[[f]]){
            filen=paste0("results/scenario",scena,"/","maize","_",methods[m],"_", varname,"_",yearsp,"_",aggtype,".RData")
            if (file.exists(filen)){
              load(filen)
              dat=melt(testing_R2)
              # dat=melt(testing_corr2)
              # df=data.frame(Score=dat$value,Variable=varname,Year=as.character(yearsp),Aggregation=typesp[t])
              newvar<-var2[match(varname,var1)]
              
              
              df=data.frame(Score=dat$value,Variable=newvar,Year=as.character(yearsp2),Method=methods[m],Aggregation=typesp[t],Trait_type=trait_types[f],Scenario=scenarios[scena])
              
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


dfscore$Variable<-factor(dfscore$Variable, levels=unique(dfscore$Variable))
dfscore$Trait_type<-factor(dfscore$Trait_type, levels=trait_types)
dfscore$Aggregation<-factor(dfscore$Aggregation,levels=typesp)
dfscore$Year<-factor(dfscore$Year, levels=c("2021","2022","2023","2021&2022","2021&2023","2022&2023","2021&2022&2023"))


g1<-ggplot(dfscore[dfscore$Trait_type=="Structural & biochemical",], aes(x=Variable, y=Score,fill=Year)) +
  geom_boxplot()+
  scale_y_continuous(limits = c(-1, 0.9) )+ # position = "right" ,
  labs(y="",x="",fill="Season") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Year), size=4,position=position_dodge(0.9)) +
  coord_cartesian(ylim =c(-0.5, 1.3))+ #min(dfscore$Score)
  # scale_fill_brewer() +
  theme(legend.position = "none",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 15,face="bold"),
        strip.text = element_text(size=19,face="bold"),axis.title=element_text(size = 19,face="bold"))+
  scale_fill_brewer(palette = "Dark2")

g2<-ggplot(dfscore[dfscore$Trait_type=="Gas exchange",], aes(x=Variable, y=Score,fill=Year)) +
  geom_boxplot()+
  # cale_y_continuous(limits = c(-1, 0.9) )+ # position = "right" ,
  labs(y="Coefficient of determination",x="",fill="Season") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Year), size=4,position=position_dodge(0.9)) +
  coord_cartesian(ylim =c(-1, 1.3))+ #min(dfscore$Score)
  # scale_fill_brewer() +
  theme(legend.position = "none",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 15,face="bold"),
        strip.text = element_text(size=19,face="bold"),axis.title=element_text(size = 17,face="bold"))+
  scale_fill_brewer(palette = "Dark2")


g3<-ggplot(dfscore[dfscore$Trait_type=="Chlorophyll fluorescence",], aes(x=Variable, y=Score,fill=Year)) +
  geom_boxplot()+
  labs(y="",x="",fill="Season") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Year), size=4,position=position_dodge(0.9)) +
  coord_cartesian(ylim =c(-1, 1.3))+ #min(dfscore$Score)
  # scale_fill_brewer() +
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 15,face="bold"),
        strip.text = element_text(size=19,face="bold"),axis.title=element_text(size = 19,face="bold"))+
  scale_fill_brewer(palette = "Dark2")


plot_grid(g1, g2, g3, hjust=-2,
          labels = c("a","b","c"), label_size=17,
          # labels = c("Anatomical","Gas exchange","Chlorophyll fluorescence"),
          nrow =3,rel_heights = c(0.8,0.8,1.1),axis= "rlbt")


ggsave("figures/SF5_yearscomb_plot_averaged_SVR.png",width = 13, height = 13)

