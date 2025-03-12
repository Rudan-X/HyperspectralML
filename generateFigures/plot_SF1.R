library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(reshape2)
library(cowplot)


setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")


files=c("SLA","ACi","CN_isotopes","fluorescence")
trait_types<-c("Structural & biochemical", "Gas exchange", "Structural & biochemical","Chlorophyll fluorescence")
ini<-TRUE
for (f in 1:length(files)){
  filen<-files[f]
  if (filen=="lightSaturatedGasExchange"){
    years=c(2021,2022)
  }else if (filen=="leaf_width"){
    years=c(2022,2023)
  }else{years=c(2021,2022,2023)}
  for (year in years){
    data <- load_traits("maize",year,filen,"none",FALSE)
    data<-data[complete.cases(data),]
    ind<-which(colnames(data)%in%c("Plot","Repeat"))
    data<-data[,-ind]
    
    if (filen=="CN_isotopes" & year==2022){
      data<-data[data$C!=0,]
    }
    
    df=melt(data)
    df$year=year
    df$trait=trait_types[f]
    ind<-which(colnames(df)%in%"variable")
    if (length(ind)==0){
      df$variable=filen
    }
    if (ini){
      dataf=df
      ini<-FALSE
    }else{
      dataf=rbind(dataf,df)
    }
  }
} 


variables<-unique(dataf$variable)
years<-c(2021,2022,2023)
cvmat<-matrix(0,length(variables),length(years))
for (y in 1:3){
  for (v in 1:length(variables)){
    var<-variables[v]
    data<-dataf$value[dataf$variable==var & dataf$year==years[y]]
    cvmat[v,y]<-sd(data) / mean(data) * 100
  }  
}
rownames(cvmat)<-variables

cvaverage<-rowMeans(cvmat)

var1<-c("SLA","C","N","CN","Vmax","Vpmax","Vpmax.Vmax","Area","C","d13C","d15N","SL","a400","gsw","iWUE","A_sat","gs_sat","ci_sat","iWUE_sat",
        "gs","sl",
        "NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res",
        "phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear",
        "maxNPQ","endNPQ","endFvFm","initialFvFm")
var2<-c("SLA","%C","%N","C:N","Vmax","Vpmax","Vpmax/Vmax","Leaf area","%C","d13C","d15N","Stomatal\n limitation","Photosynthesis\n sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","iWUE\n sat light\n 400 CO2","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
        "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n recovery\n offset","NPQ\n induction\n slope",
        "Maximum\n NPQ","Final\n NPQ","Final\n phiPSII","Fv/Fm")

ind<-match(dataf$variable,var1)
dataf$variable<-var2[ind]

dataf$year=as.factor(dataf$year)
dataf$variable<-factor(dataf$variable,levels=unique(dataf$variable))

# df=df[df$variable=="Vmax",]
# ggplot(dataf, aes(x=variable, y=value,fill=year)) +
#   geom_violin() + facet_wrap(. ~ variable,scales="free",ncol = 4) + 
#   theme(legend.text = element_text(size=14,face="plain"),legend.title=element_text(size=14,face="bold"),
#         text = element_text(size = 12,face="bold"),axis.text =  element_text(size = 13,face="bold"),
#         axis.title=element_text(size=15,face="bold"),
#         strip.text = element_text(size=13,face="bold"))+
#   labs(x="Traits",y="Measured value",fill="Season") +scale_fill_grey()

  
g1<-ggplot(dataf[dataf$trait=="Structural & biochemical",], aes(x=variable, y=value,fill=year)) +
  geom_violin() + facet_wrap(. ~ variable,scales="free",ncol = 6) +
  theme(legend.position="none",legend.text = element_text(size=14,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 12,face="bold"), axis.text.x=element_blank(),
        strip.text = element_text(size=11,face="bold"))+
  labs(x="",y="",fill="Season") +scale_fill_grey()

g2<-ggplot(dataf[dataf$trait=="Gas exchange",], aes(x=variable, y=value,fill=year)) +
  geom_violin() + facet_wrap(. ~ variable,scales="free",ncol = 7) +
  theme(legend.position="none",legend.text = element_text(size=14,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 12,face="bold"), axis.text.x=element_blank(),
        axis.title=element_text(size=16,face="bold"),
        strip.text = element_text(size=11,face="bold"))+
  labs(x="",y="Measured value",fill="Season") +scale_fill_grey()

g3<-ggplot(dataf[dataf$trait=="Chlorophyll fluorescence",], aes(x=variable, y=value,fill=year)) +
  geom_violin() + facet_wrap(. ~ variable,scales="free",ncol = 7) +
  theme(legend.position="bottom",legend.text = element_text(size=14,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 12,face="bold"),axis.text.x=element_blank(),
        strip.text = element_text(size=11,face="bold"))+
  labs(x="Traits",y="",fill="Season") +scale_fill_grey()


# plot_grid(NULL,g1, g2, g3, 
#           labels = c("","Anatomical","Gas exchange","Chlorophyll fluorescence"),
#           nrow =4,rel_heights = c(0.07,1,1,2),vjust=-0.1)

plot_grid(g1, g2, g3, 
          labels = c("a","b","c"),
          # labels = c("Anatomical","Gas exchange","Chlorophyll fluorescence"),
          nrow =3,rel_heights = c(1,1,2),axis= "rlbt",hjust=-1 )

ggsave("figures/SF1.png",width = 10, height = 12)
