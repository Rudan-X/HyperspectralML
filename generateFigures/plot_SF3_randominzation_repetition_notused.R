library(dplyr)
library(tidyr)
library(ggplot2)
library(pls)
library(doParallel)
library(caret)
library(reshape2)

rm(list = ls())
setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")
###### Result 2 random plot average only within genotype VS randomize across genotype ####

files<-c("SLA","CN_isotopes","fluorescence","ACi")
lists<-list()
lists[[1]]<-c("SLA")
lists[[2]]<-c("C","N","CN","d13C","d15N")
lists[[3]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
              "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
              "maxNPQ","endNPQ","endFvFm","initialFvFm")
lists[[4]]<-c("Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE")

trait_types<-list()
trait_types[[1]]<-c("SLA","C","N","CN","d13C","d15N")
trait_types[[2]]<-c("Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE")
trait_types[[3]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
              "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
              "maxNPQ","endNPQ","endFvFm","initialFvFm")

vars<-c("SLA","C",  "N",  "CN", "d13C","d15N","NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
        "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
        "maxNPQ","endNPQ","endFvFm","initialFvFm","Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE" )
years=c("2021","2022","2023")

nt<-c(1,5,13,7)
cor1<-array(0,c(length(vars),3,20))
rownames(cor1)<-vars
colnames(cor1)<-years

for (y in 1:3){
  year<-years[y]
  for (f in 1:4){
    filen<-files[f]
    traits<-load_traits("maize",year,filen,"plot",FALSE)
    traits<-traits[complete.cases(traits),]
    acc_ref <- read.csv(paste0("data/accessions/",year,"_maize_reference.csv"))
    ind<-match(traits$Plot,acc_ref$Plot)
    traits$acc<-acc_ref$Accession[ind]
    traits=traits[!is.na(traits$acc),]
    uniq_acc<-unique(traits$acc)
    
    for (t in 1:nt[f]){
      for (rep in 1:20){
        traits_acc<-matrix(NA,length(uniq_acc),2)
        for (a in 1:length(uniq_acc)){
          ind<-which(traits$acc%in%uniq_acc[a])
          randind<-sample(c(1,2),length(ind),replace=F)
          traits_acc[a,randind]<-as.matrix(traits[ind,colnames(traits)==lists[[f]][t]])
        }
        traits_acc=traits_acc[complete.cases(traits_acc),]
        cor1[rownames(cor1)%in%lists[[f]][t],y,rep]<-cor(traits_acc[,1],traits_acc[,2])
      }
    }
  }
}


cor2<-array(0,c(length(vars),3,20))
rownames(cor2)<-vars
colnames(cor2)<-years
for (y in 1:3){
  year<-years[y]
  for (f in 1:4){
    filen<-files[f]
    traits<-load_traits("maize",year,filen,"none",FALSE)
    traits<-traits[complete.cases(traits),]
    acc_ref <- read.csv(paste0("data/accessions/",year,"_maize_reference.csv"))
    ind<-match(traits$Plot,acc_ref$Plot)
    traits$acc<-acc_ref$Accession[ind]
    uniq_acc<-unique(traits$acc)


    for (t in 1:nt[f]){
      for (rep in 1:20){
        traits_acc<-matrix(NA,length(uniq_acc)*3,2)
        for (a in 1:length(uniq_acc)){
          ind<-which(traits$acc%in%uniq_acc[a])
          plot_a<-unique(traits$Plot[ind])
          randind<-sample(c(1,2),2,replace=F)
          for (p in 1:2){
            rep_p<-ind[traits$Plot[ind]%in%plot_a[randind[p]]]
            if (length(rep_p)>0){
              traits_acc[((a-1)*3+1):((a-1)*3+length(rep_p)),p]<-as.matrix(traits[rep_p,colnames(traits)==lists[[f]][t]])
            }
          }
        }
        traits_acc=traits_acc[complete.cases(traits_acc),]
        cor2[rownames(cor2)%in%lists[[f]][t],y,rep]<-cor(traits_acc[,1],traits_acc[,2])
      }
    }
  }
}

save(cor1,cor2,file="./results/tables/corr_random_plot.RData")

load("./results/tables/corr_random_plot.RData")
# PLOTTING
df<-melt(cor1)
colnames(df)<-c("Variable","Year","Repetition","Correlation")

df$Year<-factor(df$Year, levels=c("2021","2022","2023"))
df$Variable<-factor(df$Variable, unique(df$Variable))
my_comparisons <- list( c("2021", "2022"), c("2021", "2023"), c("2022", "2023"))

ggplot(df, aes(x=Year, y=Correlation)) +
  geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "top")


df2<-melt(cor2)
colnames(df2)<-c("Variable","Year","Repetition","Correlation")

df2$Year<-factor(df2$Year, levels=c("2021","2022","2023"))
df2$Variable<-factor(df2$Variable, unique(df2$Variable))

my_comparisons <- list( c("2021", "2022"), c("2021", "2023"), c("2022", "2023"))

ggplot(df2, aes(x=Year, y=Correlation)) +
  geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "top")

df$Type<-"Plot-averaged"
df2$Type<-"Raw data"

dftot<-rbind(df,df2)

dftot$Trait_type<-"Structural & biochemical"
dftot$Trait_type[dftot$Variable%in%trait_types[[2]]]<-"Photosynthesis-related"
dftot$Trait_type[dftot$Variable%in%trait_types[[3]]]<-"Chlorophyll fluorescence"

dftot$Trait_type<-factor(dftot$Trait_type, c("Structural & biochemical","Photosynthesis-related","Chlorophyll fluorescence"))

g1<-ggplot(dftot, aes(x=Year, y=Correlation,fill=Type)) +
  geom_boxplot() + 
  stat_compare_means(label = "p.signif",label.y.npc = "top")+
  scale_fill_grey()


finaltraits<-c("SLA","N","CN", "Vpmax","Vmax","a400","gsw",
               "phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

# dffinal<-dftot[dftot$Year=="2021",]
dffinal<-dftot[dftot$Variable%in%finaltraits,]

dffinal$Year<-factor(dffinal$Year,levels=c("2021","2022","2023"))

ggplot(dffinal, aes(x=Variable, y=Correlation,fill=Type)) +
  geom_boxplot() + 
  scale_fill_grey()+ facet_wrap(Year~Trait_type, scales = "free_x")+
  theme(axis.text =  element_text(angle = 30))+
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        strip.text = element_text(size=17,face="bold"),axis.title=element_text(size = 17,face="bold"))
ggsave("figures/SF3_randomization.png",width = 13, height = 13)

