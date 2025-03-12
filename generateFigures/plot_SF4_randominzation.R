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
###### Result 2 random plot average only within genotype VS randomize across genotype ####

files<-c("SLA","CN_isotopes","fluorescence","ACi")
lists<-list()
lists[[1]]<-c("SLA")
lists[[2]]<-c("C","N","CN","d13C","d15N")
lists[[3]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
              "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
              "maxNPQ","endNPQ","endFvFm","initialFvFm")
lists[[4]]<-c("Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE")


vars<-c("SLA","C",  "N",  "CN", "d13C","d15N","NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
        "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
        "maxNPQ","endNPQ","endFvFm","initialFvFm","Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE" )
years=c("2021","2022","2023")

# nt<-c(1,5,13,7)
# cor1<-matrix(0,length(vars),3)
# rownames(cor1)<-vars
# colnames(cor1)<-years
# for (y in 1:3){
#   year<-years[y]
#   for (f in 1:4){
#     filen<-files[f]
#     traits<-load_traits("maize",year,filen,"plot",FALSE)
#     traits<-traits[complete.cases(traits),]
#     acc_ref <- read.csv(paste0("data/accessions/",year,"_maize_reference.csv"))
#     ind<-match(traits$Plot,acc_ref$Plot)
#     traits$acc<-acc_ref$Accession[ind]
#     traits=traits[!is.na(traits$acc),]
#     uniq_acc<-unique(traits$acc)
#     
#     for (t in 1:nt[f]){
#       traits_acc<-matrix(NA,length(uniq_acc),2)
#       for (a in 1:length(uniq_acc)){
#         ind<-which(traits$acc%in%uniq_acc[a])
#         randind<-sample(c(1,2),length(ind),replace=F)
#         traits_acc[a,randind]<-as.matrix(traits[ind,colnames(traits)==lists[[f]][t]])
#       }
#       traits_acc=traits_acc[complete.cases(traits_acc),]
#       cor1[rownames(cor1)%in%lists[[f]][t],y]<-cor(traits_acc[,1],traits_acc[,2])
#     }
#   }
# }
# 
# 
# cor2<-matrix(0,length(vars),3)
# rownames(cor2)<-vars
# colnames(cor2)<-years
# for (y in 1:3){
#   year<-years[y]
#   for (f in 1:4){
#     filen<-files[f]
#     traits<-load_traits("maize",year,filen,"none",FALSE)
#     traits<-traits[complete.cases(traits),]
#     acc_ref <- read.csv(paste0("data/accessions/",year,"_maize_reference.csv"))
#     ind<-match(traits$Plot,acc_ref$Plot)
#     traits$acc<-acc_ref$Accession[ind]
#     uniq_acc<-unique(traits$acc)
#     
#     
#     for (t in 1:nt[f]){
#       traits_acc<-matrix(NA,length(uniq_acc)*3,2)
#       for (a in 1:length(uniq_acc)){
#         ind<-which(traits$acc%in%uniq_acc[a])
#         plot_a<-unique(traits$Plot[ind])
#         randind<-sample(c(1,2),2,replace=F)
#         for (p in 1:2){
#           rep_p<-ind[traits$Plot[ind]%in%plot_a[randind[p]]]
#           if (length(rep_p)>0){
#             traits_acc[((a-1)*3+1):((a-1)*3+length(rep_p)),p]<-as.matrix(traits[rep_p,colnames(traits)==lists[[f]][t]])
#           }
#         }
#       }
#       traits_acc=traits_acc[complete.cases(traits_acc),]
#       cor2[rownames(cor2)%in%lists[[f]][t],y]<-cor(traits_acc[,1],traits_acc[,2])
#     }
#   }
# }
# 
# save(cor1,cor2,file="corr_random_plot.RData")

load("corr_random_plot.RData")

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

rownames(cor1)<-var2[match(rownames(cor1),var1)]
rownames(cor2)<-var2[match(rownames(cor2),var1)]
# PLOTTING
df<-melt(cor1)
colnames(df)<-c("Variable","Year","Correlation")

df$Year<-factor(df$Year, levels=c("2021","2022","2023"))
df$Variable<-factor(df$Variable, unique(df$Variable))
my_comparisons <- list( c("2021", "2022"), c("2021", "2023"), c("2022", "2023"))

ggplot(df, aes(x=Year, y=Correlation)) +
  geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "top")


df2<-melt(cor2)
colnames(df2)<-c("Variable","Year","Correlation")

df2$Year<-factor(df2$Year, levels=c("2021","2022","2023"))
df2$Variable<-factor(df2$Variable, unique(df2$Variable))

my_comparisons <- list( c("2021", "2022"), c("2021", "2023"), c("2022", "2023"))

ggplot(df2, aes(x=Year, y=Correlation)) +
  geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "top")

df$Type<-"Plot-averaged"
df2$Type<-"Raw data"

dftot<-rbind(df,df2)

trait_types<-list()
trait_types[[1]]<-c("SLA","%C","%N","C:N")
trait_types[[2]]<-c("Vmax","Vpmax","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2")
trait_types[[3]]<-c("phiPSII\n induction\n residual","Maximum\n NPQ","Final\n phiPSII","Fv/Fm")


dftot$Trait_type<-"Structural & biochemical"
dftot$Trait_type[dftot$Variable%in%trait_types[[2]]]<-"Gas exchange"
dftot$Trait_type[dftot$Variable%in%trait_types[[3]]]<-"Chlorophyll fluorescence"

dftot$Trait_type<-factor(dftot$Trait_type, c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))

g1<-ggplot(dftot, aes(x=Year, y=Correlation,fill=Type)) +
  geom_boxplot() + 
  stat_compare_means(label = "p.signif",label.y.npc = "top")+
  scale_fill_grey()


finaltraits<-c("SLA","%N","C:N", "Vpmax","Vmax","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2",
               "phiPSII\n induction\n residual","Maximum\n NPQ","Final\n phiPSII","Fv/Fm")

dffinal<-dftot[dftot$Year=="2021",]
dffinal<-dffinal[dffinal$Variable%in%finaltraits,]


ggplot(dffinal, aes(x=Variable, y=Correlation,fill=Type)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_grey()+ facet_wrap(.~Trait_type, scales = "free_x")+
  theme(axis.text =  element_text(angle = 30))+
  geom_text(aes(label = round(Correlation,2),y=Correlation+0.01),size = 5, vjust=-0.25, position=position_dodge(width=0.9))+
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        strip.text = element_text(size=17,face="bold"),axis.title=element_text(size = 17,face="bold"))
ggsave("figures/SF4_randomization.png",width = 13, height = 8)

