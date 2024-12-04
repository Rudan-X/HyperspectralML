library(dplyr)
library(tidyr)
library(ggplot2)
library(caret)
library(reshape2)
library(ggh4x)


rm(list = ls())
setwd("C:/Users/Rudan/Documents/HyperspectralML")
source("script/Rfunctions_data.R")


# trait_ref=list()
# years<-c(2021,2022,2023)
# for (y in 1:3){
#   year<-years[y]
#   accref<-read.csv(paste0("data/accessions/",year,"_maize_reference.csv"))
#   if (y==1){
#     for (i in 1:length(accref$Accession)){
#       temp<-strsplit(accref$Accession[i],'_')
#       accref$Accession[i]<-paste0(temp[[1]][1],'_00',temp[[1]][2])
#     }
#   }
#   trait_ref[[y]]<-accref
# }
# 
# vars<-c("SLA","C",  "N",  "CN", "d13C","d15N","NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
# "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
# "maxNPQ","endNPQ","endFvFm","initialFvFm","Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE" )
# 
# year_comb=matrix(c(2021,2022,2021,2023,2022,2023),3,2,byrow = T)
# year_comb2=matrix(c(1,2,1,3,2,3),3,2,byrow = T)
# cormat<-array(0,c(length(vars),3,20))
# ncount<-array(0,c(length(vars),3,20))
# rownames(cormat)<-vars
# rownames(ncount)<-vars
# 
# for (c in 1:nrow(year_comb)){
#   for (filen in c("SLA","CN_isotopes","fluorescence","ACi")){
#     trait1=load_traits("maize",year_comb[c,1],filen,"plot",FALSE)
#     trait2=load_traits("maize",year_comb[c,2],filen,"plot",FALSE)
#     
#     acc_ref1 <- trait_ref[[year_comb2[c,1]]]
#     ind<-match(trait1$Plot,acc_ref1$Plot)
#     trait1$acc<-acc_ref1$Accession[ind]
#     
#     acc_ref2 <- trait_ref[[year_comb2[c,2]]]
#     ind<-match(trait2$Plot,acc_ref2$Plot)
#     trait2$acc<-acc_ref2$Accession[ind]
#     
#     common_acc<-intersect(trait1$acc,trait2$acc)
#     common_acc<-common_acc[complete.cases(common_acc)]
#     
#     trait1<-trait1[trait1$acc%in%common_acc,]
#     trait2<-trait2[trait2$acc%in%common_acc,]
#     
#     indkeep1<-which(!colnames(trait1)%in%c("acc","Repeat","Plot"))
#     indkeep2<-which(!colnames(trait2)%in%c("acc","Repeat","Plot"))
#     
#     for (r in 1:20){
#       traits_acc1<-matrix(NA,length(common_acc)*2,length(indkeep1))
#       traits_acc2<-matrix(NA,length(common_acc)*2,length(indkeep2))
#       for (a in 1:length(common_acc)){
#         ind<-which(trait1$acc%in%common_acc[a])
#         randind<-sample(seq(1,length(ind)),length(ind),replace=F)
#         traits_acc1[((a-1)*2+1):((a-1)*2+length(randind)),]<-as.matrix(trait1[ind[randind],indkeep1])
#         
#         ind<-which(trait2$acc%in%common_acc[a])
#         randind<-sample(seq(1,length(ind)),length(ind),replace=F)
#         traits_acc2[((a-1)*2+1):((a-1)*2+length(randind)),]<-as.matrix(trait2[ind[randind],indkeep2])
#       }
#       
#       colnames(traits_acc1)<-colnames(trait1[indkeep1])
#       colnames(traits_acc2)<-colnames(trait2[indkeep2])
#       
#       if (ncol(traits_acc1)>1){
#         traits_acc2<-traits_acc2[,match(colnames(traits_acc1),colnames(traits_acc2))]
#       }
#       
#       for (v in 1:ncol(traits_acc1)){
#         indvar<-match(colnames(traits_acc1)[v],vars)
#         cormat[indvar,c,r]<-cor(traits_acc1[,v],traits_acc2[,v],method="spearman",use="complete.obs")
#         ncount[indvar,c,r]<-min(sum(!is.na(traits_acc1[,v])),sum(!is.na(traits_acc1[,v])))
#       }
#     }
#     
#   }
# }
# 
# 
# save(cormat,ncount,file="correlation_between_seasons_spearman_plot_averaged.RData")

# load(file="correlation_traits_spearman.RData")

# cormatspearman=cormat
# load(file="correlation_traits.RData")
load(file="correlation_between_seasons_spearman_plot_averaged.RData")

# save(cormat,ncount,file="correlation_traits_spearman.RData")
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

check<-var2[match(rownames(cormat),var1)]
rownames(cormat)<-check

df<-melt(cormat)
colnames(df)<-c("Variable","Pair","Repetition","Correlation")
df$Pair[df$Pair==1]="2021 vs 2022"
df$Pair[df$Pair==2]="2021 vs 2023"
df$Pair[df$Pair==3]="2022 vs 2023"


finaltraits<-c("SLA","%N","C:N", "Vpmax","Vmax","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2",
               "phiPSII\n recovery\n offset","Maximum\n NPQ","Final\n phiPSII","Fv/Fm")

dffinal<-df[df$Variable%in%finaltraits,]

trait_types<-list()
trait_types[[1]]<-c("SLA","%C","%N","C:N")
trait_types[[2]]<-c("Vmax","Vpmax","A sat light\n 400 CO2","Gsw\n sat light\n 400 CO2")
trait_types[[3]]<-c("phiPSII\n recovery\n offset","Maximum\n NPQ","Final\n phiPSII","Fv/Fm")

dffinal$Trait_type<-"Structural & biochemical"
dffinal$Trait_type[dffinal$Variable%in%trait_types[[2]]]<-"Gas exchange"
dffinal$Trait_type[dffinal$Variable%in%trait_types[[3]]]<-"Chlorophyll fluorescence"
dffinal$Trait_type<-factor(dffinal$Trait_type, c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))





# my_comparisons <- list( c("2021 vs 2022","2021 vs 2023"), c("2021 vs 2022","2022 vs 2023"), c("2021 vs 2023","2022 vs 2023"))
# ggplot(dffinal, aes(x=Pair, y=Correlation)) +
#   geom_boxplot() + 
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y.npc = "top")



ggplot(dffinal, aes(x=Variable, y=Correlation,fill=Pair)) +
  geom_boxplot() +
  scale_fill_manual(values =c("#E7298A", "#66A61E", "#E6AB02"))+
  ggh4x::facet_grid2(~Trait_type, scales = "free_x",space="free_x")+
  theme(legend.position = "bottom",legend.text = element_text(size=18,face="plain"),legend.title=element_text(size=20,face="bold"),
        text = element_text(size = 15,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        strip.text = element_text(size=18,face="bold"),axis.title=element_text(size = 19,face="bold"))+
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Pair), size=4,position=position_dodge(0.9)) +
  coord_cartesian(ylim =c(0, 0.6))

ggsave("figures/SF7.png",width = 13, height = 7)
