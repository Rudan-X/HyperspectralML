

rm(list = ls())

list.of.packages <- c("pls","dplyr","plotrix","ggplot2","gridExtra","chemometrics","reshape2","ggpubr")
invisible(lapply(list.of.packages, library, character.only = TRUE))


setwd("C:/Users/Rudan/Documents/GitHub/HyperspectralML/")
for (file in c("R_Burnett/","myR/")){
  filen<-paste0("C:/Users/Rudan/Documents/GitHub/HyperspectralML/R/",file)
  myls <- list.files(path=filen,pattern="*.R")
  myls<-paste0(filen,myls)
  invisible(sapply(myls,FUN=source))
}

vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #

# "SLA","N","CN","Vpmax","Vmax","gsw","a400","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm"

year <- 2021
segments0 <- 5
nrep <- 20
segment <- 3
repl <- 10

agg <- "sampledHSR" #

dataset<-c("Calibration","Validation")
datatype <- c("raw_data","plot_averaged","genotype_averaged")
datalabel <- c("Raw data","Averaged by\n replicate","Averaged by\n genotype") #

notavail <- c()
ini<-TRUE
for (d in 1:length(datatype)){#
  for (inVar in vars){
    filen <- paste0("results/comparePLSR_SVR/PLSR_repeatedCV_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    
    if (file.exists(filen)){
      load(filen)
      res.rcv.random <- res.rcv.mse
      for (t in 1:2){
        if (ini){
          dfr2<-data.frame(score=apply(res.rcv.random$r2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          dfrmse<-data.frame(score=apply(res.rcv.random$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          dfrmseper<-data.frame(score=apply(res.rcv.random$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          dfnoc<-data.frame(NoC=apply(res.rcv.random$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          ini<-FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.random$r2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d]))
          dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.random$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d]))
          dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.random$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d]))
          dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.random$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d]))
        }
      }
    }else{notavail <- c(notavail,filen)}
    
    filen <- paste0("results/unseen_genotypes/PLSR_repeatedCV_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    if (file.exists(filen)){
      load(filen)
      res.rcv.genotype <- res.rcv.mse
      
      
      for (t in 1:2){
        dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.genotype$r2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes\n from seen season",datatype=datalabel[d]))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.genotype$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes\n from seen season",datatype=datalabel[d]))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.genotype$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes\n from seen season",datatype=datalabel[d]))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.genotype$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes\n from seen season",datatype=datalabel[d]))
      }
    }else{notavail <- c(notavail,filen)}
  }
}

datatypes <- c("raw_data","plot_averaged","genotype_averaged")
year_train <- c(2022,2023,2021,2023,2021,2022)
year_test <- c(2021,2021,2022,2022,2023,2023)
for (y in 1:2){
  # split <- paste0("Unseen genotypes\n from unseen ",year_train[y])
  split <- "Unseen genotypes\n from unseen season"
  for (d in 1:3){ 
    datatype <- datatypes[d]
    for (inVar in vars){
      
      filen <- paste0("results/unseen_seasons_and_genotypes/PLSR_model_",datatype,"_",agg,"_training",year_train[y],"_testing",year_test[y],"_i",segment*repl,"_",inVar,".RData")
      if (file.exists(filen)){
        load(filen)
        
        for (t in 1:2){
          
          dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT=split,datatype=datalabel[d]))
          dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT=split,datatype=datalabel[d]))
          dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT=split,datatype=datalabel[d]))
          dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],SPLIT=split,datatype=datalabel[d]))
        }
      }
    }
  }
}


stat_box_data <- function(y) {
  return( 
    data.frame(
      y = max(y),  #may need to modify this depending on your data
      label = paste(round(median(y),2), '\n')
    )
  )
}


dfr2$Trait_type <- varclass[match(dfr2$trait,var1)]
dfr2$Trait_type <- factor(dfr2$Trait_type,levels=c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))
dfr2$trait <- var2[match(dfr2$trait,var1)] # var1 and var2loaded by  myR/change_var_names.Rhttp://127.0.0.1:26459/graphics/cb741dd1-47ef-432f-b3ee-64b48f87b145.png
dfr2$trait <- factor(dfr2$trait, levels=var2[match(vars,var1)])
dfr2$SPLIT <- factor(dfr2$SPLIT,levels=unique(dfr2$SPLIT))
dfr2$datatype <- factor(dfr2$datatype,levels=datalabel)

ggplot(dfr2[dfr2$data_split=="Validation",], aes(x=trait, y=score,fill=SPLIT)) +
  geom_boxplot()+
  facet_grid(datatype ~ Trait_type,scales="free")+ #, switch = "y"
  labs(y="Squared Pearson Correlation",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=SPLIT), size=4,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"))+
  scale_fill_manual(name="Prediction scenario", values =c("#7FBF7F","#C2A5CF","#6A9FB5"))+ #,"#DDC17C"
  stat_compare_means(label = "p.signif",label.y = 0.05) +
  coord_cartesian(ylim =c(0, 1))


ggsave("Figures/SF_sec4.5_unseen_genotypes_and_seasons.png",width = 11, height = 7)


# ggplot(dfr2[dfr2$data_split=="Validation",], aes(x=trait, y=score,fill=Model)) +
#   geom_boxplot()+
#   facet_grid(datatype ~ Trait_type,scales="free", switch = "y")+
#   labs(y="Coefficient of determination",x="") +
#   stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Model), size=3,position=position_dodge(0.9)) +
#   theme(legend.position = "bottom",legend.text = element_text(size=12,face="plain"),legend.title=element_text(size=14,face="bold"),
#         text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
#         axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=16,face="bold"))+
#   scale_fill_brewer(palette = "Dark2") +
#   stat_compare_means(label = "p.signif",label.y.npc = "bottom") 
# 
# 
# ggsave("Figures/MF3_PLSR_SVR_aggregation2.png",width = 13, height = 8)


howl_palette <- c(
  "#496D89", # 沉静蓝（城堡与天空）
  "#C58E62", # 暖棕色（城堡主体与木质元素）
  "#D9AFA8", # 柔和粉（苏菲的裙子与温暖氛围）
  "#9F9D6C", # 草地绿（广阔草原场景）
  "#5A3E36"  # 深褐色（哈尔的头发与复古质感）
)
