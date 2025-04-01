

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

vars<-c("SLA","N","CN","Vpmax","Vmax","gsw","a400","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm") #

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
          dfr2<-data.frame(score=apply(res.rcv.random$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          dfrmse<-data.frame(score=apply(res.rcv.random$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          dfrmseper<-data.frame(score=apply(res.rcv.random$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          dfnoc<-data.frame(NoC=apply(res.rcv.random$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d])
          ini<-FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.random$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Random",datatype=datalabel[d]))
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
        dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.genotype$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes",datatype=datalabel[d]))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.genotype$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes",datatype=datalabel[d]))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.genotype$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes",datatype=datalabel[d]))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.genotype$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],SPLIT="Unseen genotypes",datatype=datalabel[d]))
      }
    }else{notavail <- c(notavail,filen)}
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
dfr2$trait <- var2[match(dfr2$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfr2$trait <- factor(dfr2$trait, levels=var2[match(vars,var1)])
dfr2$SPLIT <- factor(dfr2$SPLIT,levels=unique(dfr2$SPLIT))
dfr2$datatype <- factor(dfr2$datatype,levels=datalabel)

ggplot(dfr2[dfr2$data_split=="Validation",], aes(x=trait, y=score,fill=SPLIT)) +
  geom_boxplot()+
  facet_grid(datatype ~ Trait_type,scales="free")+ #, switch = "y"
  labs(y="Coefficient of determination",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=SPLIT), size=4,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"))+
  scale_fill_manual(name="Prediction scenario", values =c("#7FBF7F","#C2A5CF"))+
  stat_compare_means(label = "p.signif",label.y = -0.4) +
  coord_cartesian(ylim =c(-0.5, 1))


ggsave("Figures/MF4_unseen_genotypes.png",width = 11, height = 7)


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


