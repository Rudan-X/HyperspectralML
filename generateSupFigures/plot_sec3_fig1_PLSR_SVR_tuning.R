

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

vars<-c( "SLA","N","Vmax","a400","phiPSII_ind_res","maxNPQ") #

# "SLA","N","CN","Vpmax","Vmax","gsw","a400","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm"

year <- 2021
segments0 <- 5
nrep <- 20
segment <- 3
repl <- 10

agg <- "sampledHSR" #

dataset<-c("Calibration","Validation")
datatype <- c("raw_data")
datalabel <- c("Raw data") #


ini<-TRUE
for (d in 1:length(datatype)){#
  for (inVar in vars){
    load(paste0("results/comparePLSR_SVR/PLSR_notuned_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"_NOCdefault_",inVar,".RData"))
    res.plsr.notuned <- res.rcv.mse
    load(paste0("results/comparePLSR_SVR/PLSR_repeatedCV_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    
    load(paste0("results/comparePLSR_SVR/tuned_SVR_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
    res.rcv.svr.tuned <- res.rcv.svr
    load(paste0("results/comparePLSR_SVR/tuned_SVR_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
    
    for (t in 1:2){
      if (ini){
        dfr2<-data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d])
        dfrmse<-data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d])
        dfrmseper<-data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d])
        dfnoc<-data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d])
        ini<-FALSE
      }else{
        dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d]))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d]))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d]))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="TRUE",datatype=datalabel[d]))
      }
      
      dfr2<-rbind(dfr2,data.frame(score=apply(res.plsr.notuned$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="FALSE",datatype=datalabel[d]))
      dfrmse<-rbind(dfrmse,data.frame(score=apply(res.plsr.notuned$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="FALSE",datatype=datalabel[d]))
      dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.plsr.notuned$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="FALSE",datatype=datalabel[d]))
      dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.plsr.notuned$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",Tuning="FALSE",datatype=datalabel[d]))
      
      dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.svr.tuned$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="TRUE",datatype=datalabel[d]))
      dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.svr.tuned$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="TRUE",datatype=datalabel[d]))
      dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.svr.tuned$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="TRUE",datatype=datalabel[d]))
      dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.svr.tuned$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="TRUE",datatype=datalabel[d]))
      
      dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.svr$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="FALSE",datatype=datalabel[d]))
      dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.svr$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="FALSE",datatype=datalabel[d]))
      dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.svr$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="FALSE",datatype=datalabel[d]))
      dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.svr$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",Tuning="FALSE",datatype=datalabel[d]))
      
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
dfr2$trait <- var2[match(dfr2$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfr2$trait <- factor(dfr2$trait, levels=var2[match(vars,var1)])
dfr2$Model <- factor(dfr2$Model,levels=unique(dfr2$Model))
dfr2$datatype <- factor(dfr2$datatype,levels=datalabel)

ggplot(dfr2[dfr2$data_split=="Validation",], aes(x=trait, y=score,fill=Tuning)) +
  geom_boxplot()+
  facet_grid(Model ~ Trait_type,scales="free")+
  labs(y="Coefficient of determination",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Tuning), size=3,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"))+
  scale_fill_brewer(name="Hyperparameter tuning",palette = "Dark2") +
  stat_compare_means(label = "p.signif",label.y = -2.4) +
  coord_cartesian(ylim =c(-2.5, 0.85))


ggsave("Figures/SF_sec3.1.png",width = 10, height = 7)
