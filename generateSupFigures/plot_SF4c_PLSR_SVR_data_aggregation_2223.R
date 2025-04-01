

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


segments0 <- 5
nrep <- 20
segment <- 3
repl <- 10

agg <- "sampledHSR" #

dataset<-c("Calibration","Validation")
datatype <- c("plot_averaged","genotype_averaged", "raw_data")
datalabel <- c("Averaged per replicate","Averaged per genotype","Raw data") #



for (year in c(2022,2023)){
  ini<-TRUE
  for (inVar in vars){
    if (inVar=="N" | inVar=="CN"){
      endd <- 2
    }else{
      endd <- 3
    }  
    for (d in 1:endd){#
    
      load(paste0("results/comparePLSR_SVR/PLSR_repeatedCV_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      # load(paste0("results/comparePLSR_SVR/tuned_SVR_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
      # res.rcv.svr.tuned <- res.rcv.svr
      load(paste0("results/comparePLSR_SVR/SVR_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
      for (t in 1:2){
        if (ini){
          dfr2<-data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d])
          dfrmse<-data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d])
          dfrmseper<-data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d])
          dfnoc<-data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d])
          ini<-FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d]))
          dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d]))
          dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d]))
          dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d]))
        }
        
        # dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.svr.tuned$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="Tuned SVR",datatype=datalabel[d]))
        # dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.svr.tuned$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="Tuned SVR",datatype=datalabel[d]))
        # dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.svr.tuned$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="Tuned SVR",datatype=datalabel[d]))
        # dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.svr.tuned$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="Tuned SVR",datatype=datalabel[d]))
        # 
        dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.svr$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",datatype=datalabel[d]))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.svr$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",datatype=datalabel[d]))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.svr$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",datatype=datalabel[d]))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.svr$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="SVR",datatype=datalabel[d]))

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
  dfr2$datatype <- factor(dfr2$datatype,levels=c("Raw data","Averaged per replicate","Averaged per genotype"))
  
  
  ggplot(dfr2[dfr2$data_split=="Validation",], aes(x=trait, y=score,fill=datatype)) +
    geom_boxplot()+
    facet_grid(Model ~ Trait_type,scales="free")+ #, switch = "y"
    labs(y="Coefficient of determination",x="") +
    stat_summary(fun.data = stat_box_data,geom = "text", aes(group=datatype), size=4,position=position_dodge(0.9)) +
    theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
          text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 14,face="bold"),
          axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"))+
    scale_fill_manual(name = "Data aggregation",values=c("#DAA520","#7E88C1","#A62E38")) +
    stat_compare_means(label = "p.signif",label.y = -0.4) +
    coord_cartesian(ylim =c(-0.5, 1))
  
  
  ggsave(paste0("Figures/SF4c_PLSR_SVR_aggregation_",year,".png"),width = 10, height = 7)

}


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


