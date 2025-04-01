

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
# load genotype-averaged data
aggHSR <- c("fullHSR","sampledHSR") #
aggHSR_lab <- c("Full HSR","Sampled HSR") #
datatype<-c("Calibration","Validation")
segment <- 3

ini<-TRUE
for (npermutations in c(30,90)){
  if (npermutations==30){
    repl <- 10
    label <- "10 x repeated 3-fold CV / 30 permutations"
  }else{#
    repl <- 30
    label <- "30 x repeated 3-fold CV / 90 permutations"
  }
  for (a in 1:2){
    agg <- aggHSR[a]
    for (inVar in vars){
      load(paste0("results/compareMSE_PRESS/repeatedCV_PRESS_vs_MSE_model_o", segments0*nrep,"i",segment*repl,"_",agg,year,"_",inVar,".RData"))
      load(paste0("results/compareMSE_PRESS/permutation_PRESS_vs_MSE_model_o", segments0*nrep,"per",npermutations,"_",agg,year,"_",inVar,".RData"))
      
      for (t in 1:2){
        if (ini){
          dfr2<-data.frame(score=apply(res.rcv.press$R2[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label)
          dfrmse<-data.frame(score=apply(res.rcv.press$RMSEP[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label)
          dfrmseper<-data.frame(score=apply(res.rcv.press$RMSEPper[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label)
          dfnoc<-data.frame(NoC=apply(res.rcv.press$bestnoc, 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label)
          ini<-FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.press$R2[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
          dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.press$RMSEP[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
          dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.press$RMSEPper[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
          dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.press$bestnoc, 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
        }
        
        dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data=datatype[t],Selection="RepeatedCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        
        dfr2<-rbind(dfr2,data.frame(score=apply(res.singlecv.mse$R2[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="SingleCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.singlecv.mse$RMSEP[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="SingleCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.singlecv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="SingleCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.singlecv.mse$bestnoc, 2, mean),trait=inVar,data=datatype[t],Selection="SingleCV_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        
        
        dfr2<-rbind(dfr2,data.frame(score=apply(res.perm.press$R2[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.perm.press$RMSEP[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.perm.press$RMSEPper[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.perm.press$bestnoc, 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_PRESS",HSRtype=aggHSR_lab[a],Npermutation=label))
        
        dfr2<-rbind(dfr2,data.frame(score=apply(res.perm.mse$R2[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.perm.mse$RMSEP[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.perm.mse$RMSEPper[[t]], 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.perm.mse$bestnoc, 2, mean),trait=inVar,data=datatype[t],Selection="Permutation_MSE",HSRtype=aggHSR_lab[a],Npermutation=label))
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


dfr2$trait <- var2[match(dfr2$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfr2$trait <- factor(dfr2$trait, levels=var2[match(vars,var1)])

g <- ggplot(dfr2[dfr2$data=="Calibration",], aes(x=trait, y=score,fill=Selection)) +
  geom_boxplot()+
  facet_grid(HSRtype ~ Npermutation,scales="free", switch = "y")+
  labs(y="Coefficient of determination",x="") +
  # stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Selection), size=2,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=12,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=16,face="bold"))+
  scale_fill_brewer(palette = "Dark2") +
  stat_compare_means(label = "p.signif",label.y.npc = "bottom") 


print(g)

dev.print(file="figures/SF3_calibrationR2.png", device=png, width=2200,height = 1500,res=200)

