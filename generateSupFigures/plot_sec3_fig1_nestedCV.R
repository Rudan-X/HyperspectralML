

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


vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "a400", "gsw",  "iWUE","SL",
         "NPQ_ind_amp","NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res","maxNPQ",  "endNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")

year <- 2021
segments0 <- 5
nrep <- 20
# load genotype-averaged data
aggHSR <- c("sampledHSR") #
aggHSR_lab <- c("Full HSR") #

dataset<-c("Calibration","Validation")
datatype <- "raw_data"
segment <- 3


ini<-TRUE
for (npermutations in c(30)){
  if (npermutations==30){
    repl <- 10
    label <- "10 x repeated 3-fold CV / 30 permutations"
  }else{#
    repl <- 30
    label <- "30 x repeated 3-fold CV / 90 permutations"
  }
  for (a in 1:1){
    agg <- aggHSR[a]
    for (inVar in vars){
      load(paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      # load(paste0("results/compareMSE_PRESS/permutation_model_",datatype,"_",agg,year,"_o", segments0*nrep,"per",segment*repl,"_",inVar,".RData"))

      for (t in 1:2){
        if (ini){
          dfr2<-data.frame(score=melt(res.rcv.mse$R2[[t]])$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label)
          dfrmse<-data.frame(score=melt(res.rcv.mse$RMSEP[[t]])$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label)
          dfrmseper<-data.frame(score=melt(res.rcv.mse$RMSEPper[[t]])$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label)
          dfnoc<-data.frame(NoC=melt(res.rcv.mse$bestnoc)$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label)


          ini<-FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.mse$R2[[t]])$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label))
          dfrmse<-rbind(dfrmse,data.frame(score=melt(res.rcv.mse$RMSEP[[t]])$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label))
          dfrmseper<-rbind(dfrmseper,data.frame(score=melt(res.rcv.mse$RMSEPper[[t]])$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label))
          dfnoc<-rbind(dfnoc,data.frame(NoC=melt(res.rcv.mse$bestnoc)$value,trait=inVar,data=dataset[t],Selection=paste0(repl, "x repeated 3-fold CV"),HSRtype=aggHSR_lab[a],Npermutation=label))
        }
        # dfr2<-rbind(dfr2,data.frame(score=apply(res.perm.press$R2[[t]], 2, mean),trait=inVar,data=dataset[t],Selection=paste0(segment*repl, " permutations(PRESS)"),HSRtype=aggHSR_lab[a],Npermutation=label))
        # dfrmse<-rbind(dfrmse,data.frame(score=apply(res.perm.press$RMSEP[[t]], 2, mean),trait=inVar,data=dataset[t],Selection=paste0(segment*repl, " permutations(PRESS)"),HSRtype=aggHSR_lab[a],Npermutation=label))
        # dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.perm.press$RMSEPper[[t]], 2, mean),trait=inVar,data=dataset[t],Selection=paste0(segment*repl, " permutations(PRESS)"),HSRtype=aggHSR_lab[a],Npermutation=label))
        # dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.perm.press$bestnoc, 2, mean),trait=inVar,data=dataset[t],Selection=paste0(segment*repl, " permutations(PRESS)"),HSRtype=aggHSR_lab[a],Npermutation=label))
      }
    }
  }
}



repl <- 10

for (a in 1:1){
  agg <- aggHSR[a]
  for (inVar in vars){
    load(paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    for (t in 1:2){
      dfr2<-rbind(dfr2,data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,data=dataset[t],Selection="SingleCV",HSRtype=aggHSR_lab[a],Npermutation=label))
      dfrmse<-rbind(dfrmse,data.frame(score=melt(res.singlecv.mse$RMSEP[[t]])$value,trait=inVar,data=dataset[t],Selection="SingleCV",HSRtype=aggHSR_lab[a],Npermutation=label))
      dfrmseper<-rbind(dfrmseper,data.frame(score=melt(res.singlecv.mse$RMSEPper[[t]])$value,trait=inVar,data=dataset[t],Selection="SingleCV",HSRtype=aggHSR_lab[a],Npermutation=label))
      dfnoc<-rbind(dfnoc,data.frame(NoC=melt(res.singlecv.mse$bestnoc)$value,trait=inVar,data=dataset[t],Selection="SingleCV",HSRtype=aggHSR_lab[a],Npermutation=label))

    }
  }
}

sel <- c( # "90 permutations(PRESS)",
         "10x repeated 3-fold CV",
         "SingleCV")


dfnoc$scoretype <- "NoC"
colnames(dfnoc)[1] <- "score"
dfr2$scoretype <- "Coefficient of determination"

df <- rbind(dfnoc,dfr2)

stat_box_data <- function(y) {
  return( 
    data.frame(
      y = max(y),  #may need to modify this depending on your data
      label = paste(round(median(y),2), '\n')
    )
  )
}


df$trait <- var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=rev(var2_flat[match(vars,var1)]))
df$Selection <- factor(df$Selection,levels=sel)
df$scoretype <- factor(df$scoretype,levels=c("Coefficient of determination","NoC"))



ggplot(df[df$data=="Validation",], aes(x=trait, y=score,fill=Selection)) +
  geom_boxplot()+
  facet_grid(. ~ scoretype ,scales="free")+
  labs(y="",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Selection), size=3,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=10,face="plain"),legend.title=element_text(size=12,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=16,face="bold"))+
  scale_fill_brewer(name="Inner-loop",palette = "Paired") +
  stat_compare_means(label = "p.signif",label.y.npc = "top") +
  coord_flip()

ggsave("Figures/SF_sec3.1_nestedCV.png",width = 10, height = 12)




