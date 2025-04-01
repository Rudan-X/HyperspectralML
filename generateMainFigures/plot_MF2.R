

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

dataset<-c("Calibration","Validation")
datatype <- "raw_data"
segment <- 3

ini<-TRUE
for (a in 1:2){
  agg <- aggHSR[a]
  for (inVar in vars){
    repl <- 30
    load(paste0("results/compareMSE_PRESS/permutation_model_",datatype,"_",agg,year,"_o", segments0*nrep,"per",segment*repl,"_",inVar,".RData"))
    
    for (t in 1:2){
      if (ini){
        dfr2<-data.frame(score=apply(res.perm.press$R2[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a])
        dfrmse<-data.frame(score=apply(res.perm.press$RMSEP[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a])
        dfrmseper<-data.frame(score=apply(res.perm.press$RMSEPper[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a])
        dfnoc<-data.frame(NoC=apply(res.perm.press$bestnoc, 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a])
        ini<-FALSE
      }else{
        dfr2<-rbind(dfr2,data.frame(score=apply(res.perm.press$R2[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a]))
        dfrmse<-rbind(dfrmse,data.frame(score=apply(res.perm.press$RMSEP[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a]))
        dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.perm.press$RMSEPper[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a]))
        dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.perm.press$bestnoc, 2, mean),trait=inVar,data=dataset[t],Selection="90 permutations (PRESS)",HSRtype=aggHSR_lab[a]))
      }
      
      repl <- 10
      load(paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="10x repeated 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="10x repeated 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="10x repeated 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data=dataset[t],Selection="10x repeated 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      
      dfr2<-rbind(dfr2,data.frame(score=apply(res.singlecv.mse$R2[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="Single 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      dfrmse<-rbind(dfrmse,data.frame(score=apply(res.singlecv.mse$RMSEP[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="Single 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.singlecv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data=dataset[t],Selection="Single 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
      dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.singlecv.mse$bestnoc, 2, mean),trait=inVar,data=dataset[t],Selection="Single 3-fold CV (MSE)",HSRtype=aggHSR_lab[a]))
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


library(RColorBrewer)
paired_colors <- brewer.pal(9, "Paired")

# Example: assign 3 custom colors from "Dark2"
my_colors <- c(
  "90 permutations (PRESS)" = paired_colors[4],  # dark green
  "10x repeated 3-fold CV (MSE)" = paired_colors[5],  # purple
  "Single 3-fold CV (MSE)" = paired_colors[9]   # mustard
)


ggplot(dfr2[dfr2$data=="Validation",], aes(x=trait, y=score,fill=Selection)) +
  geom_boxplot()+
  facet_grid(HSRtype ~ .,scales="free")+ #, switch = "y"
  labs(y="Coefficient of determination",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Selection), size=4,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=12,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=16,face="bold"))+
  scale_fill_manual(name="NoC selection", values = my_colors)+
  coord_cartesian(ylim =c(0.1, 1)) +
  stat_compare_means(label = "p.signif",label.y = 0.15) 

ggsave("Figures/MF2_validationR2.png",width = 10.5, height = 7)


