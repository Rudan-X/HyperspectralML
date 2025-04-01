

rm(list = ls())

list.of.packages <- c("pls","dplyr","plotrix","ggplot2","gridExtra","chemometrics","reshape2","ggpubr","ggh4x")
invisible(lapply(list.of.packages, library, character.only = TRUE))


setwd("C:/Users/Rudan/Documents/GitHub/HyperspectralML/")
for (file in c("R_Burnett/","myR/")){
  filen<-paste0("C:/Users/Rudan/Documents/GitHub/HyperspectralML/R/",file)
  myls <- list.files(path=filen,pattern="*.R")
  myls<-paste0(filen,myls)
  invisible(sapply(myls,FUN=source))
}

vars <- c("C", "CN", "Vpmax", "Vpmax.Vmax", "Area", "d13C", "d15N", "SL", "gsw",  "iWUE",   "gs", "sl",  
          "NPQ_ind_amp", "NPQ_ind_rate", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res",  "phiPSII_ind_amp",  "phiPSII_ind_rate", "NPQ_ind_linear",  
          "endNPQ", "endFvFm",  "initialFvFm")

#, "A_sat", "gs_sat","ci_sat", "iWUE_sat"
segments0 <- 5
nrep <- 20
segment <- 3
repl <- 10

agg <- "sampledHSR" #

dataset<-c("Calibration","Validation")
datatype <- c("raw_data","plot_averaged","genotype_averaged") #
datalabel <- c("Raw data","Averaged by\n replicate","Averaged by\n genotype") #


ini<-TRUE
for (year in c(2021,2022,2023)){
  
 #,2022,2023
  for (d in 1:length(datatype)){#
    for (inVar in vars){
      filen <- paste0("results/compare_season/PLSR_repeatedCV_model_",datatype[d],"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData") 
      if (file.exists(filen)){
        load(filen)

        for (t in 1:2){
          if (ini){
            dfr2<-data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year)
            dfrmse<-data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year)
            dfrmseper<-data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year)
            dfnoc<-data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year)
            ini<-FALSE
          }else{
            dfr2<-rbind(dfr2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year))
            dfrmse<-rbind(dfrmse,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year))
            dfrmseper<-rbind(dfrmseper,data.frame(score=apply(res.rcv.mse$RMSEPper[[t]], 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year))
            dfnoc<-rbind(dfnoc,data.frame(NoC=apply(res.rcv.mse$bestnoc, 2, mean),trait=inVar,data_split=dataset[t],Model="PLSR",datatype=datalabel[d], Season=year))
          }
          
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
dfr2$trait <- var2[match(dfr2$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfr2$trait <- factor(dfr2$trait, levels=var2[match(vars,var1)])
dfr2$Model <- factor(dfr2$Model,levels=unique(dfr2$Model))
dfr2$datatype <- factor(dfr2$datatype,levels=datalabel)
dfr2$Season <- factor(dfr2$Season,levels=unique(dfr2$Season))

ggplot(dfr2[dfr2$data_split=="Validation",], aes(x=trait, y=score,fill=datatype)) +
  geom_boxplot()+
  facet_grid(Season ~ Trait_type,scales="free",space="free")+ #, switch = "y"
  # ggh4x::facet_grid2(Trait_type ~ datatype, scales = "free", independent = "x")+
  labs(y="Coefficient of determination",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=datatype), size=2.5,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 11),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(name = "Season",palette = "Dark2") +
  #stat_compare_means(label = "p.signif",label.y = -0.4) +
  coord_cartesian(ylim =c(-0.5, 1))


ggsave("Figures/SF_sec3.8_othertraits_PLSR_aggregation_season.png",width = 14, height = 9)


