

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

get_R2 <- function(meas,pred){
  R2=1-sum((meas-pred)^2)/sum((meas-mean(meas))^2)
  return(R2)
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
datatypes <- c("raw_data","plot_averaged","genotype_averaged")
datalabel <- c("Raw data","Averaged by replicate","Averaged by genotype") #


ini<-TRUE
ini2<-TRUE
for (d in 1:3){#
  datatype <- datatypes[d]
  for (inVar in vars){
    traitHSR <- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_",datatype,".csv"))
    
    if (inVar=="N" | inVar=="CN"){
      Start.wave <- 1500
      End.wave <- 2400
    }else{
      Start.wave <- 400
      End.wave <- 2400
    }
    
    if (agg=="fullHSR"){
      wv <- seq(Start.wave,End.wave,1)
    }else{
      wv <- seq(Start.wave,End.wave,by=5)
    }
    
    spec <- as.matrix(traitHSR[, which(names(traitHSR) %in%  paste0("mean_",wv))])
    traitHSR <- cbind(data.frame(y=traitHSR[,inVar]),spec)
    colnames(traitHSR)[1] <- "y"
    traitHSR <- traitHSR[!is.na(traitHSR$y),]
    
    load(paste0("results/comparePLSR_SVR/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    folds <- res.rcv.mse$folds
    for (r in 1:20){
      pred <- matrix(NA,nrow(res.rcv.mse$pred),1)
      for (f in 1:length(folds[[r]])){
        val.ind <- folds[[r]][[f]]
        cal.ind <- setdiff(seq(1,nrow(res.rcv.mse$pred)),val.ind)
        pred[val.ind] <- res.rcv.mse$pred[val.ind,(r-1)*segments0+f]
      }
      if (ini){
        df_recR2<-data.frame(score=get_R2(traitHSR$y,pred),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="Coefficient of determination")
        df_recRMSEP<-data.frame(score=round(sqrt(mean((traitHSR$y-pred)^2)),2),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="RMSE%")
        
        ini <- FALSE
      }else{
        df_recR2<- rbind(df_recR2,data.frame(score=get_R2(traitHSR$y,pred),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="Coefficient of determination"))
        df_recRMSEP<- rbind(df_recRMSEP,data.frame(score=round(sqrt(mean((traitHSR$y-pred)^2)),2),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="RMSE%"))
      }
    }
        
    t <- 2
    if (ini2){
      df_origR2<-data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="Coefficient of determination")
      df_origRMSEP<-data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="RMSE%")
      ini2<-FALSE
    }else{
      df_origR2<-rbind(df_origR2,data.frame(score=apply(res.rcv.mse$R2[[t]], 2, mean),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="Coefficient of determination"))
      df_origRMSEP<-rbind(df_origRMSEP,data.frame(score=apply(res.rcv.mse$RMSEP[[t]], 2, mean),trait=inVar,Model="PLSR",datatype=datalabel[d],metric="RMSE%"))
    }
  }
}

df_recR2$Gathering <- "Reconstruction"
df_recRMSEP$Gathering <- "Reconstruction"
df_origR2$Gathering <- "Individual folds"
df_origRMSEP$Gathering <- "Individual folds"

df <- rbind(df_origR2,df_origRMSEP,df_recR2,df_recRMSEP)

df$Trait_type <- varclass[match(df$trait,var1)]
df$Trait_type <- factor(df$Trait_type,levels=c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))
df$trait <- var2[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=var2[match(vars,var1)])
df$Model <- factor(df$Model,levels=unique(df$Model))
df$datatype <- factor(df$datatype,levels=datalabel)




stat_box_data <- function(y) {
  return( 
    data.frame(
      y = max(y),  #may need to modify this depending on your data
      label = paste(round(median(y),2), '\n')
    )
  )
}

ggplot(df, aes(x=trait, y=score,fill=Gathering)) +
  geom_boxplot()+
  ggh4x::facet_grid2(datatype ~ metric, scales = "free", independent = "y")+
  labs(y="Evaluation performance",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Gathering), size=4,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"))+
  stat_compare_means(label = "p.signif",label.y = -0.4) +
  scale_fill_brewer(name = "Metric collection") 

ggsave("Figures/SF_sec3.4b_PLSR_aggregation_metrics.png",width = 10, height = 12)


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


