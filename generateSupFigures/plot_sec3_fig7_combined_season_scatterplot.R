

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

vars<-c( "SLA","N", "CN", "d13C",  "Vpmax","Vmax", "a400", "gsw",  
         "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate","maxNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")

dtsi<-c( rep(3,3), 2,rep(1,4), rep(1,6),rep(2,3))
MLs <- c(rep("PLSR",3),rep("SVR",5),rep("PLSR",9))

segments0 <- 5
nrep <- 20
repl <- 10
segment <- 3

# load genotype-averaged data
agg <- "sampledHSR"#

dts <- c("raw_data","plot_averaged","genotype_averaged")
dtlabel <- c("Raw data","Plot averaged","Genotype averaged")

ini<-TRUE

for (i in 1:length(vars)){
  inVar <- vars[i]
  datatype <- dts[dtsi[i]]
  traitHSR1 <- read.csv(paste0("data/combined_data/traits_and_HSR2021_",datatype,".csv"))
  traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR2022_",datatype,".csv"))
  traitHSR3 <- read.csv(paste0("data/combined_data/traits_and_HSR2023_",datatype,".csv"))
  
  sametrait <- intersect(colnames(traitHSR1),colnames(traitHSR2))
  sametrait <- intersect(sametrait,colnames(traitHSR3))
  seasons <- c(rep(2021,nrow(traitHSR1)),rep(2022,nrow(traitHSR2)),rep(2023,nrow(traitHSR3)))
  traitHSR <- rbind(traitHSR1[,sametrait],traitHSR2[,sametrait])
  traitHSR <- rbind(traitHSR,traitHSR3[,sametrait])
  
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
  
  seasons <- seasons[!is.na(traitHSR$y)]
  traitHSR <- traitHSR[!is.na(traitHSR$y),]
  if (inVar=="d13C"){
    seasons <- seasons[traitHSR$y>(-30)]
    traitHSR <- traitHSR[traitHSR$y>(-30),]
  }
  if (inVar=="Vmax"){
    seasons <- seasons[traitHSR$y<75]
    traitHSR <- traitHSR[traitHSR$y<75,]
  }
  
  if (inVar=="endNPQ"){
    seasons <- seasons[traitHSR$y<0.5]
    traitHSR <- traitHSR[traitHSR$y<0.5,]
  }
  
  if (inVar=="NPQ_rel_res"){
    seasons <- seasons[traitHSR$y<0.5]
    traitHSR <- traitHSR[traitHSR$y<0.5,]
  }
  
  if (inVar=="NPQ_ind_rate"){
    seasons <- seasons[traitHSR$y<0.04]
    traitHSR <- traitHSR[traitHSR$y<0.04,]
  }
  if (MLs[i]== "PLSR"){
    
    filen=paste0("results/compare_season/PLSR_singleCV_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    
    if (file.exists(filen)){
      load(filen)
      folds <- res.singlecv.mse$folds
      t <- 2
      r <- 1
      pred <- matrix(NA,nrow(res.singlecv.mse$pred),1)
      for (f in 1:length(folds[[r]])){
        val.ind <- folds[[r]][[f]]
        cal.ind <- setdiff(seq(1,nrow(res.singlecv.mse$pred)),val.ind)
        pred[val.ind] <- res.singlecv.mse$pred[val.ind,(r-1)*segments0+f]
      }
      if (exists("res.singlecv.mse")){
        if (ini){
          dfr2<-data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar, Season=seasons)
          ini <- FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar, Season=seasons))
        }
        
      }
    }else{
      filen=paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
      load(filen)
      folds <- res.rcv.mse$folds
      t <- 2
      r <- 1
      pred <- matrix(NA,nrow(res.singlecv.mse$pred),1)
      for (f in 1:length(folds[[r]])){
        val.ind <- folds[[r]][[f]]
        cal.ind <- setdiff(seq(1,nrow(res.rcv.mse$pred)),val.ind)
        pred[val.ind] <- res.rcv.mse$pred[val.ind,(r-1)*segments0+f]
      }
      if (exists("res.singlecv.mse")){
        if (ini){
          dfr2<-data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar, Season=seasons)
          ini <- FALSE
        }else{
          dfr2<-rbind(dfr2,data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar, Season=seasons))
        }
        
      }
    }
  }else if (MLs[i]== "SVR"){
    filen=paste0("results/compare_season/SVR_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    
    if (file.exists(filen)){
      load(filen)
      folds <- res.rcv.svr$folds
      t <- 2
      r <- 1
      pred <- matrix(NA,nrow(res.rcv.svr$pred),1)
      for (f in 1:length(folds[[r]])){
        val.ind <- folds[[r]][[f]]
        cal.ind <- setdiff(seq(1,nrow(res.rcv.svr$pred)),val.ind)
        pred[val.ind] <- res.rcv.svr$pred[val.ind,(r-1)*segments0+f]
      }
      
      dfr2<-rbind(dfr2,data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar, Season=seasons))
    }
  }
}


df <- dfr2
df$trait <- var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=var2_flat[match(vars,var1)])
df$Season <- factor(df$Season,levels = unique(df$Season))


save(df,file= "results/compare_season/Sfig3_7.Rdata")

# load(file= "results/compare_season/Sfig3_7.Rdata")


corr_df <- df %>%
  group_by(trait) %>%
  summarise(
    R2 = 1-sum((Predicted - Measured)^2)/sum((Measured-mean(Measured))^2),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("R² = ", round(R2, 2)),
    x = -Inf,  # left side of the panel
    y = Inf    # top of the panel
  )



ggplot(df, aes(x=Measured, y=Predicted,color=Season)) +
  geom_point(size=0.2)+
  # ggh4x::facet_grid2(trait ~ ., scales = "free", independent = "x")+
  facet_wrap(trait ~ ., ncol=3 ,scales = "free")+
  labs(x="Measured",y="Predicted") +
  geom_smooth(method = "lm", color = "black", se = FALSE) + 
  geom_text(data = corr_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = -0.1, vjust = 1.5,
            size = 3, show.legend = FALSE) +

  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text = element_text(face = "bold")) +
  scale_color_brewer(name = "Season",palette = "Dark2") +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("Figures/SF_sec3.7_combined_seasons_scatterplot.png",width = 8, height = 11)
