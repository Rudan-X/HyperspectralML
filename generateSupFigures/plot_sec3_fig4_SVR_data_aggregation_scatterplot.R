

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
t <- 1
for (d in 1:length(datatypes)){#
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
    

    load(paste0("results/comparePLSR_SVR/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
    folds <- res.rcv.svr$folds
    for (r in 1:1){
      pred <- matrix(NA,nrow(res.rcv.svr$pred),1)
      for (f in 1:length(folds[[r]])){
        val.ind <- folds[[r]][[f]]
        cal.ind <- setdiff(seq(1,nrow(res.rcv.svr$pred)),val.ind)
        pred[val.ind] <- res.rcv.svr$pred[val.ind,(r-1)*segments0+f]
      }
      if (ini){
        df <- data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar,data_split=dataset[t],Model="SVR",datatype=datalabel[d])
        ini <- FALSE
      }else{
        df <- rbind(df,data.frame(Predicted=pred,Measured=traitHSR$y,trait=inVar,data_split=dataset[t],Model="SVR",datatype=datalabel[d]))
      }
      
    }
  }
}

df$Trait_type <- varclass[match(df$trait,var1)]
df$Trait_type <- factor(df$Trait_type,levels=c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))
df$trait <- var2[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=var2[match(vars,var1)])
df$Model <- factor(df$Model,levels=unique(df$Model))
df$datatype <- factor(df$datatype,levels=datalabel)






pos_df <- df %>%
  group_by(trait) %>%
  summarise(
    x = min(Measured, na.rm = TRUE) + (max(Measured, na.rm = TRUE)- min(Measured, na.rm = TRUE)), 
    miny = min(Predicted, na.rm = TRUE),
    maxy = max(Predicted, na.rm = TRUE),
    # y = min(Measured, na.rm = TRUE) + (max(Measured, na.rm = TRUE) - min(Measured, na.rm = TRUE)) * 0.95,
    .groups = "drop"
  )

rmse_df <- df %>%
  group_by(trait, datatype) %>%
  summarise(
    RMSE = sqrt(mean((Predicted - Measured)^2, na.rm = TRUE)),
    MeanObs = mean(Measured, na.rm = TRUE),
    R2 = cor(Measured, Predicted, use = "complete.obs")^2,
    RMSE_perc = 100 * RMSE / MeanObs,
    .groups = "drop"
  ) %>%
  left_join(pos_df, by = "trait") %>%
  mutate(
    label = paste0("R² = ", round(R2, 2), ", RMSE% = ", round(RMSE_perc, 1)),
    y_nudge = case_when(
              datatype == "Raw data" ~ 0.25,
              datatype == "Averaged by replicate" ~ 0.15,
              datatype == "Averaged by genotype" ~ 0.05),
  ) %>% 
  mutate(
    y= miny + (maxy - miny) * y_nudge
  )



ggplot(df, aes(x=Measured, y=Predicted,color=datatype)) +
  geom_point(size=0.2)+
  # ggh4x::facet_grid2(trait ~ ., scales = "free", independent = "x")+
  facet_wrap(trait ~ ., scales = "free")+
  labs(x="Measured",y="Predicted") +
  geom_smooth(method = "lm") + 
  geom_text(data = rmse_df,
            aes(x = x, y = y, label = label, color=datatype),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5,
            size = 3,show.legend = FALSE)+
  # stat_cor(aes(label = ..r.label..))+
  theme(legend.position = "bottom",legend.text = element_text(size=12,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold")) + 
  scale_color_manual(name = "Data aggregation",values=c("#DAA520","#7E88C1","#A62E38"))

ggsave("Figures/SF_sec3.4_SVR_aggregation_scatterplot.png",width = 10, height = 8)


