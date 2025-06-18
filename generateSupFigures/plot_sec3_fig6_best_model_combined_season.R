
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

# "C","d15N", "SL","iWUE","NPQ_ind_rate", "NPQ_ind_linear",, "NPQ_rel_res",  "endNPQ"
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
for (year in c(2021,2022,2023)){
  for (i in 1:length(vars)){
    inVar <- vars[i]
    datatype <- dts[dtsi[i]]
    
    if (MLs[i]== "PLSR"){
      if (year==2021){
        filen1 <- paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
      }else{
        filen1 <- paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
      }
      if (file.exists(filen1)){
        load(filen1)
        t <- 2
        if (ini){
          if (exists("res.singlecv.mse")){
            dfr2<-data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season=year)
            
            ini<-FALSE
            rm(res.singlecv.mse)
          }
        }else{
          if (exists("res.singlecv.mse")){
            dfr2<-rbind(dfr2,data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season=year))
            rm(res.singlecv.mse)
          }
        }
      }
    }else if (MLs[i]== "SVR"){
      if (year==2021){
        filen <- paste0("results/comparePLSR_SVR/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData")
      }else{
        filen <- paste0("results/compare_season/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData")
      }
      
      if (file.exists(filen)){
        load(filen)
        dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.svr$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season=year))
        rm(res.rcv.svr)
        }
    }
  }
}


year1 <- c(2021,2021,2022)
year2 <- c(2022,2023,2023)
yearlabel <- c("2021&2022","2021&2023","2022&2023")


for (y in 1:3){
  for (i in 1:length(vars)){
    inVar <- vars[i]
    datatype <- dts[dtsi[i]]
    if (MLs[i]== "PLSR"){
      filen=paste0("results/compare_season/PLSR_singleCV_model_",datatype,"_",agg,year1[y],"and",year2[y],"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
      
      if (file.exists(filen)){
        load(filen)
        if (exists("res.singlecv.mse")){
          dfr2<-rbind(dfr2,data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season=yearlabel[y]))
          rm(res.singlecv.mse)
        }
      }
    }else if (MLs[i]== "SVR"){
      filen <- paste0("results/compare_season/SVR_repeatedCV_model_",datatype,"_",agg,year1[y],"and",year2[y],"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
      if (file.exists(filen)){
        load(filen)

        dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.svr$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season=yearlabel[y]))
        rm(res.rcv.svr)
      }
    }
  }
}


for (i in 1:length(vars)){
  inVar <- vars[i]
  datatype <- dts[dtsi[i]]
  if (MLs[i]== "PLSR"){
    filen=paste0("results/compare_season/PLSR_singleCV_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    
    
    if (file.exists(filen)){
      load(filen)
      if (exists("res.singlecv.mse")){
        dfr2<-rbind(dfr2,data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season="2021&2022&2023"))
      }
    }else{
      filen=paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
      load(filen)
      if (exists("res.singlecv.mse")){
        dfr2<-rbind(dfr2,data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season="2021&2022&2023"))
      }
    }
  }else if (MLs[i]== "SVR"){
    filen=paste0("results/compare_season/SVR_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    
    if (file.exists(filen)){
      load(filen)
  
      dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.svr$R2[[t]])$value,trait=inVar,Aggregation=paste0(year," & ", MLs[i]), Season="2021&2022&2023"))
    }
  }
}


df <- dfr2

df$trait <- var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=var2_flat[match(vars,var1)])
df$Season <- factor(df$Season,levels = unique(df$Season))


df_sum <- df %>%
  group_by(trait,Season) %>%
  summarise(median_R2 = median(score),.groups = "drop")

df_sum$trait <- factor(df_sum$trait, levels=var2_flat[match(vars,var1)])
# df_sum$Aggregation <- factor(df_sum$Aggregation,levels = unique(df$Aggregation))

df_max <- df_sum %>%
  group_by(trait) %>%
  filter(median_R2 == max(median_R2)) %>%
  ungroup()

df_max$trait <- factor(df_max$trait, levels=var2_flat[match(vars,var1)])
# df_max$Aggregation <- factor(df_max$Aggregation,levels = unique(df$Aggregation))

ggplot(data=df, aes(x=Season, y=score , fill=Season)) +
  geom_violin()+
  facet_wrap(.~trait,ncol=3) +
  labs(x="",y="R\u00b2") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text = element_text(face = "bold")) +
  coord_cartesian(ylim = c(-0.4,1.1)) +
  geom_text(data = df_sum, 
            aes(x = Season, y = 0.95, label = round(median_R2, 2)),
            vjust = -0.5, size = 3, inherit.aes = FALSE,show.legend = FALSE) +
  geom_text(data = df_max,
            aes(x = Season, y = 0.95, label = round(median_R2, 2)),
            vjust = -0.5,  fontface = "bold" ,size = 3, inherit.aes = FALSE,show.legend = FALSE) +
  scale_fill_brewer(name="Season",palette = "Dark2") 

ggsave("Figures/SF_sec3.6_combined_season_singleCV.png",width = 8, height = 10)
