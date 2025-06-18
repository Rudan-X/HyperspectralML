
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
         "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm","maxNPQ")

year <- 2021
segments0 <- 5
nrep <- 20
repl <- 10
# load genotype-averaged data
agg <- "sampledHSR"#

dts <- c("raw_data","plot_averaged","genotype_averaged")
dtlabel <- c("Raw data","Plot averaged","Genotype averaged")

segment <- 3

ini<-TRUE

for (d in 1:3){
  datatype <- dts[d]
  for (inVar in vars){
    load(paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    
    t <- 2
    ML <- "PLSR"
    if (ini){
      dfr2<-data.frame(score=melt(res.rcv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(dtlabel[d]," & ", ML))

      ini<-FALSE
    }else{
      dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.mse$R2[[t]])$value,trait=inVar,Aggregation=paste0(dtlabel[d]," & ", ML)))
    }
    
    ML <- "SVR"
    load(file= paste0("results/comparePLSR_SVR/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
    dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.svr$R2[[t]])$value,trait=inVar,Aggregation=paste0(dtlabel[d]," & ", ML)))

  }
  
}

df <- dfr2

df$trait <- var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=var2_flat[match(vars,var1)])
df$Aggregation <- factor(df$Aggregation,levels = unique(df$Aggregation))

df_sum <- df %>%
  group_by(trait,Aggregation) %>%
  summarise(median_R2 = median(score),.groups = "drop")

df_sum$trait <- factor(df_sum$trait, levels=var2_flat[match(vars,var1)])
df_sum$Aggregation <- factor(df_sum$Aggregation,levels = unique(df$Aggregation))

df_max <- df_sum %>%
  group_by(trait) %>%
  filter(median_R2 == max(median_R2)) %>%
  ungroup()

df_max$trait <- factor(df_max$trait, levels=var2_flat[match(vars,var1)])
df_max$Aggregation <- factor(df_max$Aggregation,levels = unique(df$Aggregation))

ggplot(data=df, aes(x=Aggregation, y=score , fill=Aggregation)) +
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
  coord_cartesian(ylim = c(-0.4,1)) +
  geom_text(data = df_sum, 
            aes(x = Aggregation, y = 0.85, label = round(median_R2, 2)),
            vjust = -0.5, size = 3, inherit.aes = FALSE,show.legend = FALSE) +
  geom_text(data = df_max,
            aes(x = Aggregation, y = 0.85, label = round(median_R2, 2)),
            vjust = -0.5,  fontface = "bold" ,size = 3, inherit.aes = FALSE,show.legend = FALSE) +
  scale_fill_brewer(name="Aggregation & ML",palette = "Paired")

ggsave("Figures/SF_sec3.4_PLSR_SVR_aggregation.png",width = 8, height = 10)
