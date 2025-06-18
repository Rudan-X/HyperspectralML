

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


# vars<-c( "SLA","N", "CN", "d13C",  "Vpmax","Vmax", "a400", "gsw",  
#          "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate", 
#          "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm","maxNPQ")

vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "a400", "gsw",  "iWUE","SL",
         "NPQ_ind_amp","NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm","maxNPQ",  "endNPQ")

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
    if (ini){
      dfr2<-data.frame(score=melt(res.rcv.mse$R2[[t]])$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d])
      dfrmse<-data.frame(score=melt(res.rcv.mse$RMSEP[[t]])$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d])
      dfrmseper<-data.frame(score=melt(res.rcv.mse$RMSEPper[[t]])$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d])
      dfnoc<-data.frame(NoC=melt(res.rcv.mse$bestnoc)$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d])


      ini<-FALSE
    }else{
      dfr2<-rbind(dfr2,data.frame(score=melt(res.rcv.mse$R2[[t]])$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d]))
      dfrmse<-rbind(dfrmse,data.frame(score=melt(res.rcv.mse$RMSEP[[t]])$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d]))
      dfrmseper<-rbind(dfrmseper,data.frame(score=melt(res.rcv.mse$RMSEPper[[t]])$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d]))
      dfnoc<-rbind(dfnoc,data.frame(NoC=melt(res.rcv.mse$bestnoc)$value,trait=inVar,Selection="RepeatedCV",Aggregation=dtlabel[d]))
    }

  }

}



# for (a in 1:1){
#   agg <- aggHSR[a]
#   for (inVar in vars){
#     load(paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
#     for (t in 1:2){
#       dfr2<-rbind(dfr2,data.frame(score=melt(res.singlecv.mse$R2[[t]])$value,trait=inVar,Selection="SingleCV",Aggregation=datatype))
#       dfrmse<-rbind(dfrmse,data.frame(score=melt(res.singlecv.mse$RMSEP[[t]])$value,trait=inVar,Selection="SingleCV",Aggregation=datatype))
#       dfrmseper<-rbind(dfrmseper,data.frame(score=melt(res.singlecv.mse$RMSEPper[[t]])$value,trait=inVar,Selection="SingleCV",Aggregation=datatype))
#       dfnoc<-rbind(dfnoc,data.frame(NoC=melt(res.singlecv.mse$bestnoc)$value,trait=inVar,Selection="SingleCV",Aggregation=datatype))
# 
#     }
#   }
# }


dfnoc$scoretype <- "Number of component"
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


dfr2$trait <- var2_flat[match(dfr2$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfr2$trait <- factor(dfr2$trait, levels=rev(var2_flat[match(vars,var1)]))
dfr2$Aggregation <- factor(dfr2$Aggregation, levels=dtlabel)


dfnoc$trait <- var2_flat[match(dfnoc$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfnoc$trait <- factor(dfnoc$trait, levels=rev(var2_flat[match(vars,var1)]))
dfnoc$Aggregation <- factor(dfnoc$Aggregation, levels=dtlabel)


# Step 1: Compute medians per trait and aggregation
medians <- dfr2 %>%
  group_by(trait, Aggregation) %>%
  summarise(med = median(score), .groups = "drop")

# Step 2: Identify best (max) median per trait
best_medians <- medians %>%
  group_by(trait) %>%
  filter(med == max(med)) %>%
  select(trait, Aggregation) %>%
  mutate(best = TRUE)

# Step 3: Join back to full data to flag best aggregation
dfr2_labeled <- dfr2 %>%
  left_join(best_medians, by = c("trait", "Aggregation")) %>%
  mutate(best = ifelse(is.na(best), FALSE, best))

label_df <- dfr2_labeled %>%
  group_by(trait, Aggregation, best) %>%
  summarise(y = 1,# max(score),
            label = paste0(round(median(score), 2)),
            .groups = "drop")

paired_colors <- c(
  "Raw data" = "#A6CEE3",         # light blue
  "Plot averaged" = "#B2DF8A",    # light red/pink
  "Genotype averaged" =  "#FB9A99"# light green
)

paired_colors_label <- c(
  "Raw data"         = "#619CCF",  # darker blue
  "Plot averaged"    = "#7FBF4D",  # deeper pink/red
  "Genotype averaged"=  "#E26C6C"  # medium green
)

g1 <- ggplot(dfr2, aes(x = trait, y = score, fill = Aggregation)) +
  geom_boxplot() +
  facet_grid(. ~ scoretype, scales = "free_y") +
  labs(y = "", x = "") +
  geom_text(data = label_df,
            aes(x = trait, y = y, label = label, 
                fontface = ifelse(best, "bold", "plain"),
                color = Aggregation,
                group = Aggregation),
            position = position_dodge(width = 0.9),
            size = 3,show.legend = FALSE,
            inherit.aes = FALSE) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    legend.position = "none",
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x =  element_text(angle=45,hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  ) +
  scale_fill_manual(name = "Aggregation", values = paired_colors) +
  scale_color_manual(values = paired_colors_label, name = "Aggregation")+
  coord_flip(ylim = c(-1, 1.1))


stat_box_data <- function(y) {
  data.frame(
    y =  max(y),  # vertical position for label
    label = paste(round(median(y), 2), '\n')
  )
}

# Step 1: Compute medians per trait and aggregation
medians <- dfnoc %>%
  group_by(trait, Aggregation) %>%
  summarise(label =  paste0(round(median(score),0)), .groups = "drop")

medians$y <- 35
g2 <- ggplot(dfnoc, aes(x=trait, y=score,fill=Aggregation)) +
  geom_boxplot()+
  facet_grid(. ~ scoretype , scales = "free_y")+
  labs(y="",x="") +
  theme_minimal() + 
  geom_text(data = medians,
            aes(x = trait, y = y, label = label, 
                color = Aggregation,
                group = Aggregation),
            position = position_dodge(width = 0.9),
            size = 3,show.legend = FALSE,
            inherit.aes = FALSE) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "none",
        text = element_text(size = 14,face="bold"),
        axis.text =  element_text(size = 12,face="bold"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=14,face="bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5))+
  scale_fill_manual(name = "Aggregation", values = paired_colors) +
  scale_color_manual(values = paired_colors_label, name = "Aggregation")+
  coord_flip()

ggarrange(g1, g2, ncol = 2, widths = c(1.3,1), legend = "bottom", common.legend = T)

ggsave("Figures/SF_sec3.2_aggregation_repeatedCV.png",width = 10, height = 12)




