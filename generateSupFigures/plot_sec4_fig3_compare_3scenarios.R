

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

plsmethod<-"oscorespls"
maxnoc<-40

agg <- "sampledHSR"


vars<-c( "SLA","N", "CN", "d13C",  "Vpmax","Vmax", "a400", "gsw",  
         "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate","maxNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")

dts <- c("raw_data","plot_averaged","genotype_averaged")
dtsi<-c( rep(3,3), 2,rep(1,4), rep(1,6),rep(2,3))
# MLs <- c(rep("PLSR",3),rep("SVR",5),rep("PLSR",9))
MLs <- c("SVR",rep("PLSR",2),rep("SVR",5),rep("PLSR",9))
dtlabels <- c("Raw data","Plot averaged","Genotype averaged")

segments0 <- 5
nrep <- 20

segment<-3
repl<-1


ini<-TRUE

for (year in c(2021,2022,2023)){
  for (i in 1:length(vars)){
    inVar <- vars[i]
    datatype <- dts[dtsi[i]]
    if (MLs[i]=="PLSR"){
      
      load(file=paste0("results/unseen_genotypes_specific/PLSR_singleCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      res <- res.plsr
    }else if (MLs[i]=="SVR"){
      load(file = paste0("results/unseen_genotypes_specific/SVR_singleCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      res <- res.svr
    }
    if (ini){
      dfr2<-data.frame(score=melt(res$cor2[[2]])$value,trait=inVar,AggregationML=paste0(dtlabels[dtsi[i]], " & ", MLs[i]), Season=year, Scenario = "Unseen genotypes\n from seen season")
      ini<-FALSE
    }else{
      dfr2<-rbind(dfr2,data.frame(score=melt(res$cor2[[2]])$value,trait=inVar,AggregationML=paste0(dtlabels[dtsi[i]], " & ", MLs[i]), Season=year, Scenario = "Unseen genotypes\n from seen season"))
    }
  }
}



for (year_test in c(2021,2022,2023)){
  for (year_train in setdiff(c(2021,2022,2023),year_test)){
    for (i in 1:length(vars)){
      inVar <- vars[i]
      datatype <- dts[dtsi[i]]
      if (MLs[i]=="PLSR"){
        filen <- paste0("results/unseen_seasons_and_genotypes_specific/PLSR_model_scaled_",datatype,"_",agg,"_training",year_train,"_testing",year_test,"_o", segments0*nrep,"_i",segment*repl,"_",inVar,".RData")
        load(filen)
        res <- res.plsr
      }else if (MLs[i]=="SVR"){
        filen <- paste0("results/unseen_seasons_and_genotypes_specific/SVR_model_scaled_",datatype,"_",agg,"_training",year_train,"_testing",year_test,"_o", segments0*nrep,"_i",segment*repl,"_",inVar,".RData")
        load(filen)
        res <- res.svr
      }
      
      
      dfr2<-rbind(dfr2,data.frame(score=melt(res$cor2[[2]])$value,trait=inVar,AggregationML=paste0(dtlabels[dtsi[i]], " & ", MLs[i]), Season=year_test, Scenario = "Unseen genotypes\n from unseen season"))
    }
  }
}

# segment<-3
# repl<-10
# 
# for (year in c(2021,2022,2023)){
#   for (i in 1:length(vars)){
#     inVar <- vars[i]
#     datatype <- dts[dtsi[i]]
#     if (MLs[i]=="PLSR"){
#       
#       if (year==2021){
#         filen <- paste0("results/compareMSE_PRESS/repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
#       }else{
#         filen <- paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
#       }
#       load(filen)
#       res <- res.rcv.mse
#       
#     }else if (MLs[i]=="SVR"){
#       if (year==2021){
#         filen <- paste0("results/comparePLSR_SVR/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData")
#       }else{
#         filen <- paste0("results/compare_season/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData")
#       }
#       load(filen)
#       res <- res.rcv.svr
#     }
#     
#     dfr2<-rbind(dfr2,data.frame(score=melt(res$R2[[2]])$value,trait=inVar,AggregationML=paste0(dtlabels[dtsi[i]], " & ", MLs[i]), Season=year, Scenario = "Random split"))
#   }
# }






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
dfr2$AggregationML <- factor(dfr2$AggregationML, levels=dtlabels)


# Step 1: Compute medians per trait and aggregation
medians <- dfr2 %>%
  group_by(trait, Scenario, Season) %>%
  summarise(med = median(score), .groups = "drop")

# Step 2: Identify best (max) median per trait
best_medians <- medians %>%
  group_by(trait, Season) %>%
  filter(med == max(med)) %>%
  select(trait, Scenario, Season) %>%
  mutate(best = TRUE)

# Step 3: Join back to full data to flag best aggregation
dfr2_labeled <- dfr2 %>%
  left_join(best_medians, by = c("trait", "Scenario", "Season")) %>%
  mutate(best = ifelse(is.na(best), FALSE, best))

label_df <- dfr2_labeled %>%
  group_by(trait, Scenario, Season, best) %>%
  summarise(y = 1,# max(score),
            label = paste0(round(median(score), 2)),
            .groups = "drop")

paired_colors <- c(
  # "Random split" = "#7FBF7F",    
  "Unseen genotypes\n from seen season" = "#C2A5CF",
  "Unseen genotypes\n from unseen season" ="#6A9FB5"

)

paired_colors_label <- c(
  # "Random split" = "#1B7837",       # darker green
  "Unseen genotypes\n from seen season" = "#762A83",    # darker purple
  "Unseen genotypes\n from unseen season" ="#1F78B4"
)

ggplot(dfr2, aes(x = trait, y = score, fill = Scenario)) +
  geom_boxplot() +
  facet_grid(. ~ Season, scales = "free_y") +
  labs(y = "Squared Pearson Correlation", x = "") +
  geom_text(data = label_df,
            aes(x = trait, y = y, label = label, 
                fontface = ifelse(best, "bold", "plain"),
                color = Scenario,
                group = Scenario),
            position = position_dodge(width = 0.9),
            size = 3,show.legend = FALSE,
            inherit.aes = FALSE) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    legend.position = "bottom",
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x =  element_text(angle=45,hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  ) +
  scale_fill_manual(name = "Scenario", values = paired_colors) +
  scale_color_manual(values = paired_colors_label, name = "Scenario")+
  coord_flip(ylim = c(-0.3, 1.1)) +
  stat_compare_means(label = "p.signif",label.y = -0.1) 


ggsave("Figures/SF_sec4.3_compare_scenario2&3scaled.png",width = 8, height = 10)

