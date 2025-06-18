

rm(list = ls())

list.of.packages <- c("dplyr","reshape2","openxlsx")
invisible(lapply(list.of.packages, library, character.only = TRUE))


setwd("C:/Users/Rudan/Documents/GitHub/HyperspectralML/")
for (file in c("R_Burnett/","myR/")){
  filen<-paste0("C:/Users/Rudan/Documents/GitHub/HyperspectralML/R/",file)
  myls <- list.files(path=filen,pattern="*.R")
  myls<-paste0(filen,myls)
  invisible(sapply(myls,FUN=source))
}


vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "a400", "gsw",  "iWUE","SL",
         "NPQ_ind_amp","NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm","maxNPQ",  "endNPQ")

year <- 2021
segments0 <- 5
nrep <- 20
# load genotype-averaged data
aggHSR <- c("fullHSR","sampledHSR") #
aggHSR_lab <- c("Full HSR","Sampled HSR") #

dataset<-c("Calibration","Validation")
datatype <- "raw_data"






tables <- list()

for (a in 1:2){
  best_noc <- c()
  ini<-TRUE
  agg <- aggHSR[a]
  for (inVar in vars){

    load(file=paste0("results/model_sensitivity/PLSR_model_permutation_",datatype,"_",agg,year,"_",inVar,".RData"))
    t <- 2
    temp <- melt(res.press$R2[[t]])
    
    temp2 <- melt(res.press$RMSEPper[[t]])
    
    
    df0 <- cbind(temp[,2:3],temp2[,3])
    colnames(df0) <- c("NoC","R2","RMSEper")
    df0$trait <- inVar
    df0$HSRtype <- aggHSR_lab[a]
    
    
    if (ini){
      df <- df0
      ini <- F
    }else{
      df<-rbind(df,df0)
    }
    if ("bestno" %in% names(res.press)){
      best_noc <- c(best_noc, res.press$bestno)
    }else{
      best_noc <- c(best_noc, res.press$bestnoc)
    }
    rm(res.press)
  }
  
  df$trait <- long_var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
  df$trait <- factor(df$trait, levels=long_var2_flat[match(vars,var1)])
  
  
  summary_df <- df %>%
    group_by(trait,NoC, HSRtype) %>%
    summarise(
      Median_R2 =round(median(R2, na.rm = TRUE),2),
      Std_R2 = round(sd(R2, na.rm = TRUE),3),
      lowerQR_R2 = quantile(R2, 0.25, na.rm = TRUE),
      upperQR_R2 = quantile(R2, 0.75, na.rm = TRUE),
      Median_RMSEper =median(RMSEper, na.rm = TRUE),
      Mean_RMSEper = mean(RMSEper, na.rm = TRUE),
      lowerQR_RMSEper = quantile(RMSEper, 0.25, na.rm = TRUE),
      upperQR_RMSEper = quantile(RMSEper, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  
  summary_df$NoC <- as.numeric(as.character(summary_df$NoC))
  
  summary_df$trait <- factor(summary_df$trait, levels=long_var2_flat[match(vars,var1)])
  
  bestnoc_df <- data.frame(trait=levels(summary_df$trait),NoC=best_noc)
  
  bestnoc_summary <- inner_join(bestnoc_df, summary_df, by = c("trait", "NoC"))
  bestnoc_summary$trait <- factor(bestnoc_summary$trait, levels = levels(summary_df$trait))
  
  tables[[a]] <- bestnoc_summary[,-c(3,9)]
  colnames(tables[[a]]) <- gsub("RMSEper","RMSE%",colnames(tables[[a]]))
}


write.xlsx(
  x = list("Full HSR" = tables[[1]], "Sub-sampled HSR" = tables[[2]]),
  file = "results/supp_tables/Calibration_PRESS_based.xlsx"
)
write.csv(table,file = "results/supp_tables/PRESS_selection_statistics.csv", row.names = F)
