

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


vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "a400", "gsw",  "iWUE", "SL",
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
for (a in 2:2){
  NoC_final <- c()
  NoC_opt_sd <- c()
  R2_opt_sd <- c()
  R2_opt_med <- c()
  RMSE_opt_sd <- c()
  ini<-TRUE
  
  agg <- aggHSR[a]
  for (inVar in vars){
    load(file=paste0("results/model_sensitivity/PLSR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
    
    t <- 2
    temp <- reshape_mat(res.singlecv.mse$R2[[t]])
    colnames(temp) <- seq(1,ncol(temp))
    temp <- melt(temp)
    
    temp2 <- reshape_mat(res.singlecv.mse$RMSEPper[[t]])
    colnames(temp2) <- seq(1,ncol(temp2))
    temp2 <- melt(temp2)
    
    df0 <- cbind(temp[,2:3],temp2[,3])
    colnames(df0) <- c("NoC","R2","RMSE")
    df0$trait <- inVar
    df0$HSRtype <- aggHSR_lab[a]
    
    
    if (ini){
      df <- df0
      ini <- F
    }else{
      df<-rbind(df,df0)
    }
    
    NoC_final <- c(NoC_final, res.singlecv.mse$bestnoc_final)
    
    NoC_opt_sd <- c(NoC_opt_sd, sd(res.singlecv.mse$bestnoc_MSEP, na.rm = T))
    
    R2_opt_sd <- c(R2_opt_sd,sd(get_nocR(res.singlecv.mse$R2[[2]],res.singlecv.mse$bestnoc_MSEP)))
    R2_opt_med <- c(R2_opt_med, median(get_nocR(res.singlecv.mse$R2[[2]],res.singlecv.mse$bestnoc_MSEP)))
    
    RMSE_opt_sd <- c(RMSE_opt_sd,sd(get_nocR(res.singlecv.mse$RMSEPper[[2]],res.singlecv.mse$bestnoc_MSEP)))
    
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
      Median_RMSE =median(RMSE, na.rm = TRUE),
      Mean_RMSE = mean(RMSE, na.rm = TRUE),
      lowerQR_RMSE = quantile(RMSE, 0.25, na.rm = TRUE),
      upperQR_RMSE = quantile(RMSE, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  summary_df$NoC <- as.numeric(as.character(summary_df$NoC))
  
  summary_df$trait <- factor(summary_df$trait, levels=long_var2_flat[match(vars,var1)])
  
  bestnoc_df <- data.frame(trait=levels(summary_df$trait),NoC=NoC_final,NoC_opt_sd=round(NoC_opt_sd,0),R2_opt_sd=round(R2_opt_sd,3),R2_opt_med=round(R2_opt_med,2),RMSE_opt_sd=round(RMSE_opt_sd,1))
  
  bestnoc_summary <- inner_join(bestnoc_df, summary_df, by = c("trait", "NoC"))
  bestnoc_summary$trait <- factor(bestnoc_summary$trait, levels = levels(summary_df$trait))
  
  
  tables[[a]] <- bestnoc_summary[,-c(7,13)]
  colnames(tables[[a]]) <- gsub("RMSE","RMSE%",colnames(tables[[a]]))
}

write.xlsx(
  x = list("Full HSR" = tables[[1]], "Sub-sampled HSR" = tables[[2]]),
  file = "results/supp_tables/Calibration_MSE_based.xlsx"
)
write.csv(table,file = "results/supp_tables/MSE_selection_statistics.csv", row.names = F)

