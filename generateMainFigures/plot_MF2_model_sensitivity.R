

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


vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "a400", "gsw",  "iWUE", "SL",
         "NPQ_ind_amp","NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res","maxNPQ",  "endNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")


# vars<-c( ) #

year <- 2021
segments0 <- 5
nrep <- 20
# load genotype-averaged data
aggHSR <- c("fullHSR") #

aggHSR_lab <- c("Full HSR","Sampled HSR") #

dataset<-c("Calibration","Validation")
datatype <- "raw_data"


reshape_mat <- function(arr){
  reshaped_list <- lapply(1:dim(arr)[3], function(i) {
    slice <- arr[,,i]  # 5 x 20
    matrix(slice, ncol = 1)  # flatten to 100 x 1
  })
  
  # Combine each 100x1 column to get 100 x 40 matrix
  result <- do.call(cbind, reshaped_list)  # 100 x 40
  return(result)
}

get_nocR <- function(mat,noc){
  sel_R <- matrix(NA, dim(mat)[1], dim(mat)[2])
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      sel_R[i,j] <- mat[i,j,noc[1,j]]
    }
  }
  return(sel_R)
}


NoC_final <- c()

ini<-TRUE

for (a in 1:1){
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
    
  }
}


df$trait <- long_var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=long_var2_flat[match(vars,var1)])


summary_df <- df %>%
  group_by(trait,NoC, HSRtype) %>%
  summarise(
    Median_R2 =round(median(R2, na.rm = TRUE),2),
    Std_R2 = round(sd(R2, na.rm = TRUE),3),
    lower_R2 = quantile(R2, 0.25, na.rm = TRUE),
    upper_R2 = quantile(R2, 0.75, na.rm = TRUE),
    Median_RMSE =median(RMSE, na.rm = TRUE),
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    lower_RMSE = quantile(RMSE, 0.25, na.rm = TRUE),
    upper_RMSE = quantile(RMSE, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

summary_df <- summary_df %>%
  group_by(trait) %>%
  mutate(
    scale_factor = max(Median_R2, na.rm = TRUE) / max(Median_RMSE, na.rm = TRUE),
    RMSE_scaled       = Median_RMSE * scale_factor,
    RMSE_lower_scaled = lower_RMSE    * scale_factor,
    RMSE_upper_scaled = upper_RMSE    * scale_factor
  ) %>%
  ungroup()


summary_df$NoC <- as.numeric(as.character(summary_df$NoC))

summary_df$trait <- factor(summary_df$trait, levels=long_var2_flat[match(vars,var1)])

bestnoc_df <- data.frame(trait=levels(summary_df$trait),NoC=NoC_final)

bestnoc_summary <- inner_join(bestnoc_df, summary_df, by = c("trait", "NoC"))
bestnoc_summary$trait <- factor(bestnoc_summary$trait, levels = levels(summary_df$trait))


ggplot(summary_df, aes(x = NoC)) +
  
  geom_errorbar(aes(ymin = lower_R2, ymax = upper_R2), width = 0.5, color = "orange") +
  geom_point(aes(y = Median_R2),shape = 1, size = 0.7, color = "orange") +
  
  geom_errorbar(aes(ymin = RMSE_lower_scaled, ymax = RMSE_upper_scaled), width = 0.5, color = "purple") +
  geom_point(aes(y = RMSE_scaled),shape = 1, size = 0.7, color = "purple") +
  
  facet_wrap(~trait,ncol=3) + #, scales = "free_y"
  labs(x = "NoC") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, size=7),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "orange"),
        axis.title.y.right = element_text(color = "purple")) +
  scale_x_continuous(breaks = seq(0, 40, 2)) +
  geom_hline(data = bestnoc_summary, aes(yintercept = Median_R2), color ="black", linetype = "dashed") +
  geom_text(data = bestnoc_summary, aes(x = NoC+1, y = Median_R2 + 0.2, label = round(Median_R2, 2)),color = "black", hjust = 0, size=3) +
  geom_vline(data = bestnoc_summary, aes(xintercept = NoC), linetype = "dashed", color= "blue") +
  geom_text(data = bestnoc_summary, aes(y = -0.3, x = NoC + 1, label = NoC),color = "blue", hjust = 0, size=3) +
  # Left y-axis: R²
  scale_y_continuous(
    name = "R\u00b2",
    
    # Right y-axis: RMSE% (scaled back)
    sec.axis = sec_axis(
      trans = ~ . * max(summary_df$Median_RMSE, na.rm = TRUE) / max(summary_df$Median_R2, na.rm = TRUE),
      name = "RMSE (%)"
    )
  ) 
ggsave("Figures/MF2_model_sensitivity.png",width = 8, height = 12)
