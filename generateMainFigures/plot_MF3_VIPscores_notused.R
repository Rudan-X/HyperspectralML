
rm(list = ls())

list.of.packages <- c("pls","reshape2","caret", "dplyr")
invisible(lapply(list.of.packages, library, character.only = TRUE))

path0 <- "/home/mpimp-golm.mpg.de/xu2004/HyperspectralML/"
if (!dir.exists(path0)) {
  path0 <- "C:/Users/Rudan/Documents/GitHub/HyperspectralML/"
}


setwd(path0)
for (file in c("R_Burnett/","myR/")){
  filen<-paste0(path0,"R/",file)
  myls <- list.files(path=filen,pattern="*.R")
  myls<-paste0(filen,myls)
  invisible(sapply(myls,FUN=source))
}

vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "Vpmax.Vmax", "SL","a400", "gsw",  "iWUE",
         "NPQ_ind_amp","NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm","maxNPQ",  "endNPQ")


vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "Vpmax.Vmax", "SL","a400", "gsw",  "iWUE",
         "NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res", 
         "phiPSII_ind_amp", "phiPSII_ind_rate") #,"NPQ_ind_amp"

year <- 2021

agg <- "sampledHSR"
datatype <- "raw_data"

# MLmethods <- c("Stochastic Gradient Boosting",
#                "eXtreme Gradient Boosting",
#                "SVM linear kernel",
#                "SVM radial kernel",
#                "Neural Network")

MLmethods <- c("XGB",
               "SVM(linear)",
               "SVM(radial)",
               "NN")

MLs <- c("XGB","SVML","SVMR","NN")

ini<-TRUE

for (inVar in vars){
  
  load( file=paste0("results/MLmethods/model_R2_",datatype,"_",agg,year,"_",inVar,".RData"))

  for (i in 1:3){
    m <- MLs[i]
    
    temp1 <- melt(best_score[[m]]$R2)
    temp2 <- melt(best_score[[m]]$RMSEPper)
    
    df0 <- data.frame(Rsquared=temp1$value,RMSE=temp2$value,ML=m,trait=inVar)
    if (ini){
      df <- df0
      ini <- F
    }else{
      df<-rbind(df,df0)
    }
  }
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
reshape_mat <- function(arr){
  reshaped_list <- lapply(1:dim(arr)[3], function(i) {
    slice <- arr[,,i]  # 5 x 20
    matrix(slice, ncol = 1)  # flatten to 100 x 1
  })
  
  # Combine each 100x1 column to get 100 x 40 matrix
  result <- do.call(cbind, reshaped_list)  # 100 x 40
  return(result)
}

segments0 <- 5
nrep <- 20

for (inVar in vars){
  load(file=paste0("results/model_sensitivity/PLSR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
  
  t <- 2
  
  est_noc <- selectNOC_MSE(MSEP = reshape_mat(res.singlecv.mse$RMSEP[[t]]), select_strat ="hastie",sdfact=0.5,repl = nrep, segments = segments0) 
  temp <- get_nocR(res.singlecv.mse$R2[[2]],matrix(est_noc$bestnoc,1,nrep))
  temp <- melt(temp)
  temp2 <- get_nocR(res.singlecv.mse$RMSEP[[2]],matrix(est_noc$bestnoc,1,nrep))
  temp2 <- melt(temp2)
  
  df0 <- data.frame(Rsquared=temp$value,RMSE=temp2$value,ML="PLSR",trait=inVar)
  df<-rbind(df,df0)
}


df$trait <- var2_flat[match(df$trait,var1)] 
df$trait <- factor(df$trait, levels=var2_flat[match(vars,var1)])

df$ML <- factor(df$ML, levels=c(MLs,"PLSR"))

ggplot(data=df, aes(x=ML, y=Rsquared , fill=ML)) +
  geom_boxplot()+
  facet_wrap(.~trait,ncol=3) +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name="ML method",palette = "Dark2") 

df_sum <- df %>% 
  group_by(trait,ML) %>%
  summarise(medianR2 = median(Rsquared),.groups = "drop") 
# +
#   coord_cartesian(ylim=c(0,1))


ggplot(data=df_sum, aes(x=ML, y=medianR2 , fill=ML)) +
  geom_bar(stat="identity") +
  facet_wrap(.~trait,ncol=3) +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name="ML method",palette = "Dark2") +
  geom_text(aes(label=round(medianR2,2)), vjust=-0.3, size=2)


ggsave("Figures/SF_sec2.0_model_MLcomparison_repeatedCV.png",width = 8, height = 12)
