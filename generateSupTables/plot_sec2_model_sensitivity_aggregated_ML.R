

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


vars<-c( "SLA","C","N", "CN", "d13C", "d15N",  "Vpmax","Vmax", "SL","a400", "gsw",  "iWUE",
         "NPQ_ind_amp","NPQ_ind_rate", "NPQ_ind_linear", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm","maxNPQ",  "endNPQ")


# vars<-c( ) #

year <- 2021
segments0 <- 5
nrep <- 20
# load genotype-averaged data
agg <- "sampledHSR" #
datatype <- "genotype_averaged"


calc_mean <- function(mat){
  temp <- apply(mat, c(2, 3), mean)
  colnames(temp) <- seq(1,ncol(temp))
  df <- melt(temp)
  df <- df[,2:3]
  colnames(df) <- c("NoC", "score")
  return(df)
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

get_nocR <- function(mat,noc){
  sel_R <- matrix(NA, dim(mat)[1], dim(mat)[2])
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      sel_R[i,j] <- mat[i,j,noc[1,j]]
    }
  }
  return(sel_R)
}


dts <- c("raw_data","plot_averaged","genotype_averaged")
dt_label <- c("Sample","Plot","Genotype")


ini <- T
for (d in 1:length(dts)){
  datatype <- dts[d]
  for (inVar in vars){
    if (datatype=="raw_data"){
      load(file=paste0("results/model_sensitivity/PLSR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
      load(file=paste0("results/MLmethods/model_R2_",datatype,"_",agg,year,"_",inVar,".RData"))
      temp_plsr <- get_nocR(res.singlecv.mse$R2[[2]],matrix(res.singlecv.mse$bestnoc_final,1,nrep))
      temp_plsr <- melt(temp_plsr)
      temp_svr <- melt(best_score[["SVML"]]$R2)
      print(paste0("Trait: ", inVar," SVR best C: ",best_params[["SVML"]]$C))
    }else{
      load(file=paste0("results/model_aggregation/PLSR_model",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
      load(file=paste0("results/model_aggregation/SVR_model",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
      temp_plsr <- get_nocR(res.plsr$R2[[2]],matrix(res.plsr$bestnoc_final,1,nrep))
      temp_plsr <- melt(temp_plsr)
      temp_svr <- melt(res.svr$R2[[2]])
    }
    
    
    df0 <- data.frame(trait=inVar,Rsquared=temp_plsr$value,ML_agg=paste0("PLSR & ",dt_label[d]))
    
    if (ini){
      df <- df0
      ini <- F
    }else{
      df<-rbind(df,df0)
    }
    
    df0 <- data.frame(trait=inVar,Rsquared=temp_svr$value,ML_agg=paste0("SVR & ",dt_label[d]))
    df<-rbind(df,df0)

  }
}


df$trait <- var2_flat[match(df$trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
df$trait <- factor(df$trait, levels=var2_flat[match(vars,var1)])
df$ML_agg <- factor(df$ML_agg,levels = unique(df$ML_agg))

df_sum <- df %>%
  group_by(trait,ML_agg) %>%
  summarise(median_R2 = median(Rsquared),.groups = "drop")

df_sum$trait <- factor(df_sum$trait, levels=var2_flat[match(vars,var1)])
df_sum$ML_agg <- factor(df_sum$ML_agg,levels = unique(df$ML_agg))
                    
df_max <- df_sum %>%
  group_by(trait) %>%
  filter(median_R2 == max(median_R2)) %>%
  ungroup()

df_max$trait <- factor(df_max$trait, levels=var2_flat[match(vars,var1)])
df_max$ML_agg <- factor(df_max$ML_agg,levels = unique(df$ML_agg))

ggplot(data=df, aes(x=ML_agg, y=Rsquared , fill=ML_agg)) +
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
  scale_fill_brewer(name="ML & Aggregation",palette = "Paired") +
  coord_cartesian(ylim = c(-0.4,1)) +
  geom_text(data = df_sum, aes(x = ML_agg, y = 0.8, label = round(median_R2, 2)),
            vjust = -0.5, size = 3, inherit.aes = FALSE) +
  geom_text(data = df_max, aes(x = ML_agg, y = 0.8, label = round(median_R2, 2)),
          vjust = -0.5,  fontface = "bold" ,size = 3, inherit.aes = FALSE)

ggsave("Figures/SF_sec3.1_model_sensitivity_genotype_aggregated.png",width = 8, height = 12)
