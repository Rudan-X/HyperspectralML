
list.of.packages <- c("pls","reshape2","dplyr","ggplot2")
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


vars<-c(  "SLA","C","N", "CN", "d13C", "d15N", "Vpmax","Vmax", "SL","a400", "gsw",  "iWUE",
         "phiPSII_ind_res","maxNPQ", "NPQ_ind_amp",
         "NPQ_ind_rate", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res",  
         "phiPSII_ind_amp",  "phiPSII_ind_rate", "NPQ_ind_linear",  "endNPQ", "endFvFm",  "initialFvFm") #

year <- 2021

agg <- "fullHSR"
datatype <- "raw_data"

VItypes <- c("PRI","mSR705","EVI","PSRI","ARI","SIPI","SGI","NDVI","REM","NDMI-1","NDII") #

score <- matrix(0,length(vars),length(VItypes)+1)

for (i in 1:length(vars)){
  inVar <- vars[i]
  load( file=paste0("results/VI_models/",datatype,"_",agg,year,inVar,".RData"))
  med <- apply(best_R2_linear,3,median)
  std <- apply(best_R2_linear,3,sd)
  score[i,1:length(VItypes)] <- round(med,2)
}


# for (i in 1:length(vars)){
#   inVar <- vars[i]
#   load( file=paste0("results/VI_models/",datatype,"_",agg,year,inVar,".RData"))
#   med <- apply(best_R2_linear,3,median)
#   std <- apply(best_R2_linear,3,sd)
#   df0 <- data.frame(VI=VItypes,trait=inVar,Rsquared=med,Rsquared=std)
#   if (ini){
#     df <- df0
#     ini <- F
#   }else{
#     df<-rbind(df,df0)
#   }
# }


get_nocR <- function(mat,noc){
  sel_R <- matrix(NA, dim(mat)[1], dim(mat)[2])
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      sel_R[i,j] <- mat[i,j,noc[1,j]]
    }
  }
  return(sel_R)
}      


# for (inVar in vars){
#   load(file=paste0("results/model_sensitivity/PLSR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
#   
#   t <- 2
#   est_noc <- res.singlecv.mse$bestnoc_final
#   temp <- get_nocR(res.singlecv.mse$R2[[2]],matrix(est_noc,1,nrep))
#   
#   med <- mean(best_R2_linear)
#   std <- sd(best_R2_linear)
#   df0 <- data.frame(VI=VItypes,trait=inVar,Rsquared=med,Rsquared=std)
#   df<-rbind(df,df0)
# }

for (i in 1:length(vars)){
  inVar <- vars[i]
  load(file=paste0("results/model_sensitivity/PLSR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
  
  t <- 2
  est_noc <- res.singlecv.mse$bestnoc_final
  temp <- get_nocR(res.singlecv.mse$R2[[2]],matrix(est_noc,1,nrep))
  
  med <- mean(temp)
  std <- sd(temp)
  score[i,ncol(score)] <- round(med,2)
}


colnames(score) <- c(VItypes,"PLSR")
rownames(score) <- vars
ggplot(data=df, aes(x=VI, y=Rsquared , fill=VI)) +
  geom_bar(stat="identity") +
  facet_wrap(.~trait,ncol=3) +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2") +
  geom_text(aes(label=round(Rsquared,2)), vjust=-0.3, size=2)
