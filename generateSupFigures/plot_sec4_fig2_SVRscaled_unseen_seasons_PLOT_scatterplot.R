rm(list = ls())
list.of.packages <- c("pls","dplyr","reshape2","parallel","ggplot2","caret","ggpubr","ggh4x")
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



segments0 <- 1
nrep <- 1

repl<-10
segment<-3

agg <- "sampledHSR"
datatype <-  "plot_averaged"

vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #"SLA","Vmax","a400","phiPSII_ind_res","maxNPQ"
year1 <- c(2022,2023,2021,2023,2021,2022)
year2 <- c(2021,2021,2022,2022,2023,2023)
ini <- TRUE
for (y in 1:6){
  for (inVar in vars){

    # load(file= paste0("results/unseen_seasons/SVR_model_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_",inVar,".RData"))
    load(file= paste0("results/unseen_seasons/SVR_model_scaled_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_",inVar,".RData"))
    
    if (ini){
      dfcor <- data.frame(Predicted=res.rcv.svr$pred[val.ind],Measured=res.rcv.svr$meas[val.ind],Trait=inVar,Training=year1[y],Testing=year2[y],datatype="Averaged by replicate")
      ini <- FALSE
    }else{
      dfcor <- rbind(dfcor,data.frame(Predicted=res.rcv.svr$pred[val.ind],Measured=res.rcv.svr$meas[val.ind],Trait=inVar,Training=year1[y],Testing=year2[y],datatype="Averaged by replicate"))
    }
  }
}
vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") 
year1 <- list()
year1[[1]] <- c(2022,2023)
year1[[2]] <- c(2021,2023)
year1[[3]] <- c(2021,2022)

year2 <- c(2021,2022,2023)
for (y in 1:3){

  for (inVar in vars){
    load(file=paste0("results/unseen_seasons/SVR_model_scaled_",datatype,"_",agg,"_training",year1[[y]][1],"&",year1[[y]][2],"_testing",year2[y],"_",inVar,".RData"))
    dfcor <- rbind(dfcor,data.frame(Predicted=res.rcv.svr$pred[val.ind],Measured=res.rcv.svr$meas[val.ind],Trait=inVar,Training="Combined",Testing=year2[y],datatype="Averaged by replicate"))
  }

}

dfcor$Trait <- var2[match(dfcor$Trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfcor$Trait <- factor(dfcor$Trait, levels=var2[match(vars,var1)])
dfcor$Training <- factor(dfcor$Training,levels=c("2021","2022","2023","Combined"))
dfcor$Testing <- factor(dfcor$Testing,levels=unique(dfcor$Testing)) 
dfcor$datatype <- factor(dfcor$datatype,levels=c("Averaged by replicate" ))



ggplot(dfcor, aes(x=Measured, y=Predicted,color=Training)) +
  geom_point(size=0.02)+
  ggh4x::facet_grid2(Trait ~ Testing, scales = "free", independent = "x")+
  labs(y="Predicted",x="Measured") +
  geom_smooth(method = "lm",se = FALSE, linewidth = 0.5) + 
  stat_cor(aes(label = paste0("italic(r)^2 == ", signif(..r..^2, 2))),label.x.npc = 0.75,label.y.npc = 0.48,size=3)+
  # stat_cor(method="pearson",cor.coef.name = "r",label.y.npc='top',label.x.npc='left') +
  theme(legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 10,face="bold"),
        axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=13,face="bold")) +
  scale_color_manual(name = "Training season",values= c("2021"="#1B9E77","2022"="#D95F02","2023"= "#7570B3","Combined"="#4682B4")) 

ggsave("Figures/SF_sec4.2_SVR_unseen_seasons_PLOT_scatterplot.png",width = 7, height = 11)
