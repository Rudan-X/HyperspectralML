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

get_R2 <- function(meas,pred){
  R2=1-sum((meas-pred)^2)/sum((meas-mean(meas))^2)
  return(R2)
}

segments0 <- 1
nrep <- 1

repl<-10
segment<-3

agg <- "sampledHSR"

datalabel <- c("Raw data","Averaged per replicate","Averaged per genotype") #
datatypes <- c("raw_data","plot_averaged","genotype_averaged")
vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #"SLA","Vmax","a400","phiPSII_ind_res","maxNPQ"
year1 <- c(2022,2023,2021,2023,2021,2022)
year2 <- c(2021,2021,2022,2022,2023,2023)
ini <- TRUE
for (y in 1:6){
  for (d in 1:3){ 
    datatype <- datatypes[d]
    for (inVar in vars){

      filen <- paste0("results/unseen_seasons/PLSR_model_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_NOC25_",inVar,".RData")
      if (file.exists(filen)){
        load(filen)
        
        if (ini){
          dfcor <- data.frame(score=cor(res.rcv.mse$pred[val.ind],res.rcv.mse$meas[val.ind])^2,Trait=inVar,Training=year1[y],Testing=year2[y],datatype=datalabel[d])
          ini <- FALSE
        }else{
          dfcor <- rbind(dfcor,data.frame(score=cor(res.rcv.mse$pred[val.ind],res.rcv.mse$meas[val.ind])^2,Trait=inVar,Training=year1[y],Testing=year2[y],datatype=datalabel[d]))
        }
      }
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
  for (d in 1:3){ 
    datatype <- datatypes[d] 
    for (inVar in vars){
      filen <- paste0("results/unseen_seasons/PLSR_model_",datatype,"_",agg,"_training",year1[[y]][1],"&",year1[[y]][2],"_testing",year2[y],"_NOC25_",inVar,".RData")
      if (file.exists(filen)){
        load(filen)
        
        dfcor <- rbind(dfcor,data.frame(score=cor(res.rcv.mse$pred[val.ind],res.rcv.mse$meas[val.ind])^2,Trait=inVar,Training="Combined",Testing=year2[y],datatype=datalabel[d]))
      }
    }
  }
}

dfcor$Trait <- var2[match(dfcor$Trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfcor$Trait <- factor(dfcor$Trait, levels=var2[match(vars,var1)])
dfcor$Training <- factor(dfcor$Training,levels=c("2021","2022","2023","Combined"))
dfcor$Testing <- factor(dfcor$Testing,levels=unique(dfcor$Testing)) 
dfcor$datatype <- factor(dfcor$datatype,levels=datalabel)



ggplot(dfcor, aes(x=Trait,y=score,fill=datatype)) +
  geom_bar(stat = "identity", position = "dodge")+
  # theme_minimal() +
  # ggh4x::facet_grid2(Training ~ Testing, scales = "free", independent = "x")+
  facet_grid(Training ~ Testing, scales = "free")+
  labs(y="Squared Pearson Correlation",x="") +
  geom_text(aes(label = round(score,2),y=score+0.01),size = 2.5, vjust=-0.25, position=position_dodge(width=0.9))+
  theme(legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 10,face="bold"),axis.text.x  =  element_text(angle=30),
        axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=13,face="bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50")) +
  scale_fill_manual(name = "Data aggregation",values=c("#74C476","#FFB7C5","#2C7BB6")) +
  coord_cartesian(ylim =c(0, 1)) 

ggsave("Figures/SF10_PLSR_NOC25_predicting_unseen_seasons_barplot_preprocess.png",width = 12, height = 8)
