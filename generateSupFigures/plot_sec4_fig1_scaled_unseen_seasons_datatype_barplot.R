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

vars<-c( "SLA","N", "CN", "d13C",  "Vpmax","Vmax", "a400", "gsw",  
         "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate","maxNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")

dts <- c("raw_data","plot_averaged","genotype_averaged")
dtsi<-c( rep(3,3), 2,rep(1,4), rep(1,6),rep(2,3))
# MLs <- c(rep("PLSR",3),rep("SVR",5),rep("PLSR",9))
MLs <- c("SVR",rep("PLSR",2),rep("SVR",5),rep("PLSR",9))

dtlabels <- c("Raw data","Plot averaged","Genotype averaged")
agg <- "sampledHSR"

segments0 <- 1
nrep <- 1

repl<-10
segment<-3

year1 <- c(2022,2023,2021,2023,2021,2022)
year2 <- c(2021,2021,2022,2022,2023,2023)
ini <- TRUE
for (y in 1:6){
  for (i in 1:length(vars)){
    inVar <- vars[i]
    datatype <- dts[dtsi[i]]
    if (MLs[i]=="PLSR"){
      filename <- paste0("results/unseen_seasons_specific/PLSR_model_scaled_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_i",segment*repl,"_",inVar,".RData")
    }else if (MLs[i]=="SVR"){
      filename <- paste0("results/unseen_seasons_specific/SVR_model_scaled_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_i",segment*repl,"_",inVar,".RData")
    }
    load(filename)
    val.ind <- res$val.ind
    if (ini){
      dfcor <- data.frame(score=cor(res$pred[val.ind],res$meas[val.ind])^2,Trait=inVar,Training=paste0("Train ", year1[y]),Testing=paste0("Test ", year2[y]),Aggregation=paste0(dtlabels[dtsi[i]]," & ", MLs[i]))
      ini <- FALSE
    }else{
      dfcor <- rbind(dfcor,data.frame(score=cor(res$pred[val.ind],res$meas[val.ind])^2,Trait=inVar,Training=paste0("Train ", year1[y]),Testing=paste0("Test ", year2[y]),Aggregation=paste0(dtlabels[dtsi[i]]," & ", MLs[i])))
    }
  }
}

year1 <- list()
year1[[1]] <- c(2022,2023)
year1[[2]] <- c(2021,2023)
year1[[3]] <- c(2021,2022)

year2 <- c(2021,2022,2023)

for (y in 1:3){
  for (i in 1:length(vars)){
    inVar <- vars[i]
    datatype <- dts[dtsi[i]]
    if (MLs[i]=="PLSR"){
      filename <- paste0("results/unseen_seasons_specific/PLSR_model_scaled_",datatype,"_",agg,"_training",year1[[y]][1],"&",year1[[y]][2],"_testing",year2[y],"_i",segment*repl,"_",inVar,".RData")
    }else if (MLs[i]=="SVR") {
      filename <- paste0("results/unseen_seasons_specific/SVR_model_scaled_",datatype,"_",agg,"_training",year1[[y]][1],"&",year1[[y]][2],"_testing",year2[y],"_i",segment*repl,"_",inVar,".RData")
    }

    load(filename)
    val.ind <- res$val.ind
    dfcor <- rbind(dfcor,data.frame(score=cor(res$pred[val.ind],res$meas[val.ind])^2,Trait=inVar,Training="Combined seasons",Testing=paste0("Test ", year2[y]),Aggregation=paste0(dtlabels[dtsi[i]]," & ", MLs[i])))
  }
}

dfcor$Trait <- var2_flat[match(dfcor$Trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfcor$Trait <- factor(dfcor$Trait, levels=rev(var2_flat[match(vars,var1)]))
dfcor$Training <- factor(dfcor$Training,levels=c("Train 2021","Train 2022","Train 2023","Combined seasons"))
dfcor$Testing <- factor(dfcor$Testing,levels=unique(dfcor$Testing)) 

aggml_comb <- c("Raw data & PLSR","Raw data & SVR", 
                "Plot averaged & PLSR","Plot averaged & SVR",
                "Genotype averaged & PLSR","Genotype averaged & SVR" )

dfcor$Aggregation <- factor(dfcor$Aggregation,levels = aggml_comb)

paired_colors <- c(
  "Raw data & PLSR" = "#A6CEE3",
  "Raw data & SVR" = "#1F78B4", 
  "Plot averaged & PLSR" = "#B2DF8A",
  "Plot averaged & SVR" = "#33A02C",
  "Genotype averaged & PLSR" = "#FB9A99",
  "Genotype averaged & SVR" = "#E31A1C")


ggplot(dfcor, aes(x=Trait,y=score,fill=Aggregation)) +
  geom_bar(stat = "identity", position = "dodge")+
  # theme_minimal() +
  # ggh4x::facet_grid2(Training ~ Testing, scales = "free", independent = "x")+
  facet_grid(Training ~ Testing, scales = "free")+
  labs(y="Squared Pearson Correlation",x="") +
  geom_text(aes(label = round(score,2),y=score+0.01),size = 2.5, hjust=0.1, position=position_dodge(width=0.9))+
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    legend.position = "bottom",
        # legend.text = element_text(size=9,face="plain"),
        # legend.title=element_text(size=10,face="bold"),
        # text = element_text(size = 10,face="bold"),
        # axis.text =  element_text(size = 10),
        # axis.text.x  =  element_text(angle=30,hjust = 1),
        # axis.title=element_text(size=11,face="bold"), 
    strip.text = element_text(size=10,face="bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) +
  scale_fill_manual(name="Aggregation & ML",values=paired_colors) + #"#74C476","#FFB7C5","#2C7BB6"
  coord_flip(ylim =c(0, 1)) + guides(fill = guide_legend(nrow = 2))

ggsave("Figures/SF_sec4.1_scaled_predicting_unseen_seasons_barplot.png",width = 8, height = 10)
