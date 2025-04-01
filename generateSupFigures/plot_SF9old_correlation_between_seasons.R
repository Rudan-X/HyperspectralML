
rm(list = ls())

list.of.packages <- c("pls","dplyr","plotrix","ggplot2","ggpubr","gridExtra","chemometrics","reshape2")
invisible(lapply(list.of.packages, library, character.only = TRUE))


setwd("C:/Users/Rudan/Documents/GitHub/HyperspectralML/")
for (file in c("R_Burnett/","myR/")){
  filen<-paste0("C:/Users/Rudan/Documents/GitHub/HyperspectralML/R/",file)
  myls <- list.files(path=filen,pattern="*.R")
  myls<-paste0(filen,myls)
  invisible(sapply(myls,FUN=source))
}

vars<-c( "SLA","N","Vmax","a400","phiPSII_ind_res","maxNPQ") #

year1 <- c(2021,2021,2022)
year2 <- c(2022,2023,2023)
yearlabel <- c("2021&2022","2021&2023","2022&2023")



traits <- list()
ini <- TRUE
for (y in 1:3){
  traitHSR1 <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[y],"_plot_averaged.csv"))
  traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR",year2[y],"_plot_averaged.csv"))
  
  ind <- which(colnames(traitHSR1)%in%"mean_350")
  traitHSR1 <- traitHSR1[,seq(1,ind-1)]
  ind <- which(colnames(traitHSR2)%in%"mean_350")
  traitHSR2 <- traitHSR2[,seq(1,ind-1)]
  
  dfy1 <- traitHSR1[order(traitHSR1$Accession, traitHSR1$Plot), ]
  dfy2 <- traitHSR2[order(traitHSR2$Accession, traitHSR2$Plot), ]
  
  rownames(dfy1)[grep("^1",dfy1$Plot)] <- paste0(dfy1$Accession[grep("^1",dfy1$Plot)],"_Plot1")
  rownames(dfy1)[grep("^2",dfy1$Plot)] <- paste0(dfy1$Accession[grep("^2",dfy1$Plot)],"_Plot2")
  
  rownames(dfy2)[grep("^1",dfy2$Plot)] <- paste0(dfy2$Accession[grep("^1",dfy2$Plot)],"_Plot1")
  rownames(dfy2)[grep("^2",dfy2$Plot)] <- paste0(dfy2$Accession[grep("^2",dfy2$Plot)],"_Plot2")
  
  ind <- intersect(rownames(dfy1), rownames(dfy2))
  
  df1 <- dfy1[ind,]
  df2 <- dfy2[ind,]
  
  for (inVar in vars){
    t1 <- df1[,inVar]
    t2 <- df2[,inVar]
    ind <- which(is.na(t1)| is.na(t2))
    if (length(ind)>0){
      t1 <- t1[-ind]
      t2 <- t2[-ind]
    }
    
    if (ini){
      dfcor <- data.frame(year1=t1,year2=t2,Trait=inVar,Season=yearlabel[y],datatype="Averaged by replicate")
      # dfcor <- data.frame(data=c(t1,t2),Plot=c(rep("Plot1",length(t1)),rep("Plot2",length(t2))),Trait=inVar,Season=year)
      ini <- FALSE
    }else{
      dfcor <- rbind(dfcor,data.frame(year1=t1,year2=t2,Trait=inVar,Season=yearlabel[y],datatype="Averaged by replicate"))
      # dfcor <- rbind(dfcor,data.frame(data=c(t1,t2),Plot=c(rep("Plot1",length(t1)),rep("Plot2",length(t2))),Trait=inVar,Season=year))
    }
  }
}


dfcor$Trait <- var2[match(dfcor$Trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfcor$Trait <- factor(dfcor$Trait, levels=var2[match(vars,var1)])
dfcor$Season <- factor(dfcor$Season,levels=unique(dfcor$Season)) 
dfcor$datatype <- factor(dfcor$datatype,levels=c("Raw data","Averaged by replicate" ))

library(ggh4x)

ggplot(dfcor, aes(x=year1, y=year2,color=Season)) +
  geom_point(size=0.5)+
  # ggh4x::facet_grid2(. ~ Trait, scales = "free", independent = "y")+
  facet_wrap(. ~ Trait, scales = "free")+
  labs(y="Year 1",x="Year 2") +
  geom_smooth(method = "lm") + 
  stat_cor(aes(label = ..r.label..))+
  # stat_cor(method="pearson",cor.coef.name = "r",label.y.npc='top',label.x.npc='left') +
  theme(legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13,face="bold"))+
  scale_color_manual(name="Comparison",values =c("#E7298A", "#66A61E", "#E6AB02"))

ggsave("Figures/SF9_correlation_between_seasons.png",width = 9, height = 7)
