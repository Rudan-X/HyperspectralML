
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

years <- c(2021,2022,2023)
traits <- list()
ini <- TRUE
for (y in 1:3){
  year <- years[y]
  traitHSR<- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_raw_data.csv"))
  ind <- which(colnames(traitHSR)%in%"mean_350")
  traitHSR <- traitHSR[,seq(1,ind-1)]
  df_shuffled <- traitHSR[order(traitHSR$Accession, traitHSR$Plot), ]
  ind1 <- grep("^1",df_shuffled$Plot)
  ind2 <- grep("^2",df_shuffled$Plot)
  # group <- rep("Plot",nrow(df_shuffled))
  # group[ind1] <- "Plot1"
  # group[ind2] <- "Plot2"
  # 
  df1 <- df_shuffled[ind1, ]
  df2 <- df_shuffled[ind2, ]
  
  temp1 <- paste0(df1$Accession,"_",df1$Repeat)
  temp2 <- paste0(df2$Accession,"_",df2$Repeat)
  df1 <- df1[!duplicated(temp1), ]
  df2 <- df2[!duplicated(temp2), ]
  
  rownames(df1) <- paste0(df1$Accession,"_",df1$Repeat)
  rownames(df2) <- paste0(df2$Accession,"_",df2$Repeat)
  
  ind <- intersect(rownames(df1), rownames(df1))
  
  df1 <- df1[ind,]
  df2 <- df2[ind,]
  
  for (inVar in vars){
    if (any(colnames(traitHSR)%in%inVar)){
      t1 <- df1[,inVar]
      t2 <- df2[,inVar]
      ind <- which(is.na(t1)| is.na(t2))
      if (length(ind)>0){
        t1 <- t1[-ind]
        t2 <- t2[-ind]
      }
      
      if (ini){
        dfcor <- data.frame(plot1=t1,plot2=t2,Trait=inVar,Season=year,datatype="Raw data")
        # dfcor <- data.frame(data=c(t1,t2),Plot=c(rep("Plot1",length(t1)),rep("Plot2",length(t2))),Trait=inVar,Season=year)
        ini <- FALSE
      }else{
        dfcor <- rbind(dfcor,data.frame(plot1=t1,plot2=t2,Trait=inVar,Season=year,datatype="Raw data"))
        # dfcor <- rbind(dfcor,data.frame(data=c(t1,t2),Plot=c(rep("Plot1",length(t1)),rep("Plot2",length(t2))),Trait=inVar,Season=year))
      }
    }
  }
}


for (y in 1:3){
  year <- years[y]
  traitHSR<- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_plot_averaged.csv"))
  ind <- which(colnames(traitHSR)%in%"mean_350")
  traitHSR <- traitHSR[,seq(1,ind-1)]
  df_shuffled <- traitHSR[order(traitHSR$Accession), ]
  ind1 <- grep("^1",df_shuffled$Plot)
  ind2 <- grep("^2",df_shuffled$Plot)
  df1 <- df_shuffled[ind1, ]
  df2 <- df_shuffled[ind2, ]
  
  temp1 <- paste0(df1$Accession)
  temp2 <- paste0(df2$Accession)
  df1 <- df1[!duplicated(temp1), ]
  df2 <- df2[!duplicated(temp2), ]
  rownames(df1) <- temp1
  rownames(df2) <- temp2
  ind <- intersect(rownames(df1), rownames(df1))
  
  df1 <- df1[ind,]
  df2 <- df2[ind,]
  
  for (inVar in vars){
    t1 <- df1[,inVar]
    t2 <- df2[,inVar]
    ind <- which(is.na(t1)| is.na(t2))
    if (length(ind)>0){
      t1 <- t1[-ind]
      t2 <- t2[-ind]
    }
    
    if (ini){
      dfcor <- data.frame(plot1=t1,plot2=t2,Trait=inVar,Season=year,datatype="Averaged by replicate")
      # dfcor <- data.frame(data=c(t1,t2),Plot=c(rep("Plot1",length(t1)),rep("Plot2",length(t2))),Trait=inVar,Season=year)
      ini <- FALSE
    }else{
      dfcor <- rbind(dfcor,data.frame(plot1=t1,plot2=t2,Trait=inVar,Season=year,datatype="Averaged by replicate"))
      # dfcor <- rbind(dfcor,data.frame(data=c(t1,t2),Plot=c(rep("Plot1",length(t1)),rep("Plot2",length(t2))),Trait=inVar,Season=year))
    }
  }
}



dfcor$Trait <- var2[match(dfcor$Trait,var1)] # var1 and var2loaded by  myR/change_var_names.R
dfcor$Trait <- factor(dfcor$Trait, levels=var2[match(vars,var1)])
dfcor$Season <- factor(dfcor$Season,levels=unique(dfcor$Season)) 
dfcor$datatype <- factor(dfcor$datatype,levels=c("Raw data","Averaged by replicate" ))

library(ggh4x)

ggplot(dfcor, aes(x=plot1, y=plot2, color = Season)) +
  geom_point(size=0.5)+
  ggh4x::facet_grid2(Trait ~ datatype, scales = "free", independent = "x")+
  labs(y="Replicate1",x="Replicate2") +
  geom_smooth(method = "lm") + 
  stat_cor(aes(label = ..r.label..))+
  # stat_cor(method="pearson",cor.coef.name = "r",label.y.npc='top',label.x.npc='left') +
  theme(legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=13,face="bold"))+
  scale_color_brewer(name="Season",palette = "Dark2") 

ggsave("Figures/SF_sec3.5_Plot1_vs_Plot2_rawdata_vs_plotaveraged.png",width = 7, height = 12)
