rm(list = ls())

list.of.packages <- c("ggplot2","ggpubr","ggh4x")
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


datalabel <- c("Raw data","Averaged by plot","Averaged by genotype") #
datatypes <- c("raw_data", "plot_averaged", "genotype_averaged")

ini <- TRUE
for (d in 1:length(datatypes)){
  datatype <- datatypes[d]
  for (year in c(2021,2022,2023)){
    traitHSR<- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_",datatype,".csv"))
    ind <- which(colnames(traitHSR)%in%"mean_350")
    traits <- traitHSR[,seq(1,ind-1)]
    rm(traitHSR)
    
    for (inVar in vars){
      if (inVar %in% colnames(traits)){
        y <- traits[!is.na(traits[,inVar]),inVar]
        
        if (inVar=="d13C" & datatype=="plot_averaged"){
          y <- y[y>(-30)]
        }else if (inVar=="d13C" & datatype=="genotype_averaged"){
          y <- y[y>(-20)]
        }
        
        if (inVar=="Vmax" & datatype=="raw_data"){
          y <- y[y<75]
        }
        
        if (inVar=="C"){
          y <- y[y<60]
        }
        
        if (inVar=="endNPQ"){
          y <- y[y<0.5]
        }
        
        if (inVar=="NPQ_rel_res"){
          y <- y[y<0.5]
        }
        
        if (inVar=="NPQ_ind_rate"){
          y <- y[y<0.04]
        }
        
        if (ini){
          df <- data.frame(values=y,Trait=inVar,Season=year,datatype=datalabel[d])
          ini <- FALSE
        }else{
          df <- rbind(df,data.frame(values=y,Trait=inVar,Season=year,datatype=datalabel[d]))
        }
      }
    }
  }
}

df$Trait_type <- varclass[match(df$Trait,var1)]
df$Trait_type <- factor(df$Trait_type,levels=c("Structural & biochemical","Gas exchange","Chlorophyll fluorescence"))
df$Trait <- short_var2[match(df$Trait,var1)] # var1 and short_var2loaded by  myR/change_var_names.R
df$Trait <- factor(df$Trait, levels=short_var2[match(vars,var1)])
df$datatype <- factor(df$datatype,levels=datalabel)


ggplot(df[df$Trait_type=="Structural & biochemical",], aes(x=Trait, y=values,fill=datatype)) +
  geom_violin()+
  labs(y="Measured value",x="") +
  ggh4x::facet_grid2(Trait ~ Season,scales="free", independent = "x")+
  scale_fill_manual(name = "Data aggregation",values=c("#F5D491","#7E88C1","#A62E38"))+ 
  theme(legend.position = "bottom",legend.text = element_text(size=11,face="plain"),legend.title=element_text(size=11,face="bold"),
      text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 10,face="bold"),axis.text.x  = element_blank(),
      axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=13,face="bold"))

ggsave("Figures/SF_sec1.1_trait_distribution_structural.png",width = 6, height = 8)


ggplot(df[df$Trait_type=="Gas exchange",], aes(x=Trait, y=values,fill=datatype)) +
  geom_violin()+
  labs(y="Measured value",x="") +
  ggh4x::facet_grid2(Trait ~ Season,scales="free", independent = "x")+
  scale_fill_manual(name = "Data aggregation",values=c("#F5D491","#7E88C1","#A62E38"))+ 
  theme(legend.position = "bottom",legend.text = element_text(size=11,face="plain"),legend.title=element_text(size=11,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 10,face="bold"),axis.text.x  = element_blank(),
        axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=11,face="bold"))

ggsave("Figures/SF_sec1.2_trait_distribution_gasexchange.png",width = 6, height = 8)

ggplot(df[df$Trait_type=="Chlorophyll fluorescence",], aes(x=Trait, y=values,fill=datatype)) +
  geom_violin()+
  labs(y="Measured value",x="") +
  ggh4x::facet_grid2(Trait ~ Season,scales="free", independent = "x")+
  scale_fill_manual(name = "Data aggregation",values=c("#F5D491","#7E88C1","#A62E38"))+ 
  theme(legend.position = "bottom",legend.text = element_text(size=11,face="plain"),legend.title=element_text(size=12,face="bold"),
        text = element_text(size = 12,face="bold"),axis.text =  element_text(size = 10,face="bold"),axis.text.x  = element_blank(),
        axis.title=element_text(size=12,face="bold"), strip.text = element_text(size=11,face="bold"))

ggsave("Figures/SF_sec1.3_trait_distribution_fluorescence.png",width = 6, height = 12)

