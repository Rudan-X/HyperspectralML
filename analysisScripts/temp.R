
path0 <- "/home/mpimp-golm.mpg.de/xu2004/HyperspectralML/"
if (!dir.exists(path0)) {
  path0 <- "C:/Users/Rudan/Documents/GitHub/HyperspectralML/"
}

setwd(path0)

agg <- "sampledHSR"

vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #
year1 <- c(2022,2023,2021,2023,2021,2022)
year2 <- c(2021,2021,2022,2022,2023,2023)
for (y in 1:6){
  for (datatype in c("plot_averaged")){ 
    for (inVar in vars){
      traitHSR1 <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[y],"_",datatype,".csv"))
      traitHSR1$Season <- year1[y]
      traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR",year2[y],"_",datatype,".csv"))
      traitHSR2$Season <- year2[y]
      sametrait <- intersect(colnames(traitHSR1),colnames(traitHSR2))
      traitHSR <- rbind(traitHSR1[,sametrait],traitHSR2[,sametrait])
      
      if (inVar=="N" | inVar=="CN"){
        Start.wave <- 1500
        End.wave <- 2400
      }else{
        Start.wave <- 400
        End.wave <- 2400
      }
      
      if (agg=="fullHSR"){
        wv <- seq(Start.wave,End.wave,1)
      }else if(agg=="sampledHSR"){
        wv <- seq(Start.wave,End.wave,by=5)
      }
      
      spec <- as.matrix(traitHSR[, which(names(traitHSR) %in%  paste0("mean_",wv))])
      traitHSR <- cbind(data.frame(season=traitHSR$Season,genotype=traitHSR$Accession,y=traitHSR[,inVar]),spec)
      traitHSR <- traitHSR[!is.na(traitHSR$y),]
      
      val.ind <- which(traitHSR$season==year2[y])
      cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
      val.data <- traitHSR[val.ind, ]
      cal.data <- traitHSR[cal.ind, ]
      cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
      cal.data <- data.frame(y=cal.data$y, Spectra=I(cal_spec))
      
      val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
      val.data <- data.frame(y=val.data$y, Spectra=I(val_spec))
      
      
      #####################
      
      
      load(file= paste0("results/unseen_seasons/SVR_model_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_",inVar,".RData"))
      res.rcv.svr$meas <- traitHSR$y
      save(val.ind,res.rcv.svr,file= paste0("results/unseen_seasons/SVR_model_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_",inVar,".RData"))
    }
  }
}
