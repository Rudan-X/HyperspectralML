
list.of.packages <- c("pls","dplyr","reshape2","parallel","ggplot2")
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


year <- 2021

plsmethod<-"oscorespls"
maxnoc<-25

segments0 <- 5
nrep <- 20

repl<-10
segment<-3

agg <- "sampledHSR"
datatypes <- c("raw_data")

vars<-c("CN") #"N",
for (inVar in vars){
  for (d in 1:1){ # 
    datatype <- datatypes[d]
    traitHSR <- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_",datatype,".csv"))
    
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
    traitHSR <- cbind(data.frame(genotype=traitHSR$Accession,y=traitHSR[,inVar]),spec)
    
    traitHSR <- traitHSR[!is.na(traitHSR$y),]
    
    templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
    temp <- matrix(0,segments0,nrep)
    temp2 <- matrix(0,nrow(traitHSR),segments0*nrep)

    filen <- paste0("results/unseen_genotypes/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData")
    if (file.exists(filen)){
      load(filen)
      res.rcv.mse$r2 <- templ
      res.singlecv.mse$r2 <- templ
      folds <- res.rcv.mse$folds
      for (r in 1:nrep){
        print(paste0("Current repetition: ",r))
        for (f in 1:length(folds[[r]])){
          val.ind <- folds[[r]][[f]]
          cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
          val.data <- traitHSR[val.ind, ]
          cal.data <- traitHSR[cal.ind, ]
          cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
          cal.data <- data.frame(genotype=cal.data$genotype,y=cal.data$y, Spectra=I(cal_spec))
          
          val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
          val.data <- data.frame(genotype=val.data$genotype,y=val.data$y, Spectra=I(val_spec))
          
          datasets<-list(cal.data,val.data)
          
          #####################
          
          for (t in 1:2){
            
            if (t==1){
              pred.mse<-res.rcv.mse$pred[cal.ind,(r-1)*segments0+f]
              pred.smse<-res.singlecv.mse$pred[cal.ind,(r-1)*segments0+f]
            }else{
              pred.mse<-res.rcv.mse$pred[val.ind,(r-1)*segments0+f]
              pred.smse<-res.singlecv.mse$pred[val.ind,(r-1)*segments0+f]
            }
            
            res.rcv.mse$r2[[t]][f,r] <- round(cor(pred.mse, datasets[[t]]$y)^2,2)
            res.singlecv.mse$r2[[t]][f,r] <- round(cor(pred.smse, datasets[[t]]$y)^2,2)
            
          }
        }
      }
      res.rcv.mse$meas <- traitHSR$y
      res.singlecv.mse$meas <- traitHSR$y
      save( res.rcv.mse,res.singlecv.mse, file=paste0("results/unseen_genotypes/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    }

  }
}




