
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

vars<-c( "SLA","N","Vmax","a400","phiPSII_ind_res","maxNPQ","C", "CN", "d13C", "d15N",  "Vpmax",  "SL", "gsw",  "iWUE", "NPQ_ind_amp",
         "NPQ_ind_rate", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res",  "phiPSII_ind_amp",  "phiPSII_ind_rate", "NPQ_ind_linear",  "endNPQ", "endFvFm",  "initialFvFm") #


year <- 2021

plsmethod<-"oscorespls"
maxnoc<-40

segments0 <- 5
nrep <- 20

# load genotype-averaged data
aggHSR <- c("fullHSR","sampledHSR") #
datatype <- "raw_data"



for (inVar in vars){
  traitHSR <- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_",datatype,".csv"))

  if (inVar=="N" | inVar=="CN"){
    Start.wave <- 1500
    End.wave <- 2400
  }else{
    Start.wave <- 400
    End.wave <- 2400
  }

  wv <- seq(Start.wave,End.wave,1)

  spec <- as.matrix(traitHSR[, which(names(traitHSR) %in%  paste0("mean_",wv))])
  traitHSR <- cbind(data.frame(y=traitHSR[,inVar]),spec)
  colnames(traitHSR)[1] <- "y"
  traitHSR <- traitHSR[!is.na(traitHSR$y),]
  
  traitHSR <- remove_outliers(inVar,traitHSR,datatype)
  
  folds <- list()
  for (r in 1:nrep){
    print(paste0("Trait: ",inVar,", Current repetition: ",r))
    folds[[r]] <- random_split(seq(1,nrow(traitHSR)), k = segments0)
  }
  save(folds, file=paste0("results/model_sensitivity/fold_split_","_",inVar,".RData"))
}


for (agg in aggHSR){
  for (inVar in vars){
    
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
    }else{
      wv <- seq(Start.wave,End.wave,by=5)
    }
    
    spec <- as.matrix(traitHSR[, which(names(traitHSR) %in%  paste0("mean_",wv))])
    traitHSR <- cbind(data.frame(y=traitHSR[,inVar]),spec)
    colnames(traitHSR)[1] <- "y"
    traitHSR <- traitHSR[!is.na(traitHSR$y),]
    traitHSR <- remove_outliers(inVar,traitHSR,datatype)
    
    
    templ <- list(train = array(0,c(segments0,nrep,maxnoc)), test = array(0,c(segments0,nrep,maxnoc)))
    temp <- array(0, c(segments0,nrep,maxnoc))
    temp2 <- array(0,c(nrow(traitHSR),segments0*nrep,maxnoc))
    temp3 <- array(0, c(segments0,nrep))
    res.singlecv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp3, pred = temp2)
    
    
    res.singlecv.mse$MSEP <- templ
    res.singlecv.mse$PRESS <- templ
    res.singlecv.mse$bestnoc_MSEP <- matrix(0, 1,nrep)
    res.singlecv.mse$bestnoc_PRESS <- matrix(0, 1,nrep)
    
    
    load(paste0("results/model_sensitivity/fold_split_","_",inVar,".RData"))
    for (r in 1:nrep){
      print(paste0("Trait: ",inVar,", Current repetition: ",r))
      for (f in 1:length(folds[[r]])){
        val.ind <- folds[[r]][[f]]
        cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
        val.data <- traitHSR[val.ind, ]
        cal.data <- traitHSR[cal.ind, ]
        cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
        cal.data <- data.frame(y=cal.data$y, Spectra=I(cal_spec))
        
        val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
        val.data <- data.frame(y=val.data$y, Spectra=I(val_spec))
        
        datasets<-list(cal.data,val.data)
        
        
        #####################
        
        model.smse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = maxnoc,  method = plsmethod, scale = FALSE,validation = "none")
        
        for (t in 1:2){
          pred.smse<-predict(model.smse, newdata=datasets[[t]]) #, ncomp = noc.singlecv.mse
          if (t==1){
            res.singlecv.mse$pred[cal.ind,(r-1)*segments0+f,] <- pred.smse
          }else{
            res.singlecv.mse$pred[val.ind,(r-1)*segments0+f,] <- pred.smse
          }

          sq_resid <- (pred.smse - datasets[[t]]$y)^2
          res.singlecv.mse$R2[[t]][f,r,] <- round(pls::R2(model.smse,newdata=datasets[[t]],intercept=F)[[1]][1,1,],2)
          
          res.singlecv.mse$MSEP[[t]][f,r,] <- apply(sq_resid,3,mean)
          res.singlecv.mse$PRESS[[t]][f,r,] <- apply(sq_resid, 3, sum)
          
          rmsep <- sqrt(apply(sq_resid,3,mean))
          res.singlecv.mse$RMSEP[[t]][f,r,] <- rmsep
          val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
          res.singlecv.mse$RMSEPper[[t]][f,r,] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
        }
        
      }
      res.singlecv.mse$bestnoc_MSEP[1,r]<-selectNOC_MSE_original(MSEP = res.singlecv.mse$MSEP[[2]][,r,], select_strat ="hastie",sdfact=0.5, segments = segments0)$bestnoc
    }
    est_noc <- selectNOC_MSE(MSEP = reshape_mat(res.singlecv.mse$RMSEP[[2]]), select_strat ="hastie",sdfact=0.5,repl = nrep, segments = segments0) 
    
    res.singlecv.mse$bestnoc_final <- est_noc$bestnoc
    save(res.singlecv.mse, file=paste0("results/model_sensitivity/PLSR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
  }
}

      




