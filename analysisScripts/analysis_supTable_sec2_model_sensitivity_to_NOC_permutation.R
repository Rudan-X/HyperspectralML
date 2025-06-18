
list.of.packages <- c("pls","reshape2","ggplot2")
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

# vars<-c( "SLA","N","Vmax","a400","phiPSII_ind_res","maxNPQ","C", "CN", "d13C", "d15N",  "Vpmax",  "SL", "gsw",  "iWUE", "NPQ_ind_amp",
#          "NPQ_ind_rate", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res",  "phiPSII_ind_amp",  "phiPSII_ind_rate", "NPQ_ind_linear",  "endNPQ", "endFvFm",  "initialFvFm") #

vars<-c("endNPQ", "NPQ_ind_rate",  "NPQ_rel_res") #

year <- 2021

plsmethod<-"oscorespls"
maxnoc<-40

npermutations <- 100
innerprop <- 0.8
# load genotype-averaged data
aggHSR <- c("fullHSR","sampledHSR") #
datatype <- "raw_data"


for (inVar in vars){
  traitHSR <- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_",datatype,".csv"))
  print(paste0("Trait: ",inVar))
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
  for (r in 1:npermutations){
    folds[[r]] <- sample(1:nrow(traitHSR),floor(innerprop*nrow(traitHSR)))
  }
  save(folds, file=paste0("results/model_sensitivity/permutation_split_",inVar,".RData"))
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
    
    templ <- list(train = array(0,c(npermutations,maxnoc)), test = array(0,c(npermutations,maxnoc)))
    temp <- array(0, c(npermutations,maxnoc))
    temp2 <- array(0,c(nrow(traitHSR),npermutations,maxnoc))

    res.press <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, PRESS = templ, pred = temp2)
    
    
    load(file=paste0("results/model_sensitivity/permutation_split_",inVar,".RData"))
    for (r in 1:npermutations){
      print(paste0("Trait: ",inVar,", Current repetition: ",r))
      cal.ind <- folds[[r]]
      val.ind <- setdiff(seq(1,nrow(traitHSR)),cal.ind)
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
          res.press$pred[cal.ind,r,] <- pred.smse
        }else{
          res.press$pred[val.ind,r,] <- pred.smse
        }

        sq_resid <- (pred.smse - datasets[[t]]$y)^2
        res.press$R2[[t]][r,] <- round(pls::R2(model.smse,newdata=datasets[[t]],intercept=F)[[1]][1,1,],2)
        
        res.press$PRESS[[t]][r,] <- apply(sq_resid, 3, sum)
        
        rmsep <- sqrt(apply(sq_resid,3,mean))
        res.press$RMSEP[[t]][r,] <- rmsep
        val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
        res.press$RMSEPper[[t]][r,] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
      }
    }
    res.press$bestnoc <- selectNOC_PRESS(PRESS = res.press$PRESS[[2]], method = "firstMin",maxComps = maxnoc)$bestnoc
    
    save(res.press, file=paste0("results/model_sensitivity/PLSR_model_permutation_",datatype,"_",agg,year,"_",inVar,".RData"))
  }
}

      




