
list.of.packages <- c("pls","dplyr","plotrix","ggplot2","gridExtra","chemometrics","reshape2")
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

vars<-c("SLA","N") #, "Vpmax","Vmax","a400","phiPSII_ind_res","maxNPQ"

# "SLA","N","CN","Vpmax","Vmax","gsw","a400","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm"

year <- 2021

plsmethod<-"oscorespls"
maxnoc<-25

segments0 <- 5
nrep <- 20

npermutations <- 90
innerprop <- 0.7
# repl<-10
# segment<-3
# load genotype-averaged data
agg <- c("fullHSR","sampledHSR") #
datatype <- "raw_data"
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
    
    templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
    temp <- matrix(0,segments0,nrep)
    temp2 <- matrix(0,nrow(traitHSR),segments0*nrep)
    res.perm.press <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
    res.perm.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
    
    folds <- list()
    for (r in 1:nrep){
      print(paste0("Current repetition: ",r))
      folds[[r]] <- random_split(seq(1,nrow(traitHSR)), k = segments0)
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
        # start_time <- Sys.time()
        iloop.perm <- innerloop_permutations(dataset=cal.data, maxComps=maxnoc, method=plsmethod, iterations=npermutations, prop=innerprop)
        # end_time <- Sys.time()
        # end_time - start_time
        
        nComps.press<-selectNOC_PRESS(PRESS=iloop.perm$PRESS,method = "firstMin",maxComps = maxnoc)
        nComps.mse<-selectNOC_MSE_original(MSEP = iloop.perm$MSEP, select_strat ="hastie",sdfact=0.5,segments = 1/(1-innerprop)) 
        
        res.perm.press$bestnoc[f,r] <- nComps.press$bestnoc
        noc.perm.press <- nComps.press$bestnoc
        
        res.perm.mse$bestnoc[f,r] <- nComps.mse$bestnoc
        noc.perm.mse <- nComps.mse$bestnoc
        
        model.press <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.perm.press, method = plsmethod, scale = FALSE, validation = "none")
        model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.perm.mse, method = plsmethod, scale = FALSE, validation = "none")
        
        for (t in 1:2){
          pred.press<-predict(model.press, newdata=datasets[[t]], ncomp = noc.perm.press)
          pred.press2<-predict(model.press2, newdata=datasets[[t]], ncomp = noc.perm.press)
          pred.mse<-predict(model.mse, newdata=datasets[[t]], ncomp = noc.perm.mse)
         
          if (t==1){
            res.perm.press$pred[cal.ind,(r-1)*segments0+f] <- pred.press
            res.perm.mse$pred[cal.ind,(r-1)*segments0+f] <- pred.mse
          }else{
            res.perm.press$pred[val.ind,(r-1)*segments0+f] <- pred.press
            res.perm.mse$pred[val.ind,(r-1)*segments0+f] <- pred.mse
          }
          
          sq_resid <- (pred.press - datasets[[t]]$y)^2
          res.perm.press$R2[[t]][f,r] <- round(pls::R2(model.press,newdata=datasets[[t]],intercept=F)[[1]][noc.perm.press],2)
          rmsep <- round(sqrt(mean(sq_resid)),2)
          res.perm.press$RMSEP[[t]][f,r] <- rmsep
          val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
          res.perm.press$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
          
          sq_resid <- (pred.mse - datasets[[t]]$y)^2
          res.perm.mse$R2[[t]][f,r] <- round(pls::R2(model.mse,newdata=datasets[[t]],intercept=F)[[1]][noc.perm.mse],2)
          rmsep <- round(sqrt(mean(sq_resid)),2)
          res.perm.mse$RMSEP[[t]][f,r] <- rmsep
          val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
          res.perm.mse$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
        }
      }
    }
    res.perm.press$folds <- folds
    res.perm.mse$folds <- folds
   
    save(res.perm.press, res.perm.mse,file=paste0("results/compareMSE_PRESS/permutation_model_",datatype,"_",agg,year,"_o", segments0*nrep,"per",npermutations,"_",inVar,".RData"))
  }
}

      




