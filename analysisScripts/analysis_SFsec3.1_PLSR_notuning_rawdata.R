
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

vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","maxNPQ") #

# "SLA","N","CN","Vpmax","Vmax","gsw","a400","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm"

year <- 2021

plsmethod<-"oscorespls"
maxnoc<-15

segments0 <- 5
nrep <- 20
repl<-10
segment<-3

agg <-"sampledHSR"


for (datatype in c("raw_data")){
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
    
    load(paste0("results/comparePLSR_SVR/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    folds <- res.rcv.mse$folds
    
    templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
    temp <- matrix(0,segments0,nrep)
    temp2 <- matrix(0,nrow(traitHSR),segments0*nrep)
    res.rcv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
    
    for (r in 1:nrep){
      print(paste0("Current repetition: ",r))
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
        
        
        model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, plsmethod = method, scale = FALSE,validation = "none")
        noc.rcv.mse <- model.mse$ncomp
        
        for (t in 1:2){
         
          pred.mse<-predict(model.mse, newdata=datasets[[t]], ncomp = noc.rcv.mse)

          if (t==1){
            res.rcv.mse$pred[cal.ind,(r-1)*segments0+f] <- pred.mse
          }else{
            res.rcv.mse$pred[val.ind,(r-1)*segments0+f] <- pred.mse
          }
          
          
          sq_resid <- (pred.mse - datasets[[t]]$y)^2
          res.rcv.mse$R2[[t]][f,r] <- round(pls::R2(model.mse,newdata=datasets[[t]],intercept=F)[[1]][noc.rcv.mse],2)
          rmsep <- round(sqrt(mean(sq_resid)),2)
          res.rcv.mse$RMSEP[[t]][f,r] <- rmsep
          val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
          res.rcv.mse$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
        
        }
      }
    }
    res.rcv.mse$folds <- folds

    noc.rcv.mse <- "default"
    save( res.rcv.mse, file=paste0("results/comparePLSR_SVR/PLSR_notuned_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_NOC",noc.rcv.mse,"_",inVar,".RData"))
  }
}



      




