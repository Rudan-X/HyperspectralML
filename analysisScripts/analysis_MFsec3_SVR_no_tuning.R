rm(list = ls())
list.of.packages <- c("pls","dplyr","reshape2","parallel","ggplot2","caret")
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


year <- 2023
plsmethod<-"oscorespls"
maxnoc<-25

segments0 <- 5
nrep <- 20

repl<-10
segment<-3

agg <- "sampledHSR"


vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","maxNPQ") 

datatypes <- c("plot_averaged","genotype_averaged","raw_data") #
for (inVar in vars){
  if (inVar=="N" | inVar=="CN"){
    endd <- 2
  }else{
    endd <- 3
  }
  for (d in 1:endd){
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
    
    load(paste0("results/comparePLSR_SVR/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    folds <- res.rcv.mse$folds
      
    templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
    temp <- matrix(0,segments0,nrep)
    temp2 <- matrix(0,nrow(traitHSR),segments0*nrep)
    res.rcv.svr <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
    
    
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
        
        model.svr <- train(as.formula("y~Spectra"),data = cal.data,method="svmLinear", preProcess=NULL,
                       trControl=trainControl(method="none"))
        
        for (t in 1:2){
          pred.svr<-predict(model.svr, newdata=datasets[[t]])
         
          if (t==1){
            res.rcv.svr$pred[cal.ind,(r-1)*segments0+f] <- pred.svr
          }else{
            res.rcv.svr$pred[val.ind,(r-1)*segments0+f] <- pred.svr
          }

          sq_resid <- (pred.svr - datasets[[t]]$y)^2
          res.rcv.svr$R2[[t]][f,r] <- 1-sum(sq_resid)/sum((datasets[[t]]$y-mean(datasets[[t]]$y))^2)
          rmsep <- round(sqrt(mean(sq_resid)),2)
          res.rcv.svr$RMSEP[[t]][f,r] <- rmsep
          val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
          res.rcv.svr$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
        }
      }
    }
    res.rcv.svr$folds <- folds
    save(res.rcv.svr,file= paste0("results/comparePLSR_SVR/SVR_model_",datatype,"_",agg,year,"_o", segments0*nrep,"_",inVar,".RData"))
  }
}




