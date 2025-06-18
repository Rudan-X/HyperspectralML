
list.of.packages <- c("pls","dplyr","plotrix","ggplot2","reshape2")
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

# "C", "CN", "Vpmax", "Vpmax.Vmax", "Area", "d13C", "d15N", "SL", "gsw",  "iWUE", "A_sat", "gs_sat","ci_sat", "iWUE_sat",   "gs", "sl",
# vars <- c("phiPSII_ind_res","maxNPQ", 
#           "NPQ_ind_amp", "NPQ_ind_rate", "NPQ_rel_amp", "NPQ_rel_rate", "NPQ_rel_res",  "phiPSII_ind_amp",  "phiPSII_ind_rate", "NPQ_ind_linear",  
#           "endNPQ", "endFvFm",  "initialFvFm")

vars <- c("d13C")

plsmethod<-"oscorespls"
maxnoc<-40

segments0 <- 5
nrep <- 20
repl<-10
segment<-3

agg <-"sampledHSR"

for (year in c(2023)){
  for (datatype in c("plot_averaged")){ #"genotype_averaged","plot_averaged"
    for (inVar in vars){
      
      traitHSR <- read.csv(paste0("data/combined_data/traits_and_HSR",year,"_",datatype,".csv"))
      
      if (inVar%in%colnames(traitHSR)){
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
        res.rcv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
        res.singlecv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
        
        folds <- list()
        for (r in 1:nrep){
          print(paste0(inVar," Current repetition: ",r))
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
            start_time <- Sys.time()
            iloop.rcv <- innerloop_repeatedCV(dataset=datasets[[1]],formula=as.formula("y~Spectra"),maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
            end_time <- Sys.time()
            end_time - start_time
            
            iloop.singlecv <- innerloop_repeatedCV(dataset=cal.data,formula=as.formula("y~Spectra"),maxComps =maxnoc,method=plsmethod,iterations=1,segments=segment)
            nComps.singlemse<-selectNOC_MSE_original(MSEP = iloop.singlecv$MSEP, select_strat ="hastie",sdfact=0.5, segments = segment)
            
            nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
            
            res.rcv.mse$bestnoc[f,r] <- nComps.mse$bestnoc
            noc.rcv.mse <- nComps.mse$bestnoc
            
            noc.singlecv.mse <- nComps.singlemse$bestnoc
            res.singlecv.mse$bestnoc[f,r] <- nComps.singlemse$bestnoc # not used
            
            model.mse <- plsr(as.formula("y~Spectra"),data=datasets[[1]], ncomp = noc.rcv.mse, method = plsmethod, scale = FALSE,validation = "none")
            model.smse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.singlecv.mse,  method = plsmethod, scale = FALSE,validation = "none")
            
            for (t in 1:2){
              pred.mse<-predict(model.mse, newdata=datasets[[t]], ncomp = noc.rcv.mse)
              pred.smse<-predict(model.smse, newdata=datasets[[t]], ncomp = noc.singlecv.mse)
              if (t==1){
                
                res.rcv.mse$pred[cal.ind,(r-1)*segments0+f] <- pred.mse
                res.singlecv.mse$pred[cal.ind,(r-1)*segments0+f] <- pred.smse
              }else{
                
                res.rcv.mse$pred[val.ind,(r-1)*segments0+f] <- pred.mse
                res.singlecv.mse$pred[val.ind,(r-1)*segments0+f] <- pred.smse
              }
              
              
              sq_resid <- (pred.mse - datasets[[t]]$y)^2
              res.rcv.mse$R2[[t]][f,r] <- round(pls::R2(model.mse,newdata=datasets[[t]],intercept=F)[[1]][noc.rcv.mse],2)
              rmsep <- round(sqrt(mean(sq_resid)),2)
              res.rcv.mse$RMSEP[[t]][f,r] <- rmsep
              val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
              res.rcv.mse$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
              
              sq_resid <- (pred.smse - datasets[[t]]$y)^2
              res.singlecv.mse$R2[[t]][f,r] <- round(pls::R2(model.smse,newdata=datasets[[t]],intercept=F)[[1]][noc.singlecv.mse],2)
              rmsep <- round(sqrt(mean(sq_resid)),2)
              res.singlecv.mse$RMSEP[[t]][f,r] <- rmsep
              val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
              res.singlecv.mse$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
            }
          }
        }
        
        res.rcv.mse$folds <- folds
        
        save(res.rcv.mse,res.singlecv.mse, file=paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_",agg,year,"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      }
    }
  }
}


      




