
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

vars<-c("Vmax","a400","phiPSII_ind_res","maxNPQ") #"SLA","N",

year <- 2023
plsmethod<-"oscorespls"
maxnoc<-25

segments0 <- 5
nrep <- 20
repl<-10
segment<-3

agg <-"sampledHSR"
year1 <- c(2021,2021,2022)
year2 <- c(2022,2023,2023)
for (y in 1:3){
  for (datatype in c("plot_averaged")){ #,"plot_averaged","genotype_averaged"
    for (inVar in vars){
      
      traitHSR1 <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[y],"_",datatype,".csv"))
      traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR",year2[y],"_",datatype,".csv"))

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
          
          iloop.rcv <- innerloop_repeatedCV(dataset=cal.data,formula=as.formula("y~Spectra"),maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
          
          nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
          
          res.rcv.mse$bestnoc[f,r] <- nComps.mse$bestnoc
          noc.rcv.mse <- nComps.mse$bestnoc
          
          
          model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.rcv.mse,  method = plsmethod, scale = FALSE,validation = "none")
          
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
      # save(res.rcv.mse, file=paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      save(res.rcv.mse, file=paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_",agg,year1[y],"and",year2[y],"_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
    }
  }
}


########################################################################################################################


for (y in 1:1){
  for (datatype in c("plot_averaged")){ #,"plot_averaged","genotype_averaged"
    for (inVar in vars){

      traitHSR1 <- read.csv(paste0("data/combined_data/traits_and_HSR2021_",datatype,".csv"))
      traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR2022_",datatype,".csv"))
      traitHSR3 <- read.csv(paste0("data/combined_data/traits_and_HSR2023_",datatype,".csv"))
      
      sametrait <- intersect(colnames(traitHSR1),colnames(traitHSR2))
      sametrait <- intersect(sametrait,colnames(traitHSR3))
      traitHSR <- rbind(traitHSR1[,sametrait],traitHSR2[,sametrait])
      traitHSR <- rbind(traitHSR,traitHSR3[,sametrait])
      
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
          
          iloop.rcv <- innerloop_repeatedCV(dataset=cal.data,formula=as.formula("y~Spectra"),maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
          
          nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
          
          res.rcv.mse$bestnoc[f,r] <- nComps.mse$bestnoc
          noc.rcv.mse <- nComps.mse$bestnoc
          
          
          model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.rcv.mse,  method = plsmethod, scale = FALSE,validation = "none")
          
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
      save(res.rcv.mse, file=paste0("results/compare_season/PLSR_repeatedCV_model_",datatype,"_allyears_o", segments0*nrep,"i",segment*repl,"_",inVar,".RData"))
      
    }
  }
}



