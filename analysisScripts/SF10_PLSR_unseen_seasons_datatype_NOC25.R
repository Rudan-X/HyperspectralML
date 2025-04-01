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



segments0 <- 1
nrep <- 1

repl<-10
segment<-3

plsmethod<-"oscorespls"
maxnoc <- 25

agg <- "sampledHSR"
datatypes <- c("plot_averaged","genotype_averaged","raw_data")

vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #
year1 <- c(2022,2023,2021,2023,2021,2022)
year2 <- c(2021,2021,2022,2022,2023,2023)

for (y in 1:6){
  for (inVar in vars){
    if (inVar=="N" | inVar=="CN"){
      endd <- 2
    }else{
      endd <- 3
    }
    for (d in 1:endd){ 
      datatype <- datatypes[d]
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
      
      
      templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
      temp <- matrix(0,segments0,nrep)
      temp2 <- matrix(0,nrow(traitHSR),1)
      res.rcv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
      
      
      val.ind <- which(traitHSR$season==year2[y])
      cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
      val.data <- traitHSR[val.ind, ]
      cal.data <- traitHSR[cal.ind, ]
      cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
      cal.data <- data.frame(y=cal.data$y, Spectra=I(cal_spec))
      
      val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
      val.data <- data.frame(y=val.data$y, Spectra=I(val_spec))
      
      datasets<-list(cal.data,val.data)
      
      #####################
      
      # iloop.rcv <- innerloop_repeatedCV(dataset=cal.data,maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
      # 
      # nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
      # 
      res.rcv.mse$bestnoc <- maxnoc
      noc.rcv.mse <- maxnoc
      model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.rcv.mse, plsmethod = method, scale = FALSE,validation = "none")
      
      for (t in 1:2){
        pred.mse<-predict(model.mse, newdata=datasets[[t]], ncomp = noc.rcv.mse)
        
        if (t==1){
          res.rcv.mse$pred[cal.ind,] <- pred.mse
        }else{
          res.rcv.mse$pred[val.ind,] <- pred.mse
        }
        
        sq_resid <- (pred.mse - datasets[[t]]$y)^2
        res.rcv.mse$R2[[t]] <- 1-sum(sq_resid)/sum((datasets[[t]]$y-mean(datasets[[t]]$y))^2)
        rmsep <- round(sqrt(mean(sq_resid)),2)
        res.rcv.mse$RMSEP[[t]] <- rmsep
        val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
        res.rcv.mse$RMSEPper[[t]] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
      }
      res.rcv.mse$meas <- traitHSR$y
      save(val.ind,res.rcv.mse,file= paste0("results/unseen_seasons/PLSR_model_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_NOC25_",inVar,".RData"))
    }
  }
}

################################# COMBINATION OF PAIRS ##############################################
vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #
year1 <- list()
year1[[1]] <- c(2022,2023)
year1[[2]] <- c(2021,2023)
year1[[3]] <- c(2021,2022)

year2 <- c(2021,2022,2023)
for (y in 1:3){
  for (inVar in vars){
    if (inVar=="N" | inVar=="CN"){
      endd <- 2
    }else{
      endd <- 3
    }
    for (d in 1:endd){ 
      datatype <- datatypes[d]
      traitHSR1a <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[[y]][1],"_",datatype,".csv"))
      traitHSR1b <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[[y]][2],"_",datatype,".csv"))
      sametrait <- intersect(colnames(traitHSR1a),colnames(traitHSR1b))
      
      traitHSR1 <- rbind(traitHSR1a[,sametrait],traitHSR1b[,sametrait])
      traitHSR1$Season <- paste0(year1[[y]][1],"&",year1[[y]][2])
      
      
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
      
      
      templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
      temp <- matrix(0,segments0,nrep)
      temp2 <- matrix(0,nrow(traitHSR),1)
      res.rcv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
      
      
      val.ind <- which(traitHSR$season==year2[y])
      cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
      val.data <- traitHSR[val.ind, ]
      cal.data <- traitHSR[cal.ind, ]
      cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
      cal.data <- data.frame(y=cal.data$y, Spectra=I(cal_spec))
      
      val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
      val.data <- data.frame(y=val.data$y, Spectra=I(val_spec))
      
      datasets<-list(cal.data,val.data)
      
      #####################
      
      # iloop.rcv <- innerloop_repeatedCV(dataset=cal.data,maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
      # 
      # nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
      # 
      res.rcv.mse$bestnoc <- maxnoc
      noc.rcv.mse <- maxnoc
      model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.rcv.mse, plsmethod = method, scale = FALSE,validation = "none")
      
      
      for (t in 1:2){
        pred.mse<-predict(model.mse, newdata=datasets[[t]], ncomp = noc.rcv.mse)
        
        if (t==1){
          res.rcv.mse$pred[cal.ind,] <- pred.mse
        }else{
          res.rcv.mse$pred[val.ind,] <- pred.mse
        }
        
        sq_resid <- (pred.mse - datasets[[t]]$y)^2
        res.rcv.mse$R2[[t]] <- 1-sum(sq_resid)/sum((datasets[[t]]$y-mean(datasets[[t]]$y))^2)
        rmsep <- round(sqrt(mean(sq_resid)),2)
        res.rcv.mse$RMSEP[[t]] <- rmsep
        val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
        res.rcv.mse$RMSEPper[[t]] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
      }
      res.rcv.mse$meas <- traitHSR$y
      save(val.ind,res.rcv.mse,file= paste0("results/unseen_seasons/PLSR_model_",datatype,"_",agg,"_training",unique(traitHSR1$Season),"_testing",year2[y],"_NOC25_",inVar,".RData"))
    }
  }
}
