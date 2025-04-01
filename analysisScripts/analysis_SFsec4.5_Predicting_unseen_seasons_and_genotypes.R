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



segments0 <- 5
nrep <- 10

repl<-10
segment<-3

plsmethod<-"oscorespls"
maxnoc <- 25

agg <- "sampledHSR"
datatypes <- c("plot_averaged","genotype_averaged","raw_data")

vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #
year_train <- c(2022,2023,2021,2023,2021,2022)
year_test <- c(2021,2021,2022,2022,2023,2023)

for (y in 2:2){
  for (inVar in vars){
    if (inVar=="N" | inVar=="CN"){
      endd <- 2
    }else{
      endd <- 3
    }
    for (d in 1:endd){ 
      datatype <- datatypes[d]
      traitHSR_train <- read.csv(paste0("data/combined_data/traits_and_HSR",year_train[y],"_",datatype,".csv"))
      traitHSR_train$Season <- year_train[y]
      traitHSR_test <- read.csv(paste0("data/combined_data/traits_and_HSR",year_test[y],"_",datatype,".csv"))
      traitHSR_test$Season <- year_test[y]
      
      uniq_acc <- unique(c(traitHSR_train$Accession,traitHSR_test$Accession))
      sametrait <- intersect(colnames(traitHSR_train),colnames(traitHSR_test))
      traitHSR <- rbind(traitHSR_train[,sametrait],traitHSR_test[,sametrait])
      
      
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
      
      traitHSR_train <- traitHSR[traitHSR$season==year_train[y],]
      traitHSR_test <- traitHSR[traitHSR$season==year_test[y],]
      
      templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
      temp <- matrix(0,segments0,nrep)
      res.rcv.mse <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, 
                          pred_test = matrix(0,nrow(traitHSR_test),segments0*nrep),
                          pred_train = matrix(0,nrow(traitHSR_train),segments0*nrep))
      
      train_folds <- list()
      test_folds <- list()
      for (r in 1:nrep){
        print(paste0("Current repetition: ",r))
        folds <- genotype_season_split(test_dataset = traitHSR_test,train_dataset = traitHSR_train, k = segments0)
        train_folds[[r]] <- folds$trainfold
        test_folds[[r]] <- folds$testfold
        for (f in 1:length(train_folds[[r]])){
          val.ind <- test_folds[[r]][[f]]
          cal.ind <- train_folds[[r]][[f]]
          val.data <- traitHSR_test[val.ind, ]
          cal.data <- traitHSR_train[cal.ind, ]
          cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
          cal.data <- data.frame(genotype=cal.data$genotype,y=cal.data$y, Spectra=I(cal_spec))
          
          val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
          val.data <- data.frame(genotype=val.data$genotype,y=val.data$y, Spectra=I(val_spec))
          
          datasets<-list(cal.data,val.data)
          
          #####################
          iloop.rcv <- innerloop_repeatedCV_genotype_season(dataset=cal.data,maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
        
          nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
          
          res.rcv.mse$bestnoc[f,r] <- nComps.mse$bestnoc
          noc.rcv.mse <- nComps.mse$bestnoc
          
          model.mse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.rcv.mse, plsmethod = method, scale = FALSE,validation = "none")
          
          for (t in 1:2){
            pred.mse<-predict(model.mse, newdata=datasets[[t]], ncomp = noc.rcv.mse)
            if (t==1){
              res.rcv.mse$pred_train[cal.ind,(r-1)*segments0+f] <- pred.mse
            }else{
              res.rcv.mse$pred_test[val.ind,(r-1)*segments0+f] <- pred.mse
            }
            
            sq_resid <- (pred.mse - datasets[[t]]$y)^2
            res.rcv.mse$R2[[t]][f,r] <- round(cor(pred.mse,datasets[[t]]$y)^2,2)
            rmsep <- round(sqrt(mean(sq_resid)),2)
            res.rcv.mse$RMSEP[[t]][f,r] <- rmsep
            val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
            res.rcv.mse$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
          }
        }
      }
      res.rcv.mse$meas <- traitHSR$y
      save(res.rcv.mse,file= paste0("results/unseen_seasons_and_genotypes/PLSR_model_",datatype,"_",agg,"_training",year_train[y],"_testing",year_test[y],"_i",segment*repl,"_",inVar,".RData"))
    }
  }
}

################################# COMBINATION OF PAIRS ##############################################
vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #
year_train <- list()
year_train[[1]] <- c(2022,2023)
year_train[[2]] <- c(2021,2023)
year_train[[3]] <- c(2021,2022)

year_test <- c(2021,2022,2023)
for (y in 1:3){
  for (inVar in vars){
    if (inVar=="N" | inVar=="CN"){
      endd <- 2
    }else{
      endd <- 3
    }
    for (d in 1:endd){ 
      datatype <- datatypes[d]
      traitHSR1a <- read.csv(paste0("data/combined_data/traits_and_HSR",year_train[[y]][1],"_",datatype,".csv"))
      traitHSR1b <- read.csv(paste0("data/combined_data/traits_and_HSR",year_train[[y]][2],"_",datatype,".csv"))
      sametrait <- intersect(colnames(traitHSR1a),colnames(traitHSR1b))
      
      traitHSR1 <- rbind(traitHSR1a[,sametrait],traitHSR1b[,sametrait])
      traitHSR1$Season <- paste0(year_train[[y]][1],"&",year_train[[y]][2])
      
      
      traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR",year_test[y],"_",datatype,".csv"))
      traitHSR2$Season <- year_test[y]
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
      
      
      val.ind <- which(traitHSR$season==year_test[y])
      cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
      val.data <- traitHSR[val.ind, ]
      cal.data <- traitHSR[cal.ind, ]
      cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
      cal.data <- data.frame(y=cal.data$y, Spectra=I(cal_spec))
      
      val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
      val.data <- data.frame(y=val.data$y, Spectra=I(val_spec))
      
      datasets<-list(cal.data,val.data)
      
      #####################
      
      iloop.rcv <- innerloop_repeatedCV(dataset=cal.data,maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
      
      nComps.mse<-selectNOC_MSE(MSEP = iloop.rcv$MSEP, select_strat ="hastie",sdfact=0.5,repl = repl, segments = segment) 
      
      res.rcv.mse$bestnoc <- nComps.mse$bestnoc
      noc.rcv.mse <- nComps.mse$bestnoc
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
      save(val.ind,res.rcv.mse,file= paste0("results/unseen_seasons/PLSR_model_",datatype,"_",agg,"_training",unique(traitHSR1$Season),"_testing",year_test[y],"_i",segment*repl,"_",inVar,".RData"))
    }
  }
}
