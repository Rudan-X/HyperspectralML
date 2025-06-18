
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
maxnoc<-40

agg <- "sampledHSR"

vars<-c( "SLA","N", "CN", "d13C",  "Vpmax","Vmax", "a400", "gsw",  
         "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate","maxNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")

dts <- c("raw_data","plot_averaged","genotype_averaged")
dtsi<-c( rep(3,3), 2,rep(1,4), rep(1,6),rep(2,3))
MLs <- c(rep("PLSR",3),rep("SVR",5),rep("PLSR",9))


segments0 <- 5
nrep <- 20

segment<-3
repl<-1

year_train <- c(2022,2023,2021,2023,2021,2022)
year_test <- c(2021,2021,2022,2022,2023,2023)

for (y in 1:6){
  for (i in 1:length(vars)){
    inVar <- vars[i]
    datatype <- dts[dtsi[i]]
    
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
    traitHSR <- remove_outliers(inVar, traitHSR, datatype)
    
    
    traitHSR_train <- traitHSR[traitHSR$season==year_train[y],]
    traitHSR_test <- traitHSR[traitHSR$season==year_test[y],]
    
    templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
    temp <- matrix(0,segments0,nrep)
    
    
    if (MLs[i]=="PLSR"){
      res.plsr <- list(R2 = templ, cor2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, 
                          pred_test = matrix(0,nrow(traitHSR_test),segments0*nrep),
                          pred_train = matrix(0,nrow(traitHSR_train),segments0*nrep))
    }else if (MLs[i]=="SVR"){
      
      res.svr <- list(R2 = templ, cor2 = templ, RMSEP = templ, RMSEPper = templ, 
                          pred_test = matrix(0,nrow(traitHSR_test),segments0*nrep),
                          pred_train = matrix(0,nrow(traitHSR_train),segments0*nrep))
    }
    
    train_folds <- list()
    test_folds <- list()
    for (r in 1:nrep){
      print(paste0("Testing",year_test[y],", Training",year_train[y], ", Var:  ",inVar, " ,Current repetition: ",r))
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
        
        ####################
        if (MLs[i]=="PLSR"){
          iloop.singlecv <- innerloop_repeatedCV_genotype_season(dataset=cal.data,maxComps =maxnoc,method=plsmethod,iterations=repl,segments=segment,parallel=FALSE)
          

          nComps.singlemse<-selectNOC_MSE_original(MSEP = iloop.singlecv$MSEP, select_strat ="hastie",sdfact=0.5, segments = segment)
          
          res.plsr$bestnoc[f,r] <- nComps.singlemse$bestnoc
          noc.singlecv.mse <- nComps.singlemse$bestnoc
          
          model.smse <- plsr(as.formula("y~Spectra"),data=cal.data, ncomp = noc.singlecv.mse, method = plsmethod, scale = FALSE,validation = "none")
          
          for (t in 1:2){
            
            pred.smse<-predict(model.smse, newdata=datasets[[t]], ncomp = noc.singlecv.mse)
            if (t==1){
              res.plsr$pred_train[cal.ind,(r-1)*segments0+f] <- pred.smse
            }else{
              res.plsr$pred_test[val.ind,(r-1)*segments0+f] <- pred.smse
            }
            
            sq_resid <- (pred.smse - datasets[[t]]$y)^2
            res.plsr$R2[[t]][f,r] <- round(pls::R2(model.smse,newdata=datasets[[t]],intercept=F)[[1]][noc.singlecv.mse],2)
            res.plsr$cor2[[t]][f,r] <- round(cor(pred.smse,datasets[[t]]$y)^2,2)
            rmsep <- round(sqrt(mean(sq_resid)),2)
            res.plsr$RMSEP[[t]][f,r] <- rmsep
            val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
            res.plsr$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
          }
        }else if (MLs[i]=="SVR"){
          model.svr <- train(as.formula("y~Spectra"),data = cal.data,method="svmLinear", preProcess=NULL,
                             trControl=trainControl(method="none"))
          
          for (t in 1:2){
            pred.svr<-predict(model.svr, newdata=datasets[[t]])
            
            if (t==1){
              res.svr$pred_train[cal.ind,(r-1)*segments0+f] <- pred.svr
            }else{
              res.svr$pred_test[val.ind,(r-1)*segments0+f] <- pred.svr
            }
            
            sq_resid <- (pred.svr - datasets[[t]]$y)^2
            res.svr$R2[[t]][f,r] <- 1-sum(sq_resid)/sum((datasets[[t]]$y-mean(datasets[[t]]$y))^2)
            res.svr$cor2[[t]][f,r] <- round(cor(pred.svr,datasets[[t]]$y)^2,2)
            rmsep <- round(sqrt(mean(sq_resid)),2)
            res.svr$RMSEP[[t]][f,r] <- rmsep
            val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
            res.svr$RMSEPper[[t]][f,r] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
          }
        }
      }
    }
    
    if (MLs[i]=="PLSR"){
      res.plsr$train_folds <- train_folds
      res.plsr$test_folds <- test_folds
      res.plsr$measured <- traitHSR$y
      filen <- paste0("results/unseen_seasons_and_genotypes_specific/PLSR_model_",datatype,"_",agg,"_training",year_train[y],"_testing",year_test[y],"_o", segments0*nrep,"_i",segment*repl,"_",inVar,".RData")
      save(res.plsr, file=filen)
      
    }else if (MLs[i]=="SVR"){
      
      res.svr$train_folds <- train_folds
      res.svr$test_folds <- test_folds
      res.svr$measured <- traitHSR$y
      filen <- paste0("results/unseen_seasons_and_genotypes_specific/SVR_model_",datatype,"_",agg,"_training",year_train[y],"_testing",year_test[y],"_o", segments0*nrep,"_i",segment*repl,"_",inVar,".RData")
      save(res.svr, file=filen)
    }
  }
}


