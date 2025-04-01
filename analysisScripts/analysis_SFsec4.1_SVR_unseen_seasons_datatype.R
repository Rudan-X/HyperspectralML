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
      res.rcv.svr <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
      
      
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
      
      # # trans <- BoxCoxTrans(datasets[[1]]$y)
      # # datasets[[1]]$y_trans <- predict(trans, datasets[[1]]$y)
      datasets[[1]]$y_trans <- log(datasets[[1]]$y)
      # 
      # # trans <- BoxCoxTrans(datasets[[2]]$y)
      # # datasets[[2]]$y_trans <-predict(trans, datasets[[2]]$y)
      datasets[[2]]$y_trans <- log(datasets[[2]]$y)
      
      MLmethod <- "svmRadial"
      
      if (any(colnames(datasets[[t]])%in%"y_trans")){
        model.svr <- train(as.formula("y_trans~Spectra"),data = datasets[[1]],method=MLmethod, preProcess = c("center", "scale"),
                           trControl=trainControl(method="none"))
      }else{
        model.svr <- train(as.formula("y~Spectra"),data = datasets[[1]],method=MLmethod, preProcess = c("center", "scale"),
                           trControl=trainControl(method="none"))
      }
      
      
      for (t in 1:2){
        pred.svr<-predict(model.svr, newdata=datasets[[t]])
        if (any(colnames(datasets[[t]])%in%"y_trans")){
          pred.svr<-exp(pred.svr)
        }
        
        if (t==1){
          res.rcv.svr$pred[cal.ind,] <- pred.svr
        }else{
          res.rcv.svr$pred[val.ind,] <- pred.svr
        }
        
        cor(pred.svr , datasets[[t]]$y)^2
        plot(pred.svr , datasets[[t]]$y)
        
        sq_resid <- (pred.svr - datasets[[t]]$y)^2
        res.rcv.svr$R2[[t]] <- 1-sum(sq_resid)/sum((datasets[[t]]$y-mean(datasets[[t]]$y))^2)
        rmsep <- round(sqrt(mean(sq_resid)),2)
        res.rcv.svr$RMSEP[[t]] <- rmsep
        val_data_range <- range(datasets[[t]]$y, na.rm = TRUE)
        res.rcv.svr$RMSEPper[[t]] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
      }
      res.rcv.svr$meas <- traitHSR$y
      save(val.ind,res.rcv.svr,file= paste0("results/unseen_seasons/SVR_model_",datatype,"_",agg,"_training",year1[y],"_testing",year2[y],"_",inVar,".RData"))
    }
  }
}

################################# COMBINATION OF PAIRS ##############################################
vars<-c("SLA","Vmax","a400","phiPSII_ind_res","initialFvFm") #,"N"
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
      traitHSR1a$Season <- year1[[y]][1]
      
      traitHSR1b <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[[y]][2],"_",datatype,".csv"))
      traitHSR1b$Season <- year1[[y]][2]
      
      sametrait <- intersect(colnames(traitHSR1a),colnames(traitHSR1b))
      
      traitHSR1 <- rbind(traitHSR1a[,sametrait],traitHSR1b[,sametrait])
      # traitHSR1$Season <- paste0(year1[[y]][1],"&",year1[[y]][2])
      
      
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
      res.rcv.svr <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)
      
      
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
      if (inVar=="SLA" & year2[[y]]!=2021){
        ind1a <- which(traitHSR$season==year1[[y]][1])
        ind1b <- which(traitHSR$season==year1[[y]][2])
       
        datasets[[1]]$y_trans <- datasets[[1]]$y
        datasets[[2]]$y_trans <- datasets[[2]]$y
        
        datasets[[1]]$y_trans[ind1a] <- scale(datasets[[1]]$y[ind1a])
        datasets[[1]]$y_trans[ind1b] <- scale(datasets[[1]]$y[ind1b])
        datasets[[2]]$y_trans <- scale(datasets[[2]]$y)
      }
      
      
      #####################
      MLmethod <- "svmLinear"
      if (any(colnames(datasets[[1]])%in%"y_trans")){
        model.svr <- train(as.formula("y_trans~Spectra"),data = datasets[[1]],method=MLmethod, preProcess = c("center", "scale"),
                           trControl=trainControl(method="none"))
      }else{
        model.svr <- train(as.formula("y~Spectra"),data = datasets[[1]],method=MLmethod, preProcess = c("center", "scale"),
                           trControl=trainControl(method="none"))
      }
      
      model.svr <- train(as.formula("y~Spectra"),data = cal.data,method="svmLinear", preProcess = c("center", "scale"),
                         trControl=trainControl(method="none"))
      
      for (t in 1:2){
        pred.svr<-predict(model.svr, newdata=datasets[[t]])
        
        if (any(colnames(datasets[[t]])%in%"y_trans")){
          mu <- mean(datasets[[2]]$y)
          sigma <- sd(datasets[[2]]$y)
          y <- pred.svr * sigma + mu
        }
        
        
        if (t==1){
          res.rcv.svr$pred[cal.ind,] <- pred.svr
        }else{
          res.rcv.svr$pred[val.ind,] <- pred.svr
        }
        
        if (any(colnames(datasets[[t]])%in%"y_trans")){
          observed <- datasets[[t]]$y_trans
        }else{
          observed <- datasets[[t]]$y
        }
        sq_resid <- (pred.svr - observed)^2
        res.rcv.svr$R2[[t]] <- 1-sum(sq_resid)/sum((observed-mean(observed))^2)
        rmsep <- round(sqrt(mean(sq_resid)),2)
        res.rcv.svr$RMSEP[[t]] <- rmsep
        val_data_range <- range(observed, na.rm = TRUE)
        res.rcv.svr$RMSEPper[[t]] <- (rmsep/(val_data_range[2]-val_data_range[1]))*100
      }
      res.rcv.svr$meas <- traitHSR$y
      save(val.ind,res.rcv.svr,file= paste0("results/unseen_seasons/SVR_model_",datatype,"_",agg,"_training",unique(traitHSR1$Season),"_testing",year2[y],"_",inVar,".RData"))
    }
  }
}
