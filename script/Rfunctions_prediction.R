get_XYdata <- function(species, year, excelfile, agg, fullrange) {
  filen <- paste("data/traits/", year, "_", excelfile, "_", species, ".csv", sep = "")
  ydata <- read.csv(filen)

  if (excelfile == 'CN_isotopes' && species == 'maize') {
    ydata <- ydata[complete.cases(ydata), ]
    ydata$CN <- as.numeric(ydata$CN)
  }

  # temp=ydata[,3:ncol(ydata)]
  # check <- colSums(temp)
  # ydata <- ydata[, check != 0]

  nvar <- ncol(ydata) - 2

  filen <- paste("data/reflectance/", year, "_hyperspectral_reflectance_", species, ".csv", sep = "")
  Xdata <- read.csv(filen)

  colnames(ydata)[1:2] <- colnames(Xdata)[1:2]

  Xdata$Plot <- as.character(Xdata$Plot)
  ydata$Plot <- as.character(ydata$Plot)

  if (species != "maize" && excelfile != "fluorescence") {
    Xdata$Repeat <- as.character(Xdata$Repeat)
    ydata$Repeat <- as.character(ydata$Repeat)
  }

  if (!fullrange) {
    Xdata <- cbind(Xdata[, 1:2], Xdata[, 152:(ncol(Xdata) - 100)])
  }

  data <- merge(ydata, Xdata, by = c("Plot", "Repeat"))
  # data <- data[, -2]

  if (excelfile == "ACi" && species == "barley") {
    check <- data$Jmax < 1000
    data <- data[check, ]
  }

  if (agg) {
    data <- data %>%
      group_by(Plot) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }

  data <- na.omit(data)

  if (agg==FALSE){
    data[,1]=paste0(data[,1],"_",data[,2])
  }
  if (species=="maize"){
    data <- data[, -2]
  }else if(species=="barley" && excelfile=="fluorescence"){
    data <- data[, -2]
  }


  return(list(data,nvar))
}


create_6folds <- function(species, year,scena){
  filen=paste0("data/splits/accgroups_",species,year,"rep_1.csv")
  acc_ref=read.csv(filen)
  folds=list()
  if (scena==1){
    f=1
    for (p in 2:3){
      for (r in 1:3){
        folds[[f]]=paste0(acc_ref[,p],"_",r)
        f=f+1
      }
    }
  }else if (scena==2){
    # folds0 <- createFolds(acc_ref$Accession, k = 6)
    for (f in 1:6){
      # test_ind=folds0[[f]]
      # temp=c(acc_ref[test_ind,2],acc_ref[test_ind,3])
      # temp=rep(temp,each=3)
      # repeats=rep(1:3,length(temp)/3)
      # folds[[f]]=paste0(temp,"_",repeats)
      # write.csv(folds[[f]],file=paste0("data/splits/",species,"/",year,"/scenario2_6folds_testfold",f,".csv"))

      ## Create inner test fold
      # train_ind=setdiff(1:nrow(acc_ref),test_ind)
      # train_acc=acc_ref[train_ind,]
      # innerfolds0 <- createFolds(train_acc$Accession, k = 5)
      # innerfolds=list()
      for (i in 1:5){
        # itest_ind=innerfolds0[[i]]
        # temp=c(train_acc[itest_ind,2],train_acc[itest_ind,3])
        # temp=rep(temp,each=3)
        # repeats=rep(1:3,length(temp)/3)
        # innerfolds[[i]]=paste0(temp,"_",repeats)
        # write.csv(innerfolds[[i]],file=paste0("data/splits/",species,"/",year,"/scenario2_otestfold",f,"_itestfold",i,".csv"))
      }
      temp=read.csv(file=paste0("data/splits/",species,"/",year,"/scenario2_6folds_testfold",f,".csv"))
      folds[[f]]=temp$x
    }
  }
  return(folds)
}

get_test_train_data_5fold <- function(data,varname, nvar,folds,f) {
  # trait=data[,1+i]
  # lower <- mean(trait) - 2*sd(trait)
  # upper <- mean(trait) + 2*sd(trait)
  # removeind=which(trait>=upper | trait<=lower)
  # if (length(removeind)!=0){
  #   data=data[-removeind,]
  # }
  ind_test=which(data$Plot %in% folds[[f]])
  ind_train=setdiff(1:nrow(data),ind_test)
  data_test=data[ind_test,]
  data_train=data[ind_train,]
  cat("fold: ",f, "Testing number: ", nrow(data_test),"\n")
  yind=which(colnames(data_train)%in%varname)

  data_train_i=cbind(data_train[,yind],data_train[,(2+nvar):ncol(data_train)])
  colnames(data_train_i)=c("y",colnames(data_train)[(2+nvar):ncol(data_train)])
  data_test_i=cbind(data_test[,yind],data_test[,(2+nvar):ncol(data_test)])
  colnames(data_test_i)=c("y",colnames(data_test)[(2+nvar):ncol(data_test)])
  train_plot=data$Plot[ind_train]
  test_plot=data$Plot[ind_test]

  return(list(data_train_i, data_test_i,train_plot,test_plot))
}

get_test_train_data_rep2fold <- function(data, scena, species, year, nvar, rep, p,i) {
  accgroups <- read.csv(paste("data/splits_2fold/accgroups_", species, year, "rep_", rep, ".csv", sep = ""), header = TRUE)
  accgroups <- data.frame(lapply(accgroups, as.character), stringsAsFactors = FALSE)

  if (scena == 1) {
    if (p==0){
      train_c="group1"
      test_c="group2"
    }else{
      train_c="group2"
      test_c="group1"
    }
    temp <- data.frame(Plot = accgroups[, train_c])
    data_train <- merge(data, temp, by = "Plot")

    temp <- data.frame(Plot = accgroups[, test_c])
    data_test <- merge(data, temp, by = "Plot")
  } else if (scena == 2) {
    if (p==0){
      train_c=1
      test_c=2
    }else{
      train_c=2
      test_c=1
    }
    kfolds <- read.csv(paste("data/splits_2fold/kfolds_", species,year, "_rep_", rep, ".csv", sep = ""), header = TRUE)
    kfolds=kfolds+1 # from python, the index start with 0
    temp0 <- data.frame(Plot = accgroups[kfolds[, train_c], 2:3])
    n <- nrow(temp0)
    temp <- data.frame(Plot = rep(NA, n * 2))
    temp[1:n, 1] <- temp0[, 1]
    temp[(n + 1):(n * 2), 1] <- temp0[, 2]
    data_train <- merge(data, temp, by = "Plot")

    temp0 <- data.frame(Plot = accgroups[kfolds[, test_c], 2:3])
    n <- nrow(temp0)
    temp <- data.frame(Plot = rep(NA, n * 2))
    temp[1:n, 1] <- temp0[, 1]
    temp[(n + 1):(n * 2), 1] <- temp0[, 2]
    data_test <- merge(data, temp, by = "Plot")
  }
  data_train_i=cbind(data_train[,1+i],data_train[,(2+nvar):ncol(data_train)])
  colnames(data_train_i)=c("y",colnames(data_train)[(2+nvar):ncol(data_train)])
  data_test_i=cbind(data_test[,1+i],data_test[,(2+nvar):ncol(data_test)])
  colnames(data_test_i)=c("y",colnames(data_test)[(2+nvar):ncol(data_test)])

  train_plot=data_train$Plot
  test_plot=data_test$Plot
  return(list(data_train_i, data_test_i,train_plot,test_plot))
}


get_model_tuning <- function(species,year,method, train_data,train_plot,scena, folds,f,tuning_cv,parallel) {
  innerfolds_train=list()
  innerfolds_test=list()
  if (scena!=3){
    for (i in 1:5){
      if (scena==1){
        train_folds=folds[setdiff(1:length(folds),f)]
        test_fold_i=train_folds[[i]]
        test_ind=which(train_plot %in% test_fold_i)
      }else if (scena==2){
        temp=read.csv(paste0("data/splits/",species,"/",year,"/scenario2_otestfold",f,"_itestfold",i,".csv"))
        test_ind=which(train_plot %in% temp$x)
      }
      innerfolds_test[[i]]=test_ind
      innerfolds_train[[i]]=setdiff(1:length(train_plot),test_ind)
    }

    if (tuning_cv=="LOO"){
      my_control <- trainControl(method="LOOCV")
    }else if (tuning_cv=="5fold"){
      my_control <- trainControl(method="cv",number=5, index=innerfolds_train,indexOut=innerfolds_test)
    }else{
      my_control <- trainControl(method="cv",number=1)
    }
  }else{
    my_control <- trainControl(method="cv",number=5)
  }

  if (parallel==TRUE){
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
  }

  if (method=="svm"){
    tuned_model <- train(y ~ .,data = train_data,method="svmLinear", preProcess=NULL,trControl=my_control,
                         tuneGrid = expand.grid(C = seq(0.25, 1, length = 4))) #if tuneGrid not specified, only C=1 is checked
  }else if (method=="PLSR"){
    tuned_model <- train(y ~ .,data = train_data,method="pls", preProcess=NULL,trControl=my_control,tuneGrid = expand.grid(ncomp = 1:25))
  }else if (method=="svmPoly"){
    tuned_model <- train(y ~ .,data = train_data,method="svmPoly", preProcess=NULL,trControl=my_control) # degree,scale,C are run automatically
  }else{
    tuned_model <- train(y ~ .,data = train_data,method=method, preProcess=NULL,trControl=my_control)
  }
  if (parallel==TRUE){
    stopCluster(cl)
  }
  return(tuned_model)
}


get_model_tuning_rep2fold <- function(method, species, scena, year, train_data,train_plot, nvar, rep, p, i,tuning_cv,parallel) {

  if (scena == 1) {
    my_control <- trainControl(method="cv",number=5)
    if (method=="svm"){
      tuned_model <- train(y ~ .,data = train_data,method="svmLinear", preProcess=NULL,trControl=my_control,
                           tuneGrid = expand.grid(C = seq(0.25, 1, length = 4))) #if tuneGrid not specified, only C=1 is checked
    }
  } else if (scena == 2) {
    innerfolds_train=list()
    innerfolds_test=list()
    for (irep in 0:4) {
      filen <- paste("data/splits_2fold/tuning/train_ifolds_", species, year, "_rep_", rep, "_ifolds_", irep, "_iter_", p, ".csv", sep = "")
      ifold_train <- read.csv(filen, header = TRUE)
      ifold_train <- c(ifold_train[,1], ifold_train[,2])

      filen <- paste("data/splits_2fold/tuning/test_ifolds_", species, year, "_rep_", rep, "_ifolds_", irep, "_iter_", p, ".csv", sep = "")
      ifold_test <- read.csv(filen, header = TRUE)
      ifold_test <- c(ifold_test[,1], ifold_test[,2])

      innerfolds_train[[irep+1]]=which(train_plot%in%ifold_train)
      innerfolds_test[[irep+1]]=which(train_plot%in%ifold_test)
    }
  }

  if (parallel==TRUE){
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
  }

  if (tuning_cv=="LOO"){
    my_control <- trainControl(method="LOOCV")
  }else if (tuning_cv=="5fold" && scena==2){
    my_control <- trainControl(method="cv",number=5, index=innerfolds_train,indexOut=innerfolds_test)
  }else if (tuning_cv=="5fold" && scena==1){
    my_control <- trainControl(method="cv",number=5)
  }

  if (method=="svm"){
    tuned_model <- train(y ~ .,data = train_data,method="svmLinear", preProcess=NULL,trControl=my_control,
                         tuneGrid = expand.grid(C = seq(0.25, 1, length = 4))) #if tuneGrid not specified, only C=1 is checked
  }else if (method=="PLSR"){
    tuned_model <- train(y ~ .,data = train_data,method="pls", preProcess=NULL,trControl=my_control,tuneGrid = expand.grid(ncomp = 1:25))
  }else if (method=="svmPoly"){
    tuned_model <- train(y ~ .,data = train_data,method="svmPoly", preProcess=NULL,trControl=my_control) # degree,scale,C are run automatically
  }else{
    tuned_model <- train(y ~ .,data = train_data,method=method, preProcess=NULL,trControl=my_control)
  }
  if (parallel==TRUE){
    stopCluster(cl)
  }
  return(tuned_model)
}

get_R2 <- function(meas,pred){
  R2=1-sum((meas-pred)^2)/sum((meas-mean(meas))^2)
  return(R2)
}

predict_scenario1and2_6fold <- function(species, year, scena, method, excelfile, tuning,parallel,spec_ind) {
  agg <- FALSE
  fullrange <- TRUE

  cat("LOADING DATA \n")
  # Data preparation
  res <- get_XYdata(species, year, excelfile, agg, fullrange)
  data=res[[1]]
  if (sum("Repeat" %in% colnames(data))>0){
    ind<-which(colnames(data)=="Repeat")
    data<-data[,-ind]
  }

  nvar0=res[[2]]
  varnames=colnames(data)[2:(1+nvar0)]
  if (length(spec_ind)!=0){
    spec_ind<-match(spec_ind,varnames)
    varnames<-varnames[spec_ind]
    nvar<-length(varnames)
  }

  if (species=="maize" && excelfile=="ACi"){
    removeind<-which(data$Plot%in%"2122_1")
    data<-data[-removeind,]
  }
  # Get outer CV
  kfolds<-create_6folds(species,year,scena)

  if (method=="PLSR"|method=="svm" | method=="rf"){
    ntuned=1
  }else if (method=="glmnet"){
    ntuned=2
  }else if (method=="svmPoly"){
    ntuned=3
  }


  if (tuning!="none") {
    tuned_param=array(0,dim=c(length(kfolds),nvar,ntuned))
    colnames(tuned_param) <- varnames
  }

  for (i in 1:nvar) {
    varname=varnames[i]
    cat("Folder: ", excelfile, ", Variable: ", varname, ", Year: ", year, "\n")

    for (f in 1:6){
      result <- get_test_train_data_5fold(data,varname,nvar0,kfolds,f)
      data_train_i <- result[[1]]
      data_test_i <- result[[2]]
      train_plot <- result[[3]]
      test_plot <- result[[4]]
      ytrain <- data.frame(meas = data_train_i$y)
      ytest <- data.frame(meas = data_test_i$y)

      if ((species=="maize" && varname=="N") | (species=="maize" &&varname=="CN")){
        ind=(1500-350+2):(ncol(data_test_i)-100)
        data_test_i=cbind(data_test_i[,1],data_test_i[,ind])
        data_train_i=cbind(data_train_i[,1],data_train_i[,ind])
        colnames(data_test_i)[1]="y"
        colnames(data_train_i)[1]="y"
      }
      start_time <- Sys.time()
      model=get_model_tuning(species,year,method,data_train_i,train_plot,scena,kfolds,f,tuning,parallel)
      end_time <- Sys.time()
      end_time - start_time

      predtrain=predict(model, newdata =data_train_i)
      predtest=predict(model, newdata =data_test_i)
      tuned_param[f, i,] <- unlist(model$bestTune)

      restrain <- cbind(ytrain, predtrain)
      colnames(restrain) <- c("meas", "pred")
      restest <- cbind(ytest, predtest)
      colnames(restest) <- c("meas", "pred")
      cat("Training R-squared: ", round(cor(restrain$meas, restrain$pred)^2, 2), ", Correlation: ", round(cor(restrain$meas, restrain$pred), 2), "\n")
      cat("Testing R-squared: ", round(cor(restest$meas, restest$pred)^2, 2), ", Correlation: ", round(cor(restest$meas, restest$pred), 2), "\n")

      testing_filen <- paste("results/scenario", scena, "/", "testing_", species, "_", method,"_", tuning, "_", varname, year, "_testfold", f, ".txt", sep = "")
      write.csv(restest, file = testing_filen, row.names = FALSE)

      training_filen <- paste("results/scenario", scena, "/", "training_", species, "_", method,"_", tuning, "_", varname, year,  "_testfold", f, ".txt", sep = "")
      write.csv(restrain, file = training_filen, row.names = FALSE)
    }
  }

  if (tuning!="none") {
    # tuned_param_filen <- paste("results/scenario", scena, "/", "tuned_param", species, "_", excelfile, year,"_", method, "_", tuning, ".txt", sep = "")
    # write.csv(tuned_param, file = tuned_param_filen, row.names = FALSE)
    tuned_param_filen <- paste("results/scenario", scena, "/", "tuned_param", species, "_", excelfile, year,"_", method, "_", tuning, ".RData", sep = "")
    save(tuned_param, file = tuned_param_filen)
  }
}



predict_scenario1and2_rep2fold <- function(species, year, scena, method, excelfile, tuning,parallel) {
  agg <- TRUE
  fullrange <- TRUE

  # Data preparation
  res <- get_XYdata(species, year, excelfile, agg, fullrange)
  data=res[[1]]
  nvar=res[[2]]
  varnames=colnames(data)[2:(1+nvar)]

  if (species=="maize" && excelfile=="ACi"){
    removeind=which(data$Plot%in%"2122")
    data=data[-removeind,]
  }

  if (method=="PLSR"|method=="svm"){
    ntuned=1
  }else if (method=="glmnet"){
    ntuned=2
  }else if (method=="svmPoly"){
    ntuned=3
  }


  if (tuning!="none") {
    tuned_param=array(0,dim=c(10,nvar,ntuned))
    colnames(tuned_param) <- varnames
  }

  for (i in 1:nvar) {
    varname=colnames(data)[1 + i]
    cat("Folder: ", excelfile, ", Variable: ", varname, ", Year: ", year, "\n")
    count=1
    for (rep in 1:5) {
      for (p in 0:1) {
        # Obtain training and testing data
        result <- get_test_train_data_rep2fold(data, scena, species, year, nvar, rep, p,i)
        data_train_i <- result[[1]]
        data_test_i <- result[[2]]
        train_plot<- result[[3]]
        test_plot<- result[[4]]
        ytrain <- data.frame(meas = data_train_i$y)
        ytest <- data.frame(meas = data_test_i$y)

        if (varname=="N" | varname=="CN"){
          ind=(1500-350+2):(ncol(data_test_i)-100)
          data_test_i=cbind(data_test_i[,1],data_test_i[,ind])
          data_train_i=cbind(data_train_i[,1],data_train_i[,ind])
          colnames(data_test_i)[1]="y"
          colnames(data_train_i)[1]="y"
        }

        model=get_model_tuning_rep2fold(method, species, scena, year, data_train_i,train_plot, nvar, rep, p, i,tuning,parallel)
        predtrain=predict(model, newdata =data_train_i)
        predtest=predict(model, newdata =data_test_i)
        restrain <- cbind(ytrain, predtrain)
        colnames(restrain) <- c("meas", "pred")
        restest <- cbind(ytest, predtest)
        colnames(restest) <- c("meas", "pred")
        cat("Training R-squared: ", round(cor(restrain$meas, restrain$pred)^2, 2), ", Correlation: ", round(cor(restrain$meas, restrain$pred), 2), "\n")
        cat("Testing R-squared: ", round(cor(restest$meas, restest$pred)^2, 2), ", Correlation: ", round(cor(restest$meas, restest$pred), 2), "\n")


        testing_filen <- paste("results/scenario", scena, "_2fold/", "testing_", species, "_", method, "_", tuning, "_", colnames(data)[1 + i], year, "_rep", rep, "_iter", p, ".txt", sep = "")
        write.csv(restest, file = testing_filen, row.names = FALSE)

        training_filen <- paste("results/scenario", scena, "_2fold/", "training_", species, "_", method, "_", tuning, "_", colnames(data)[1 + i], year, "_rep", rep, "_iter", p, ".txt", sep = "")
        write.csv(restrain, file = training_filen, row.names = FALSE)

        tuned_param[count, i,] <- unlist(model$bestTune)
        count=count+1
      }
    }
  }

  if (tuning!="none") {
    tuned_param_filen <- paste("results/scenario", scena, "_2fold/", "tuned_param", species, "_", excelfile, year,"_", method, "_", tuning, ".RData", sep = "")
    save(tuned_param, file = tuned_param_filen)
  }
}

predict_scenario3 <- function(species, method, excelfile, tuning,agg, parallel) {
  scena=3
  fullrange <- TRUE

  cat("LOADING DATA \n")
  # Data preparation
  res <- get_XYdata(species, 2021, excelfile, agg, fullrange)
  data=res[[1]]
  nvar=res[[2]]
  varnames=colnames(data)[2:(1+nvar)]
  data21=data

  res <- get_XYdata(species, 2022, excelfile, agg, fullrange)
  data=res[[1]]
  data22=data


  if (method=="PLSR"|method=="svm"){
    ntuned=1
  }else if (method=="glmnet"){
    ntuned=2
  }else if (method=="svmPoly"){
    ntuned=3
  }


  if (tuning!="none") {
    tuned_param=array(0,dim=c(2,nvar,ntuned))
  }

  for (i in 1:nvar) {
    varname=varnames[i]
    cat("Folder: ", excelfile, ", Variable: ", varname,  "\n")
    count=1
    # use the whole dataset
    for (testyear in c(2021,2022)){
      if (testyear==2021){
        data_train=data22
        data_test=data21
      }else{
        data_train=data21
        data_test=data22
      }
      data_train_i=cbind(data_train[,1+i],data_train[,(2+nvar):ncol(data_train)])
      colnames(data_train_i)=c("y",colnames(data_train)[(2+nvar):ncol(data_train)])
      data_test_i=cbind(data_test[,1+i],data_test[,(2+nvar):ncol(data_test)])
      colnames(data_test_i)=c("y",colnames(data_test)[(2+nvar):ncol(data_test)])
      ytrain <- data.frame(meas = data_train_i$y)
      ytest <- data.frame(meas = data_test_i$y)

      model=get_model_tuning(species,method,data_train_i,c(),scena,kfolds,f,tuning,parallel)
      tuned_param[count, i,] <- unlist(model$bestTune)
      count=count+1
      predtrain=predict(model, newdata =data_train_i)
      predtest=predict(model, newdata =data_test_i)
      restrain <- cbind(ytrain, predtrain)
      colnames(restrain) <- c("meas", "pred")
      restest <- cbind(ytest, predtest)
      colnames(restest) <- c("meas", "pred")
      cat("Training R-squared: ", round(cor(restrain$meas, restrain$pred)^2, 2), ", Correlation: ", round(cor(restrain$meas, restrain$pred), 2), "\n")
      cat("Testing R-squared: ", round(cor(restest$meas, restest$pred)^2, 2), ", Correlation: ", round(cor(restest$meas, restest$pred), 2), "\n")

      if (agg){
        add="_agg"
      }else{
        add=""
      }
      testing_filen <- paste("results/scenario", scena,add ,"/", "whole_data_testing_", species, "_", method,"_", tuning, "_", varname,"_testingyear",testyear, ".txt", sep = "")
      write.csv(restest, file = testing_filen, row.names = FALSE)

      training_filen <- paste("results/scenario", scena, add, "/", "whole_data_training_", species, "_", method,"_", tuning, "_", varname,"_testingyear",testyear, ".txt", sep = "")
      write.csv(restrain, file = training_filen, row.names = FALSE)
    }
  }

  if (tuning!="none") {
    tuned_param_filen <- paste("results/scenario", scena,add , "/", "whole_data_tuned_param", species, "_", excelfile,"_", method, "_", tuning, ".RData", sep = "")
    save(tuned_param, file = tuned_param_filen)
  }
}


