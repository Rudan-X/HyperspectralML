## data dimensions##
# CN_isotopes: original_2021:599, 2022:1869
# ACi: original 2021: 382, 2022: 441
# light sat: original 2021: 831, 2022: 1681


load_traits <- function(species, year, excelfile, agg_type,dataonly){
  filen <- paste("data/traits/", year, "_", excelfile, "_", species, ".csv", sep = "")
  ydata <- read.csv(filen)
  
  if (excelfile == 'CN_isotopes' && species == 'maize') {
    ydata <- ydata[complete.cases(ydata), ]
    ydata$CN <- as.numeric(ydata$CN)
    ind<-which(is.na(ydata),arr.ind = T)
    if (dim(ind)[1]!=0){
      ydata<-ydata[-ind[,2],]
    }
  }
  colnames(ydata)[c(1,2)]<-c("Plot","Repeat")
  ydata$Plot <- as.character(ydata$Plot)
  if (species != "maize" && excelfile != "fluorescence") {
    ydata$Repeat <- as.character(ydata$Repeat)
  }
  if (agg_type=="plot") {
    ydata <- ydata %>%
      group_by(Plot) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }else if (agg_type=="acc"){
    acc_ref <- read.csv(paste0("data/accessions/",year,"_",species,"_reference.csv"))
    ind<-match(ydata$Plot,acc_ref$Plot)
    ydata$accession<-acc_ref$Accession[ind]
    ydata <- ydata %>%
      group_by(accession) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }
  if (dataonly){
    ind=which(colnames(ydata)%in%c("Plot","Repeat","accession"))
    return(ydata[,-ind])
  }else{return(ydata)}
}

 
get_XYdata <- function(species, year, excelfile, agg_type) {
  # Load traits
  ydata<-load_traits(species,year,excelfile,"none",FALSE)
  
  # Load reflectance
  filen <- paste("data/reflectance/", year, "_hyperspectral_reflectance_", species, ".csv", sep = "")
  Xdata <- read.csv(filen)
  ind<-which(is.na(Xdata),arr.ind = T)
  if (dim(ind)[1]!=0){
    Xdata<-Xdata[-ind[,2],]
  }
  
  
  # Some traits like CN_isotopes come already with averaged plots
  pooledplot<-as.numeric(ydata$Plot[ydata$Repeat%in%"Pooled"])
  if (length(pooledplot>0)){
    Xdata$Repeat[Xdata$Plot%in%pooledplot]="Pooled"
  }
  
  data <- merge(ydata, Xdata, by = c("Plot", "Repeat"))
  # data <- data[, -2]
  
  if (excelfile == "ACi" && species == "barley") {
    check <- data$Jmax < 1000
    data <- data[check, ]
  }
  # 
  # if (species=="maize" && excelfile=="ACi"){
  #   removeind=which(data$Plot%in%"2122")
  #   data=data[-removeind,]
  # }
  
  acc_ref <- read.csv(paste0("data/accessions/",year,"_",species,"_reference.csv"))
  if (year==2021){
    temp<-strsplit(acc_ref$Accession,split='_')
    numb<- unlist(temp)[2*(1:length(temp))  ]
    accessions<-paste0("SSA_00",numb)
  }else{accessions<-acc_ref$Accession}
  
  
  if (agg_type=="plot") {
    data <- data %>%
      group_by(Plot) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
    ind<-match(data$Plot,acc_ref$Plot)
    
    if (sum(colnames(data)%in%"Repeat")==0){
      temp=cbind(data$Plot,accessions[ind])
      data<-cbind(temp,data[,2:ncol(data)])
      colnames(data)[1:2]<-c("Plot","Accession")
    }else{
      data$Repeat<-accessions[ind]
      colnames(data)[colnames(data)=="Repeat"]="Accession"
    }
    
  }else if (agg_type=="acc"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Accession<-accessions[ind]
    data <- data %>%
      group_by(Accession) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }else if (agg_type=="none"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Repeat<-accessions[ind]
    colnames(data)[colnames(data)=="Repeat"]="Accession"
  }
  data<-data[complete.cases(data),]
  return(data)
}

get_both_proxies <- function(species, year, excelfile, agg_type,add_pooled) {
  # Load traits
  ydata<-load_traits(species,year,excelfile,"none",FALSE)
  
  Ndata<-load_traits(species,year,"CN_isotopes","none",FALSE)
  pooledplot<-as.numeric(Ndata$Plot[Ndata$Repeat%in%"Pooled"])
  
  # Load reflectance
  if (exists("Xdatatot")){
    rm(Xdatatot)
  }
  
  for (filen in c("SLA","CN_isotopes","fluorescence")){
    Xdata<-load_traits(species,year,filen,"none",FALSE)
    
    ind<-which(is.na(Xdata),arr.ind = T)
    if (dim(ind)[1]!=0){
      Xdata<-Xdata[-ind[,2],]
    }
    
    if (add_pooled){
      Xdata$Repeat[Xdata$Plot%in%pooledplot]="Pooled"
    }
    
    if (!exists("Xdatatot")){
      Xdatatot<-Xdata
    }else{
      Xdatatot<-merge(Xdatatot,Xdata)
    }
  }
  
  # ind<-match(c("C","CN","d13C"),colnames(Xdatatot))
  # Xdatatot<-Xdatatot[,-ind]
    
  data <- merge(ydata, Xdatatot, by = c("Plot", "Repeat"))
  
  
  acc_ref <- read.csv(paste0("data/accessions/",year,"_",species,"_reference.csv"))
  if (year==2021){
    temp<-strsplit(acc_ref$Accession,split='_')
    numb<- unlist(temp)[2*(1:length(temp))  ]
    accessions<-paste0("SSA_00",numb)
  }else{accessions<-acc_ref$Accession}
  
  
  if (agg_type=="plot") {
    data <- data %>%
      group_by(Plot) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
    ind<-match(data$Plot,acc_ref$Plot)
    
    if (sum(colnames(data)%in%"Repeat")==0){
      temp=cbind(data$Plot,accessions[ind])
      data<-cbind(temp,data[,2:ncol(data)])
      colnames(data)[1:2]<-c("Plot","Accession")
    }else{
      data$Repeat<-accessions[ind]
      colnames(data)[colnames(data)=="Repeat"]="Accession"
    }
    
  }else if (agg_type=="acc"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Accession<-accessions[ind]
    data <- data %>%
      group_by(Accession) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }else if (agg_type=="none"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Repeat<-accessions[ind]
    colnames(data)[colnames(data)=="Repeat"]="Accession"
  }
  return(data)
}

get_organic_proxies <- function(species, year, excelfile, agg_type,add_pooled) {
  # Load traits
  ydata<-load_traits(species,year,excelfile,"none",FALSE)
  
  Ndata<-load_traits(species,year,"CN_isotopes","none",FALSE)
  pooledplot<-as.numeric(Ndata$Plot[Ndata$Repeat%in%"Pooled"])
  
  # Load reflectance
  if (exists("Xdatatot")){
    rm(Xdatatot)
  }
  
  for (filen in c("SLA","CN_isotopes")){
    Xdata<-load_traits(species,year,filen,"none",FALSE)
    
    ind<-which(is.na(Xdata),arr.ind = T)
    if (dim(ind)[1]!=0){
      Xdata<-Xdata[-ind[,2],]
    }
    
    if (add_pooled){
      Xdata$Repeat[Xdata$Plot%in%pooledplot]="Pooled"
    }
    
    if (!exists("Xdatatot")){
      Xdatatot<-Xdata
    }else{
      Xdatatot<-merge(Xdatatot,Xdata)
    }
  }
  
  data <- merge(ydata, Xdatatot, by = c("Plot", "Repeat"))
  
  
  acc_ref <- read.csv(paste0("data/accessions/",year,"_",species,"_reference.csv"))
  if (year==2021){
    temp<-strsplit(acc_ref$Accession,split='_')
    numb<- unlist(temp)[2*(1:length(temp))  ]
    accessions<-paste0("SSA_00",numb)
  }else{accessions<-acc_ref$Accession}
  
  
  if (agg_type=="plot") {
    data <- data %>%
      group_by(Plot) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
    ind<-match(data$Plot,acc_ref$Plot)
    
    if (sum(colnames(data)%in%"Repeat")==0){
      temp=cbind(data$Plot,accessions[ind])
      data<-cbind(temp,data[,2:ncol(data)])
      colnames(data)[1:2]<-c("Plot","Accession")
    }else{
      data$Repeat<-accessions[ind]
      colnames(data)[colnames(data)=="Repeat"]="Accession"
    }
    
  }else if (agg_type=="acc"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Accession<-accessions[ind]
    data <- data %>%
      group_by(Accession) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }else if (agg_type=="none"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Repeat<-accessions[ind]
    colnames(data)[colnames(data)=="Repeat"]="Accession"
  }
  return(data)
}

get_fluorescence_proxies <- function(species, year, excelfile, agg_type,add_pooled) {
  # Load traits
  ydata<-load_traits(species,year,excelfile,"none",FALSE)
  
  Ndata<-load_traits(species,year,"CN_isotopes","none",FALSE)
  pooledplot<-as.numeric(Ndata$Plot[Ndata$Repeat%in%"Pooled"])
  
  # Load reflectance
  if (exists("Xdatatot")){
    rm(Xdatatot)
  }
  
  for (filen in c("fluorescence")){
    Xdata<-load_traits(species,year,filen,"none",FALSE)
    
    ind<-which(is.na(Xdata),arr.ind = T)
    if (dim(ind)[1]!=0){
      Xdata<-Xdata[-ind[,2],]
    }
    
    if (add_pooled){
      Xdata$Repeat[Xdata$Plot%in%pooledplot]="Pooled"
    }
    
    if (!exists("Xdatatot")){
      Xdatatot<-Xdata
    }else{
      Xdatatot<-merge(Xdatatot,Xdata)
    }
  }
  
  data <- merge(ydata, Xdatatot, by = c("Plot", "Repeat"))
  
  
  acc_ref <- read.csv(paste0("data/accessions/",year,"_",species,"_reference.csv"))
  if (year==2021){
    temp<-strsplit(acc_ref$Accession,split='_')
    numb<- unlist(temp)[2*(1:length(temp))  ]
    accessions<-paste0("SSA_00",numb)
  }else{accessions<-acc_ref$Accession}
  
  
  if (agg_type=="plot") {
    data <- data %>%
      group_by(Plot) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
    ind<-match(data$Plot,acc_ref$Plot)
    
    if (sum(colnames(data)%in%"Repeat")==0){
      temp=cbind(data$Plot,accessions[ind])
      data<-cbind(temp,data[,2:ncol(data)])
      colnames(data)[1:2]<-c("Plot","Accession")
    }else{
      data$Repeat<-accessions[ind]
      colnames(data)[colnames(data)=="Repeat"]="Accession"
    }
    
  }else if (agg_type=="acc"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Accession<-accessions[ind]
    data <- data %>%
      group_by(Accession) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  }else if (agg_type=="none"){
    ind<-match(data$Plot,acc_ref$Plot)
    data$Repeat<-accessions[ind]
    colnames(data)[colnames(data)=="Repeat"]="Accession"
  }
  return(data)
}

get_R2 <- function(meas,pred){
  R2=1-sum((meas-pred)^2)/sum((meas-mean(meas))^2)
  return(R2)
}

mycreateFolds <- function(strat_id, k) {
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }	
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))	
  return(folds)
}

stat_box_data <- function(y) {
  return( 
    data.frame(
      y = max(y),  #may need to modify this depending on your data
      label = paste(round(median(y), 2), '\n')
    )
  )
}


samp_seen_lines_inner <- function(data,acctot,outfold,infold,inner_flag){
  stopsamp<-FALSE
  restart<-TRUE
  while (restart){
    while (stopsamp==FALSE){
      # randomly sample 20% testing data which belong to accessions already 
      # included in the training set
      test_ind <- sample(1:nrow(data),size=nrow(data)/outfold)
      train_ind<-setdiff(1:nrow(data),test_ind)
      acctot_test<-acctot[test_ind]
      acctot_train<-acctot[train_ind]
      
      if (all(acctot_test%in%acctot_train)){
        stopsamp<-TRUE
      }
    }
    data_train<-data[train_ind,]
    data_test<-data[test_ind,]
    
    if (inner_flag){
      ## inner CV fold
      inner_folds_test <- list()
      inner_folds_train <- list()
      
      goodsamp<-TRUE
      i<-0
      while (goodsamp & i<5){
        stopsamp0<-FALSE
        start_time <- Sys.time()
        while (stopsamp0==FALSE){
          test_ind_inner <- sample(1:length(train_ind),size=length(train_ind)/infold)
          train_ind_inner <- setdiff(1:length(train_ind),test_ind_inner)
          acctot_test_inner<-acctot_train[test_ind_inner]
          acctot_train_inner<-acctot_train[train_ind_inner]
          end_time <- Sys.time()
          time<-difftime(end_time,start_time,units = "mins")[[1]]
          if (all(acctot_test_inner%in%acctot_train_inner) | time>2){
            stopsamp0<-TRUE
          }
        }
        
        if (all(acctot_test_inner%in%acctot_train_inner) | time<2){ # possible to sample inner testing fold within 2 minutes
          i<-i+1
          inner_folds_test[[i]]<-test_ind_inner
          inner_folds_train[[i]]<-train_ind_inner
          goodsamp=TRUE
        }else{
          print("Not complete inner fold sampling. Outer sampling restarted")
          goodsamp=FALSE
          next
        }
      }
      if (goodsamp){
        restart<-FALSE
      }
    }else{
      restart<-FALSE
      inner_folds_test <- c()
      inner_folds_train <- c()
    }
  }
  return(list(data_train,data_test,train_ind,test_ind,inner_folds_test,inner_folds_train))
}

mycreateInnerFolds <- function(data_Train,data_Test, k) {
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }	
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))	
  return(folds)
}

prediction_scenario1 <- function(method,species,years,excelfile,variables, agg_type,tuned){
  for (year in years){
    data<-get_XYdata(species, year, excelfile,agg_type) 
    ind=which(colnames(data)=="Repeat")
    if (length(ind)>0){
      data<-data[,-ind]
    }
    if (excelfile=="lightSaturatedGasExchange" & all(years==c(2021,2022,2023))){
      ind<-which(colnames(data)%in%c("gs_sat","iWUE_sat"))
      data<-data[,-ind]
    }
    if (year==years[1]){
      datayears<-data
    }else{
      datayears<-rbind(datayears,data)
    }
  }
  
  data0<-datayears
  
  for (varname in variables){
    data<-data0
    accessions<-unique(data$Accession)
    acctot<-data$Accession
    
    yind<-which(colnames(data)%in%varname)
    
    if (varname=="N" | varname=="CN"){
      hstart<-1500
      hend<-2400
    }else{
      hstart<-350
      hend<-2500
    }
    
    hstart<-paste0("mean_",hstart)
    hstart<-which(colnames(data)%in%hstart)
    
    hend<-paste0("mean_",hend)
    hend<-which(colnames(data)%in%hend)
    data<-cbind(data[,yind],data[,hstart:hend])
    ind=which(is.na(data[,1]))
    if (length(ind)>0){data=data[-ind,]}
    
    colnames(data)[1]<-"y"
    nrep<-100
    training_corr2<-matrix(0,1,nrep)
    testing_corr2<-matrix(0,1,nrep)
    training_R2<-matrix(0,1,nrep)
    testing_R2<-matrix(0,1,nrep)
    test_result<-list()
    train_result<-list()
    
    
    if (agg_type=="acc"){
      ofold=20
      ifold=20
    }else if (agg_type=="plot" & excelfile=="fluorescence" & all(years==c(2022))){
      ofold=13
      ifold=13
    }else if (agg_type=="plot" & excelfile=="fluorescence" & all(years==c(2023))){
      ofold=13
      ifold=13
    }else if (agg_type=="plot" & excelfile=="CN_isotopes" & all(years==c(2023))){
      ofold=13
      ifold=13
    }else if (agg_type=="plot" & excelfile=="SLA" & all(years==c(2023))){
      ofold=15
      ifold=15
    }else if (agg_type=="plot" & excelfile=="lightSaturatedGasExchange" & all(years==c(2023))){
      ofold=13
      ifold=13
    }else{
      ofold=10
      ifold=10
    }
    for (r in 1:100){
      print(paste0("Current repetition: ",r))
      if (tuned){
        res<-samp_seen_lines_inner(data,acctot,ofold,ifold,TRUE)
        data_train<-res[[1]]
        data_test<-res[[2]]
        train_ind<-res[[3]]
        test_ind<-res[[4]]
        inner_folds_test<-res[[5]]
        inner_folds_train<-res[[6]]
        my_control <- trainControl(method="cv",number=5,index=inner_folds_train,indexOut=inner_folds_test)
        if (method=="PLSR"){   
          model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,
                         trControl=my_control,tuneGrid = expand.grid(ncomp = 1:50))
        }else if (method=="SVR"){
          model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,
                         trControl=my_control,tuneGrid = expand.grid(C = seq(0.25, 1, length = 4)))
        }
      }else{
        res<-samp_seen_lines_inner(data,acctot,ofold,ifold,FALSE)
        data_train<-res[[1]]
        data_test<-res[[2]]
        if (method=="PLSR"){   
          model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,trControl=trainControl(method="none"))
          # model<- plsr(y ~ ., ncomp = 2, data = data_train)
          
        }else if (method=="SVR"){
          model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
        }
        
      }

        
      pred_train<-predict(model,data_train)
      training_corr2[r]<-postResample(pred_train, data_train$y)[2]
      training_R2[r]<-get_R2(data_train$y,pred_train)
      pred_test<-predict(model,data_test)
      testing_corr2[r]<-postResample(pred_test, data_test$y)[2]
      testing_R2[r]<-get_R2(data_test$y,pred_test)
      train_result[[r]]=cbind(pred_train, data_train$y)
      test_result[[r]]=cbind(pred_test, data_test$y)
    }
    
    if (length(years)>1){
      yearsp=c()
      for (y in years){
        yearsp=paste0(yearsp,"_",y)
      }
    }else{
      yearsp=paste0("_",years)
    }
    
    if (tuned & method=="PLSR"){
      method_lab="PLSR"
    }else if(tuned & method=="SVR"){
      method_lab="SVRtuned"
    }else if(!tuned & method=="PLSR"){
      method_lab="PLSRnotuned"
    }else if(!tuned & method=="SVR"){
      method_lab="SVR"
    }
    filename=paste0("results/scenario1/",species,"_",method_lab,"_", varname,yearsp,"_",agg_type,".RData")
    
    save(train_result,test_result,testing_R2,testing_corr2,training_R2,training_corr2,file=filename)
    print(paste0("Variable: ",varname, yearsp, ",aggregated by: ", agg_type,", median testing R2: ", median(testing_R2)))
  }
}

prediction_scenario2 <- function(method,species,years,excelfile,variables, agg_type,tuned){
  for (year in years){
    data<-get_XYdata(species, year, excelfile,agg_type) 
    ind=which(colnames(data)=="Repeat")
    if (length(ind)>0){
      data<-data[,-ind]
    }
    if (year==years[1]){
      datayears<-data
    }else{
      datayears<-rbind(datayears,data)
    }
  }
  
  data0<-datayears
  
  for (varname in variables){
    data<-data0
    accessions<-unique(data$Accession)
    acctot<-data$Accession
    
    yind<-which(colnames(data)%in%varname)
    
    if (varname=="N" | varname=="CN"){
      hstart<-1500
      hend<-2400
    }else{
      hstart<-350
      hend<-2500
    }

    hstart<-paste0("mean_",hstart)
    hstart<-which(colnames(data)%in%hstart)
    
    hend<-paste0("mean_",hend)
    hend<-which(colnames(data)%in%hend)
    data<-cbind(data[,yind],data[,hstart:hend])
    colnames(data)[1]<-"y"
    nrep<-20
    training_corr2<-matrix(0,5,nrep)
    testing_corr2<-matrix(0,5,nrep)
    training_R2<-matrix(0,5,nrep)
    testing_R2<-matrix(0,5,nrep)
    test_result=list()
    train_result=list()
    for (r in 1:nrep){
      accfolds <- mycreateFolds(accessions, k = 5)
      print(paste0("Current repetition: ",r))
      train_result[[r]]=list()
      test_result[[r]]=list()
      for (f in 1:length(accfolds)){
        test_acc<-accessions[accfolds[[f]]]
        train_acc<-accessions[setdiff(1:length(accessions),accfolds[[f]])]
        test_ind<-which(acctot%in%test_acc)
        train_ind<-setdiff(1:nrow(data),test_ind)
        data_train<-data[train_ind,]
        data_test<-data[test_ind,]
        
        if (tuned){
          acctot_train<-acctot[train_ind]
          ## inner CV fold
          inner_accfolds <- mycreateFolds(train_acc, k = 5)
          inner_folds_train <- list()
          inner_folds_test  <-list()
          for (i in 1:length(inner_accfolds)){
            inner_folds_test[[i]]<-which(acctot_train%in%train_acc[inner_accfolds[[i]]])
            inner_folds_train[[i]]<-setdiff(1:length(acctot_train),inner_folds_test[[i]])
          }

          my_control <- trainControl(method="cv",number=5,index=inner_folds_train,indexOut=inner_folds_test)
          if (method=="PLSR"){
            model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,
                           trControl=my_control,tuneGrid = expand.grid(ncomp = 1:25))
          }else if (method=="SVR"){
            model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,
                           trControl=my_control,tuneGrid = expand.grid(C = seq(0.25, 1, length = 4)))
          }

        }else{
          if (method=="PLSR"){
            model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,trControl=trainControl(method="none"))
          }else if (method=="SVR"){
            model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
          }
        }

        pred_train<-predict(model,data_train)
        training_corr2[f,r]<-postResample(pred_train, data_train$y)[2]
        training_R2[f,r]<-get_R2(data_train$y,pred_train)
        pred_test<-predict(model,data_test)
        testing_corr2[f,r]<-postResample(pred_test, data_test$y)[2]
        testing_R2[f,r]<-get_R2(data_test$y,pred_test)
        train_result[[r]][[f]]=cbind(pred_train, data_train$y)
        test_result[[r]][[f]]=cbind(pred_test, data_test$y)
      }
    }
    if (length(years)>1){
      yearsp=c()
      for (y in years){
        yearsp=paste0(yearsp,"_",y)
      }
    }else{
      yearsp=paste0("_",years)
    }
    
    if (tuned & method=="PLSR"){
      method_lab="PLSR"
    }else if(tuned & method=="SVR"){
      method_lab="SVRtuned"
    }else if(!tuned & method=="PLSR"){
      method_lab="PLSRnotuned"
    }else if(!tuned & method=="SVR"){
      method_lab="SVR"
    }
    filename=paste0("results/scenario2/",species,"_",method_lab,"_", varname,yearsp,"_",agg_type,".RData")
    save(train_result,test_result,testing_R2,testing_corr2,training_R2,training_corr2,file=filename)
    print(paste0("Variable: ",varname, yearsp, ",aggregated by: ", agg_type,", median testing R2: ", median(testing_R2)))
  }
}


prediction_scenario2_permutation <- function(method,species,years,excelfile,variables, agg_type){
  for (year in years){
    data<-get_XYdata(species, year, excelfile,agg_type) 
    ind=which(colnames(data)=="Repeat")
    if (length(ind)>0){
      data<-data[,-ind]
    }
    if (year==years[1]){
      datayears<-data
    }else{
      datayears<-rbind(datayears,data)
    }
  }
  
  data0<-datayears
  
  for (varname in variables){
    data<-data0
    accessions<-unique(data$Accession)
    acctot<-data$Accession
    
    yind<-which(colnames(data)%in%varname)
    
    if (varname=="N" | varname=="CN"){
      hstart<-1500
      hend<-2400
    }else{
      hstart<-350
      hend<-2500
    }
    
    hstart<-paste0("mean_",hstart)
    hstart<-which(colnames(data)%in%hstart)
    
    hend<-paste0("mean_",hend)
    hend<-which(colnames(data)%in%hend)
    data<-cbind(data[,yind],data[,hstart:hend])
    colnames(data)[1]<-"y"
    nrep<-20
    training_corr2<-matrix(0,5,nrep)
    testing_corr2<-matrix(0,5,nrep)
    training_R2<-matrix(0,5,nrep)
    testing_R2<-matrix(0,5,nrep)
    test_result=list()
    train_result=list()
    for (r in 1:nrep){
      accfolds <- mycreateFolds(accessions, k = 5)
      print(paste0("Current repetition: ",r))
      train_result[[r]]=list()
      test_result[[r]]=list()
      for (f in 1:length(accfolds)){
        test_acc<-accessions[accfolds[[f]]]
        train_acc<-accessions[setdiff(1:length(accessions),accfolds[[f]])]
        test_ind<-which(acctot%in%test_acc)
        train_ind<-setdiff(1:nrow(data),test_ind)
        data_train<-data[train_ind,]
        data_test<-data[test_ind,]
        
        ## Permutation of data_test
        
        feat_test<-data_test[,2:ncol(data_test)]
        for (j in 1:nrow(feat_test)){
          feat_test[j,]<-sample(feat_test[j,],replace=FALSE)
        }
        
        data_test[2:ncol(data_test)]<-feat_test

        
        if (method=="PLSR"){
          acctot_train<-acctot[train_ind]
          ## inner CV fold
          inner_accfolds <- mycreateFolds(train_acc, k = 5)
          inner_folds_train <- list()
          inner_folds_test  <-list()
          for (i in 1:length(inner_accfolds)){
            inner_folds_test[[i]]<-which(acctot_train%in%train_acc[inner_accfolds[[i]]])
            inner_folds_train[[i]]<-setdiff(1:length(acctot_train),inner_folds_test[[i]])
          }
          # 
          # if (parallel){
          #   print("tuning in parallel")
          #   cl <- makePSOCKcluster(5)
          #   registerDoParallel(cl)
          # }
          my_control <- trainControl(method="cv",number=5,index=inner_folds_train,indexOut=inner_folds_test)
          model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,
                         trControl=my_control,tuneGrid = expand.grid(ncomp = 1:25))
          # if (parallel){
          #   stopCluster(cl)
          # }
        }else if (method=="SVR"){
          model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
        }
        
        pred_train<-predict(model,data_train)
        training_corr2[f,r]<-postResample(pred_train, data_train$y)[2]
        training_R2[f,r]<-get_R2(data_train$y,pred_train)
        pred_test<-predict(model,data_test)
        testing_corr2[f,r]<-postResample(pred_test, data_test$y)[2]
        testing_R2[f,r]<-get_R2(data_test$y,pred_test)
        train_result[[r]][[f]]=cbind(pred_train, data_train$y)
        test_result[[r]][[f]]=cbind(pred_test, data_test$y)
      }
    }
    if (length(years)>1){
      yearsp=c()
      for (y in years){
        yearsp=paste0(yearsp,"_",y)
      }
    }else{
      yearsp=paste0("_",years)
    }
    filename=paste0("results/scenario2/permutated_",species,"_",method,"_", varname,yearsp,"_",agg_type,".RData")
    save(train_result,test_result,testing_R2,testing_corr2,training_R2,training_corr2,file=filename)
    print(paste0("Variable: ",varname, yearsp, ",aggregated by: ", agg_type,", median testing R2: ", median(testing_R2)))
  }
}


prediction_scenario3_unseen <- function(method,species,testing_year,excelfile,variables, agg_type,tuned){
  
  data<-get_XYdata(species, testing_year, excelfile,agg_type) 
  ind=which(colnames(data)=="Repeat")
  if (length(ind)>0){
    data<-data[,-ind]
  }
  data_test0<-data
  
  years=list()
  years[[1]]=2021
  years[[2]]=2022
  years[[3]]=2023
  years[[4]]=c(2021,2022)
  years[[5]]=c(2021,2023)
  years[[6]]=c(2022,2023)
  
  check=c()
  for (y in years){
    if (sum(testing_year==y)==0){
      check=c(check,TRUE)
    }else{check=c(check,FALSE)}
  }
  
  training_years<-years[check]
  data_train0<-list()
  
  for (r in 1:length(training_years)){
    years<-training_years[[r]]
    for (year in years){
      data<-get_XYdata(species, year, excelfile,agg_type) 
      ind=which(colnames(data)=="Repeat")
      if (length(ind)>0){
        data<-data[,-ind]
      }
      if (excelfile=="lightSaturatedGasExchange" & all(years==c(2021,2022,2023))){
        ind<-which(colnames(data)%in%c("gs_sat","iWUE_sat"))
        data<-data[,-ind]
      }
      
      if (year==years[1]){
        datayears<-data
      }else{
        datayears<-rbind(datayears,data)
      }
    }
    data_train0[[r]]<-datayears
  }
  
  yearsp<-c()
  for(k in 1:length(training_years)){
    temp<-training_years[[k]]
    check<-temp[1]
    if (length(temp)>1){
      for (j in 2:length(temp)){
        check<-paste0(check,"&",temp[j])
      }
    }
    yearsp[k]=check
  }
  
  
  for (varname in variables){
    for (r in 1:length(training_years)){
      print(paste0("Training on year: ", training_years[r]))
      data_train_year<-data_train0[[r]]
      data_test_year<-data_test0
      
      
      acctot_train_year<-data_train_year$Accession
      uniq_acc_train_year<-unique(acctot_train_year)
      acctot_test_year<-data_test_year$Accession
      
      
      yind_train<-which(colnames(data_train_year)%in%varname)
      yind_test<-which(colnames(data_test_year)%in%varname)
      
      if (varname=="N" | varname=="CN"){
        hstart<-1500
        hend<-2400
      }else{
        hstart<-350
        hend<-2500
      }
      
      hstart1<-paste0("mean_",hstart)
      hstart1<-which(colnames(data_train_year)%in%hstart1)
      hend1<-paste0("mean_",hend)
      hend1<-which(colnames(data_train_year)%in%hend1)
      
      hstart2<-paste0("mean_",hstart)
      hstart2<-which(colnames(data_test_year)%in%hstart2)
      hend2<-paste0("mean_",hend)
      hend2<-which(colnames(data_test_year)%in%hend2)
      
      
      data_train_year<-cbind(data_train_year[,yind_train],data_train_year[,hstart1:hend1])
      data_test_year<-cbind(data_test_year[,yind_test],data_test_year[,hstart2:hend2])
      
      colnames(data_train_year)[1]<-"y"
      colnames(data_test_year)[1]<-"y"
      
      data_train_year<-data_train_year[,match(colnames(data_test_year),colnames(data_train_year))]
      
      nrep<-20
      training_corr<-matrix(0,5,nrep)
      testing_corr<-matrix(0,5,nrep)
      testing_count<-matrix(0,5,nrep)
      training_R2<-matrix(0,5,nrep)
      testing_R2<-matrix(0,5,nrep)
      test_result=list()
      train_result=list()
      for (rep in 1:nrep){
        accfolds <- mycreateFolds(uniq_acc_train_year, k = 5)
        print(paste0("Current repetition: ",rep))
        train_result[[rep]]=list()
        test_result[[rep]]=list()
        for (f in 1:5){
          print(paste0("Current fold: ",f))

          test_acc<-uniq_acc_train_year[accfolds[[f]]]
          train_acc<-uniq_acc_train_year[setdiff(1:length(uniq_acc_train_year),accfolds[[f]])]

          test_ind<-which(acctot_test_year%in%test_acc)
          train_ind<-which(acctot_train_year%in%train_acc)
          data_train<-data_train_year[train_ind,]
          data_test<-data_test_year[test_ind,]

          acctot_train_set<-acctot_train_year[train_ind]
          acctot_test_set<-acctot_test_year[test_ind]
          
          if (tuned){
            inner_accfolds <- mycreateFolds(train_acc, k = 5)
            inner_folds_test <- list()
            inner_folds_train <- list()
            
            for (j in 1:length(inner_accfolds)){
              inner_folds_test[[j]]<-which(acctot_train_set%in%train_acc[inner_accfolds[[j]]])
              inner_folds_train[[j]]<-setdiff(1:length(acctot_train_set),inner_folds_test[[j]])
            }
            
            my_control <- trainControl(method="cv",number=5,index=inner_folds_train,indexOut=inner_folds_test)
            if (method=="PLSR"){
              model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,
                             trControl=my_control,tuneGrid = expand.grid(ncomp = 1:50))
            }else if (method=="SVR"){
              model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,
                             trControl=my_control,tuneGrid = expand.grid(C = seq(0.25, 1, length = 4)))
            }
          }else{
            if (method=="PLSR"){
              model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,trControl=trainControl(method="none"))
            }else{
              model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
            }
          }
          
          test_acc22<-intersect(test_acc,acctot_test_year)
          data_test_final<-data_test_year[acctot_test_year%in%test_acc22,]
          
          pred_train<-predict(model,data_train)
          training_corr[f,rep]<-cor(pred_train, data_train$y)
          training_R2[f,rep]<-get_R2(data_train$y,pred_train)
          pred_test<-predict(model,data_test_final)
          testing_corr[f,rep]<-cor(pred_test, data_test_final$y)
          testing_R2[f,rep]<-get_R2(data_test_final$y,pred_test)
          testing_count[f,rep]<-length(pred_test)
          
          train_result[[rep]][[f]]=cbind(pred_train, data_train$y)
          test_result[[rep]][[f]]=cbind(pred_test, data_test_final$y)
        }
      }
      
      
      # colnames(testing_corr2)<-yearsp
      
      if (tuned & method=="PLSR"){
        method_lab="PLSR"
      }else if(tuned & method=="SVR"){
        method_lab="SVRtuned"
      }else if(!tuned & method=="PLSR"){
        method_lab="PLSRnotuned"
      }else if(!tuned & method=="SVR"){
        method_lab="SVR"
      }
      
      filename=paste0("results/scenario3/unseen_",species,"_",method_lab,"_", varname,"_testY",testing_year,"_trainY",yearsp[r],"_",agg_type,".RData")
      save(test_result,train_result,testing_R2,testing_corr,training_R2,training_corr,file=filename)
      print(paste0("Variable: ",varname, testing_year, ",aggregated by: ", agg_type,", median testing correlation: ", median(testing_corr)))
    }
    
  }
}

prediction_scenario3_seen <- function(method,species,testing_year,excelfile,variables, agg_type,tuned){
  
  data<-get_XYdata(species, testing_year, excelfile,agg_type) 
  ind=which(colnames(data)=="Repeat")
  if (length(ind)>0){
    data<-data[,-ind]
  }
  data_test0<-data
  
  years=list()
  years[[1]]=2021
  years[[2]]=2022
  years[[3]]=2023
  years[[4]]=c(2021,2022)
  years[[5]]=c(2021,2023)
  years[[6]]=c(2022,2023)
  
  check=c()
  for (y in years){
    if (sum(testing_year==y)==0){
      check=c(check,TRUE)
    }else{check=c(check,FALSE)}
  }
  
  training_years<-years[check]
  data_train0<-list()
  
  for (r in 1:length(training_years)){
    years<-training_years[[r]]
    for (year in years){
      data<-get_XYdata(species, year, excelfile,agg_type) 
      ind=which(colnames(data)=="Repeat")
      if (length(ind)>0){
        data<-data[,-ind]
      }
      if (excelfile=="lightSaturatedGasExchange" & all(years==c(2021,2022,2023))){
        ind<-which(colnames(data)%in%c("gs_sat","iWUE_sat"))
        data<-data[,-ind]
      }
      
      if (year==years[1]){
        datayears<-data
      }else{
        datayears<-rbind(datayears,data)
      }
    }
    data_train0[[r]]<-datayears
  }
  
  yearsp<-c()
  for(k in 1:length(training_years)){
    temp<-training_years[[k]]
    check<-temp[1]
    if (length(temp)>1){
      for (j in 2:length(temp)){
        check<-paste0(check,"&",temp[j])
      }
    }
    yearsp[k]=check
  }
  
  if (agg_type=="plot" & excelfile=="fluorescence" & all(years==c(2022))){
    ofold=13
    ifold=13
  }else if (agg_type=="plot" & excelfile=="fluorescence" & all(years==c(2023))){
    ofold=13
    ifold=13
  }else if (agg_type=="plot" & excelfile=="CN_isotopes" & all(years==c(2023))){
    ofold=13
    ifold=13
  }else if (agg_type=="plot" & excelfile=="SLA" & all(years==c(2023))){
    ofold=15
    ifold=15
  }else{
    ofold=10
    ifold=10
  }
  for (varname in variables){
    for (r in 1:length(training_years)){
      print(paste0("Training on year: ", training_years[r]))
      data_train_year<-data_train0[[r]]
      data_test_year<-data_test0
      
      
      acctot_train_year<-data_train_year$Accession
      uniq_acc_train_year<-unique(acctot_train_year)
      acctot_test_year<-data_test_year$Accession
      
      
      yind_train<-which(colnames(data_train_year)%in%varname)
      yind_test<-which(colnames(data_test_year)%in%varname)
      
      if (varname=="N" | varname=="CN"){
        hstart<-1500
        hend<-2400
      }else{
        hstart<-350
        hend<-2500
      }
      
      hstart1<-paste0("mean_",hstart)
      hstart1<-which(colnames(data_train_year)%in%hstart1)
      hend1<-paste0("mean_",hend)
      hend1<-which(colnames(data_train_year)%in%hend1)
      
      hstart2<-paste0("mean_",hstart)
      hstart2<-which(colnames(data_test_year)%in%hstart2)
      hend2<-paste0("mean_",hend)
      hend2<-which(colnames(data_test_year)%in%hend2)
      
      
      data_train_year<-cbind(data_train_year[,yind_train],data_train_year[,hstart1:hend1])
      data_test_year<-cbind(data_test_year[,yind_test],data_test_year[,hstart2:hend2])
      
      colnames(data_train_year)[1]<-"y"
      colnames(data_test_year)[1]<-"y"
      
      data_train_year<-data_train_year[,match(colnames(data_test_year),colnames(data_train_year))]
      
      nrep<-100
      training_corr<-matrix(0,1,nrep)
      testing_corr<-matrix(0,1,nrep)
      training_R2<-matrix(0,1,nrep)
      testing_R2<-matrix(0,1,nrep)
      test_result=list()
      train_result=list()
      for (rep in 1:100){
        print(paste0("Current repetition: ",rep))
        if (tuned){
          res<-samp_seen_lines_inner(data_train_year,acctot_train_year,ofold,ifold,TRUE)
          data_train<-res[[1]]
          data_test<-res[[2]]
          train_ind<-res[[3]]
          test_ind<-res[[4]]
          inner_folds_test<-res[[5]]
          inner_folds_train<-res[[6]]
          my_control <- trainControl(method="cv",number=5,index=inner_folds_train,indexOut=inner_folds_test)
          if (method=="PLSR"){   
            model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,
                           trControl=my_control,tuneGrid = expand.grid(ncomp = 1:50))
          }else if (method=="SVR"){
            model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,
                           trControl=my_control,tuneGrid = expand.grid(C = seq(0.25, 1, length = 4)))
          }
        }else{
          res<-samp_seen_lines_inner(data_train_year,acctot_train_year,ofold,ifold,FALSE)
          data_train<-res[[1]]
          data_test<-res[[2]]
          train_ind<-res[[3]]
          test_ind<-res[[4]]
          if (method=="PLSR"){   
            model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,trControl=trainControl(method="none"))
            # model<- plsr(y ~ ., ncomp = 2, data = data_train)
            
          }else if (method=="SVR"){
            model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
          }
          
        }
  
        test_acc22<-intersect(acctot_train_year[test_ind],acctot_test_year)
        data_test_final<-data_test_year[acctot_test_year%in%test_acc22,]
        
        pred_train<-predict(model,data_train)
        training_corr[rep]<-cor(pred_train, data_train$y)
        training_R2[rep]<-get_R2(data_train$y,pred_train)
        pred_test<-predict(model,data_test_final)
        testing_corr[rep]<-cor(pred_test, data_test_final$y)
        testing_R2[rep]<-get_R2(data_test_final$y,pred_test)
        
        train_result[[rep]]=cbind(pred_train, data_train$y)
        test_result[[rep]]=cbind(pred_test, data_test_final$y)
      }
      
      
      # colnames(testing_corr2)<-yearsp
      
      if (tuned & method=="PLSR"){
        method_lab="PLSR"
      }else if(tuned & method=="SVR"){
        method_lab="SVRtuned"
      }else if(!tuned & method=="PLSR"){
        method_lab="PLSRnotuned"
      }else if(!tuned & method=="SVR"){
        method_lab="SVR"
      }
      
      filename=paste0("results/scenario3/seen_",species,"_",method_lab,"_", varname,"_testY",testing_year,"_trainY",yearsp[r],"_",agg_type,".RData")
      save(test_result,train_result,testing_R2,testing_corr,training_R2,training_corr,file=filename)
      print(paste0("Variable: ",varname, testing_year, ",aggregated by: ", agg_type,", median testing correlation: ", median(testing_corr)))
    }
    
  }
}




prediction_scenario3_tuning <- function(method,species,testing_year,excelfile,variables, agg_type){
  
  data<-get_XYdata(species, testing_year, excelfile,agg_type) 
  ind=which(colnames(data)=="Repeat")
  if (length(ind)>0){
    data<-data[,-ind]
  }
  data_test0<-data
  
  years=list()
  years[[1]]=2021
  years[[2]]=2022
  years[[3]]=2023
  years[[4]]=c(2021,2022)
  years[[5]]=c(2022,2023)
  years[[6]]=c(2023,2021)
  years[[7]]=c(2021,2022,2023)
  
  check=c()
  for (y in years){
    if (sum(testing_year==y)==0){
      check=c(check,TRUE)
    }else{check=c(check,FALSE)}
  }
  
  training_years<-years[check]
  data_train0<-list()
  for (r in 1:length(training_years)){
    years<-training_years[[r]]
    for (year in years){
      data<-get_XYdata(species, year, excelfile,agg_type) 
      ind=which(colnames(data)=="Repeat")
      if (length(ind)>0){
        data<-data[,-ind]
      }
      if (excelfile=="lightSaturatedGasExchange" & all(years==c(2021,2022,2023))){
        ind<-which(colnames(data)%in%c("gs_sat","iWUE_sat"))
        data<-data[,-ind]
      }
      
      if (year==years[1]){
        datayears<-data
      }else{
        datayears<-rbind(datayears,data)
      }
    }
    
    data_train0[[r]]<-datayears
  }
  
  
  nrep<-20
  training_R2<-matrix(0,3,nrep)
  testing_R2<-matrix(0,3,nrep)
  training_corr2<-matrix(0,3,nrep)
  testing_corr2<-matrix(0,3,nrep)
  test_result=list()
  train_result=list()
  for (varname in variables){
    for (t in 1:length(training_years)){
      data<-data_train0[[t]]
      data_test<-data_test0
      accessions<-unique(data$Accession)
      acctot<-data$Accession
      
      yind_train<-which(colnames(data)%in%varname)
      yind_test<-which(colnames(data_test)%in%varname)
      
      if (varname=="N" | varname=="CN"){
        hstart<-1500
        hend<-2400
      }else{
        hstart<-350
        hend<-2500
      }
      
      hstart<-paste0("mean_",hstart)
      hstart<-which(colnames(data)%in%hstart)
      
      hend<-paste0("mean_",hend)
      hend<-which(colnames(data)%in%hend)
      
      data<-cbind(data[,yind_train],data[,hstart:hend])
      data_test<-cbind(data_test[,yind_test],data_test[,hstart:hend])
      
      colnames(data)[1]<-"y"
      colnames(data_test)[1]<-"y"
      data_train<-data
      
      ##
      
      for (r in 1:nrep){
        test_out_ind1<-sample(1:nrow(data_train),size=nrow(data_test)/10)
        train_in_ind1<-setdiff(1:nrow(data_train),test_out_ind1)
        
        test_out_ind2<-sample(1:nrow(data_test),size=nrow(data_test)/10)
        train_in_ind2<-setdiff(1:nrow(data_test),test_out_ind1)
        
        
        data_train<-data[train_ind,]
        data_test<-data[test_ind,]
        
        
        if (method=="PLSR"){
          acctot_train<-acctot[train_ind]
          ## inner CV fold
          inner_accfolds <- mycreateFolds(train_acc, k = 5)
          inner_folds_train <- list()
          inner_folds_test  <-list()
          for (i in 1:length(inner_accfolds)){
            inner_folds_test[[i]]<-which(acctot_train%in%train_acc[inner_accfolds[[i]]])
            inner_folds_train[[i]]<-setdiff(1:length(acctot_train),inner_folds_test[[i]])
          }
          # 
          # if (parallel){
          #   print("tuning in parallel")
          #   cl <- makePSOCKcluster(5)
          #   registerDoParallel(cl)
          # }
          my_control <- trainControl(method="cv",number=5,index=inner_folds_train,indexOut=inner_folds_test)
          model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,
                         trControl=my_control,tuneGrid = expand.grid(ncomp = 1:25))
          # if (parallel){
          #   stopCluster(cl)
          # }
        }else if (method=="SVR"){
          model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
        }
      }
      
      ##
      
      
      if (method=="PLSR"){
        model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,trControl=trainControl(method="none"))
      }else{
        model <- train(y ~ .,data = data_train,method="svmLinear", preProcess=NULL,trControl=trainControl(method="none"))
      }
      
      pred_train<-predict(model,data_train)
      training_corr2[t,r]<-postResample(pred_train, data_train$y)[2]
      training_R2[t,r]<-get_R2(data_train$y,pred_train)
      pred_test<-predict(model,data_test)
      testing_corr2[t,r]<-postResample(pred_test, data_test$y)[2]
      testing_R2[t,r]<-get_R2(data_test$y,pred_test)
      train_result[[r]][[t]]=cbind(pred_train, data_train$y)
      test_result[[r]][[t]]=cbind(pred_test, data_test$y)
    }
    
    yearsp<-c()
    for(k in 1:length(training_years)){
      temp<-training_years[[k]]
      check<-temp[1]
      if (length(temp)>1){
        for (j in 2:length(temp)){
          check<-paste0(check,"&",temp[j])
        }
      }
      yearsp[k]=check
    }
    
    colnames(testing_corr2)<-yearsp
    
    filename=paste0("results/scenario3/",species,"_",method,"_", varname,"_",testing_year,"_",agg_type,".RData")
    save(testing_R2,testing_corr2,training_R2,training_corr2,file=filename)
    print(paste0("Variable: ",varname, testing_year, ",aggregated by: ", agg_type,", median testing R2: ", median(testing_R2)))
  }
}



rewrite_light_sat2023 <- function(){
  setwd("C:/Users/Rudan/Documents/HyperspectralML")
  filen <- "data/traits/2023_AQcurves&lightSaturatedGasExchange_maize.csv"
  ydata <- read.csv(filen)
  newydata<-ydata[abs(ydata$PARi-1800)<5,]
  # uniq<-unique(newydata$plot_id)
  # count<-matrix(0,1,length(uniq))
  # for (i in 1:length(uniq)){
  #   count[i]<-sum(newydata$plot_id%in%uniq[i])
  # }
  write.csv(newydata,)
}