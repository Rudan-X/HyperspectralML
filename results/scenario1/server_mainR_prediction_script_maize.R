library(dplyr)
library(tidyr)
library(ggplot2)
library(pls)
library(doParallel)
library(caret)

setwd("/work/xu2/CAPITALIZE_ML/")

source("script/prediction_functions.R")

prediction=function(species,scena,method,tuning,outerCV,parallel,file_ind){
  spec_ind<-list()
  if (species=="maize"){
    files=c("ACi","CN_isotopes","fluorescence") #,"SLA","fluorescence","lightSaturatedGasExchange"
    years=c(2021,2021,2021,2021)#,2022,2022,2022
    spec_ind[[1]]=c("Vpmax","Vmax")
    spec_ind[[2]]=c("N","CN")
    spec_ind[[3]]=c("phiPSII_ind_res","maxNPQ")
  }else{
    files=c("ACi","CN_isotopes" ,"fluorescence") #"LeafArea","lightSaturatedGasExchange"
    years=rep(2022,4)
    spec_ind[[1]]=c("Vcmax","Jmax","TPU")
    spec_ind[[2]]=c("N","CN")
    spec_ind[[3]]=c("phiPSII_ind_res","maxNPQ")
  }
  for (i in file_ind){
    if (outerCV=="6fold"){
      predict_scenario1and2_6fold(species, years[i], scena, method, files[i], tuning,parallel,spec_ind[[i]])
    }else if (outerCV=="rep2fold"){
      predict_scenario1and2_rep2fold(species, years[i], scena, method, files[i], tuning,parallel)
    }
  }
}

method="svmPoly"
scena=1
paral=TRUE
prediction("maize",scena,method,"5fold","6fold",paral,2:3)
method="rf"
prediction("maize",scena,method,"5fold","6fold",paral,2:3)
