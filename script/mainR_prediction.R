library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pls)
library(caret)
library(readxl)
library(e1071)
library(reshape2)
library(randomForest)
setwd("/home/rudan/Documents/CAPITALISE_ML")
species <- "barley"
year <- 2022
excelfile <- "fluorescence"
agg=FALSE
fullrange=TRUE
tuning="5fold"
plsrmethod="oscorespls"
method="glmnet"
scena=3

source("script/prediction_functions.R")

#######################Generating results##############################
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
prediction("maize",scena,method,"5fold","6fold",FALSE,1)

prediction("barley",scena,method,"5fold","6fold",FALSE,1)



# prediction("barley",1,"glmnet","5fold","rep2fold",1)

# Predict for scenario 3

files=c("SLA","fluorescence","lightSaturatedGasExchange")
# methods=c("PLSR","glmnet","svm")
# for (scena in 1:2){
#   for (method in methods){
#     prediction("maize",scena,method,"5fold","rep2fold")
#   }
# }

## Plot comparison between methods

## Plot comparison between methods
methods=list()
if (species=="barley"){
  methods[[1]]<-c("PLSR","svm","glmnet") #
  methods[[2]]<-c("5fold","5fold","5fold")#
  methodslab<-c("PLSR","SVM","GLMNET") #"PLSR",
}else{
  methods[[1]]=c("PLSR","PLSR","svm","glmnet") #
  methods[[2]]=c("LOO","5fold","5fold","5fold")#
  methodslab=c("PLSR","PLSR","SVM","GLMNET") #"PLSR",
}


for (scena in 1:2){
  p1=plot_method_comparison(methods,species, year,scena,"testing",1,methodslab)
  p2=plot_method_comparison(methods,species, year,scena,"testing",2,methodslab)
  p3=plot_method_comparison(methods,species, year,scena,"testing",3,methodslab)
  p4=plot_method_comparison(methods,species, year,scena,"testing",4,methodslab)

  figure=ggarrange(p1, p2, p3,p4,ncol = 2, nrow = 2,common.legend = TRUE, legend = "bottom")
  annotate_figure(figure,
                  bottom = text_grob("Traits"),
                  left = text_grob("R2",  rot = 90))

  filen=paste0("figures/webinar/boxplot_",species,"_scena",scena,"compare_methods.png")
  dev.print(file=filen, device=png, width=1600,height = 900,res=160)
}




## Compare different outer fold
method="svm"
tuning="5fold"
species="barley"
year=2022
p1=plot_outerCV_comparison(species, year,"testing",method,tuning,1)
p2=plot_outerCV_comparison(species, year,"testing",method,tuning,2)
p3=plot_outerCV_comparison(species, year,"testing",method,tuning,3)
p4=plot_outerCV_comparison(species, year,"testing",method,tuning,4)

figure=ggarrange(p1, p2, p3,p4,ncol = 2, nrow = 2,common.legend = TRUE, legend = "bottom")
annotate_figure(figure,
                bottom = text_grob("Traits"),
                left = text_grob("R2",  rot = 90))

filen=paste0("figures/webinar/Comparison_outerCV_",species,"_",method,"_",tuning,".png")
dev.print(file=filen, device=png, width=1600,height = 900,res=160)



### Plot scenario3 comparison between original and agg data of PLSR and SVM

p1=plot_outerCV_scena3(species,"testing","5fold",1)
p2=plot_outerCV_scena3(species,"testing","5fold",2)

figure=ggarrange(p1, p2,ncol = 2, nrow = 1,common.legend = TRUE, legend = "bottom")
annotate_figure(figure,
                bottom = text_grob("Traits"),
                left = text_grob("Correlation",  rot = 90))

filen="figures/webinar/Comparison_outerCV_scena3.png"
dev.print(file=filen, device=png, width=1600,height = 900,res=160)

###
files=c("SLA","ACi","CN_isotopes","fluorescence","lightSaturatedGasExchange" ) #,"fluorescence","lightSaturatedGasExchange","SLA","fluorescence","lightSaturatedGasExchange"
years=c(2021,2021,2021, 2021,2021 ) 


for (i in 1:3){
  plot_scenario1and2("maize", years[i], 1, "svm", files[i], '5fold')  
}



## Plot previous 2folds results

for (i in 1:length(files)){
  plot_2folds_results("maize", years[i], 1, "PLSR", files[i], 'LOO')  
}



##
