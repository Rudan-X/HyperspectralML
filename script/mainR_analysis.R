library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pls)
library(caret)
library(readxl)
library(e1071)
library(reshape2)
setwd("/home/rudan/Documents/CAPITALISE_ML")

source("script/prediction_functions.R")
source("script/plotting_functions.R")

## Plot comparison between methods
methods=list()

methods[[1]]<-c("PLSR","svm") #
methods[[2]]<-c("5fold","5fold")#
# methodslab<-c("PLSR","SVM") #"PLSR",


################ Figure 2 and 3 ########################
datat<-"Original"#"Original""Averaged"
plot_main_goodfit(methods,"maize",2021,"testing",c(.05, .1),datat)
filen=paste0("figures/publication/","maize","_",datat,"_boxplot_good_fit.png")
dev.print(file=filen, device=png, width=3200,height = 1400,res=200)



datat<-"Averaged" # "Original"
plot_main_goodfit(methods,"barley",2022,"testing",c(.05, .1),datat)
filen=paste0("figures/publication/","barley","_",datat,"_boxplot_good_fit.png")
dev.print(file=filen, device=png, width=3000,height = 1400,res=200)



################ scatter plot ##############################
datat<-"Original"#"Original""Averaged"
plot_scatter_goodfit(methods,"maize","testing",c(.05, .1),datat)
filen=paste0("figures/publication/","maize","_",datat,"_boxplot_good_fit.png")
dev.print(file=filen, device=png, width=3200,height = 1400,res=200)


################ Figure 4 Robustness ########################
# Comparison between original and averaged data
source("script/plotting_functions.R")
plot_main_robustness(methods,"testing")
filen<-"figures/publication/Figure4_robustness.png"
dev.print(file=filen, device=png, width=3200,height = 1600,res=200)


################# SUPPLEMENTARY FIGURES #####################
source("script/plotting_functions.R")
datat<- "Original"
plot_supp_badfit(methods,"maize", 2021,"testing",c(0.03,0.1),"Original")
filen=paste0("figures/publication/","maize","_",datat,"_boxplot_bad_fit.png")
dev.print(file=filen, device=png, width=3200,height = 1400,res=200)

plot_supp_badfit(methods,"maize", 2022,"testing",c(0.03,0.1),"Original")
filen=paste0("figures/publication/","maize","_",datat,"2022_boxplot_bad_fit.png")
dev.print(file=filen, device=png, width=3200,height = 1400,res=200)

plot_supp_badfit(methods,"barley", 2022,"testing",c(0.03,0.1),"Original")
filen=paste0("figures/publication/","barley","_",datat,"_boxplot_bad_fit.png")
dev.print(file=filen, device=png, width=3200,height = 1400,res=200)


methods=list()

methods[[1]]<-c("PLSR","svm","glmnet","rf","svmPoly") #
methods[[2]]<-rep("5fold",5)#


source("script/plotting_functions.R")
plot_supp_MLmethods(methods,"testing")
filen<-"figures/publication/supp_MLmethods.png"
dev.print(file=filen, device=png, width=3200,height = 1600,res=200)
