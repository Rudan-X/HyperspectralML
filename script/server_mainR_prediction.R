library(dplyr)
library(tidyr)
library(ggplot2)
library(pls)
library(doParallel)
library(caret)

# setwd("/work/xu2/HyperspectralML/")
setwd("C:/Users/Rudan/Documents/HyperspectralML")

source("script/Rfunctions_data.R")

yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=c(2021,2022)
yearsused[[4]]=2023
yearsused[[5]]=c(2021,2023)
yearsused[[6]]=c(2022,2023)
yearsused[[7]]=c(2021,2022,2023)
types=c("none","plot","acc") #c("none","plot","acc")

files=c("CN_isotopes","ACi","SLA","fluorescence","lightSaturatedGasExchange","ACimaxCO2") #,"leaf_width"
vars=list()
# vars[[1]]=c("N","CN") # "N","CN","d15N","C","d13C"
# vars[[2]]=c("Vpmax","Vmax","gsw","a400") #"Vpmax","Vmax","gsw","a400","SL","iWUE"
# vars[[3]]=c("SLA")
# vars[[4]]=c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

vars<-list()
vars[[1]]<-c("d15N","C","d13C")
vars[[2]]<-c("Vpmax.Vmax","SL","iWUE")#,"gs_sat","ci_sat","iWUE_sat"
vars[[4]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","NPQ_ind_linear","endNPQ")


# vars[[4]]=c("A_sat","ci_sat") #"gs_sat","iWUE_sat"
# vars[[5]]=c("A_maxCO2")



# vars[[1]]=c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res","phiPSII_ind_amp",
#             "phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear","maxNPQ", "endNPQ","endFvFm","initialFvFm")

# start y=4, t=2

## Prediction 
for (method in c("SVR")){ #,"SVR""PLSR"
  for (t in 2:2){ #length(types)
    aggtype=types[t]
    for (y in 1:1){ #length(yearsused)
      years<-yearsused[[y]]
      for (f in 1:2){
        prediction_scenario1(method,"maize",years,files[f],vars[[f]], aggtype,FALSE)

      }
    }
  }
}



### Scenario 3 ####

yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=2023


for (method in c("SVR","PLSR")){ #
  for (t in 2:2){ #length(types)
    aggtype=types[t]
    for (y in 1:3){ #length(yearsused)
      years<-yearsused[[y]]
      for (f in 1:4){
        if (method=="SVR"){
          tuning<-F
        }else{
          tuning<-T
        }
        prediction_scenario3_seen(method,"maize",years,files[f],vars[[f]], aggtype,tuning)
        prediction_scenario3_unseen(method,"maize",years,files[f],vars[[f]], aggtype,tuning)
      }
    }
  }
}



files=c("CN_isotopes","ACi","SLA","lightSaturatedGasExchange","fluorescence") #
vars=list()
vars[[1]]=c("N","CN") # "N",,"CN","d15N","C","d13C"
vars[[2]]=c("Vpmax","Vmax","a400") #"Vpmax","Vmax","gsw","a400","SL","iWUE"
vars[[3]]=c("SLA")
vars[[4]]=c("A_sat")
vars[[5]]=c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

files=c("fluorescence") #
vars=list()
vars[[1]]=c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res","phiPSII_ind_amp",
                        "phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear","maxNPQ", "endNPQ","endFvFm","initialFvFm")
yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=2023


for (method in c("PLSR")){ #"SVR",
  for (t in 1:1){ #length(types)
    aggtype=types[t]
    for (y in 1:1){ #length(yearsused)
      years<-yearsused[[y]]
      for (f in 1:4){
        prediction_scenario3(method,"maize",years,files[f],vars[[f]], aggtype)
      }
    }
  }
}

##### PROXIS ##########

## Prediction 
# for (method in c("SVR")){
#   for (t in 3:3){ #length(types)
#     aggtype=types[t]
#     for (y in 1:1){ #length(yearsused)
#       years<-yearsused[[y]]
#       prediction_scenario1_organic(method,"maize",years,"ACimaxCO2","A_maxPAR", aggtype)
#     }
#   }
# }
yearsused=list()
yearsused[[1]]=2021
yearsused[[2]]=2022
yearsused[[3]]=c(2021,2022)

file<-"lightSaturatedGasExchange"
var<-"A_sat"
types=c("none","plot","acc") #c("none","plot","acc")
for (method in c("SVR")){
  for (t in 1:3){ #length(types)
    aggtype=types[t]
    for (y in 1:3){ #length(yearsused)
      years<-yearsused[[y]]
      # prediction_scenario2_organic(method,"maize",years,file,var, aggtype) #,"ACimaxCO2","A_maxCO2","ACi","Vmax"
      # prediction_scenario2_fluorescence(method,"maize",years,file,var, aggtype)
      # prediction_scenario2_both_proxies(method,"maize",years,file,var, aggtype)
      prediction_scenario2_proxies_and_reflectance(method,"maize",years,file,var, aggtype)
    }
  }
}


for (method in c("SVR")){
  for (t in 1:3){ #length(types)
    aggtype=types[t]
    for (y in 1:3){ #length(yearsused)
      years<-yearsused[[y]]
      prediction_scenario2_organic_and_reflectance(method,"maize",years,"ACi","Vmax", aggtype) #,"ACimaxCO2","A_maxCO2"
    }
  }
}
