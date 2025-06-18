vegitation_indices <- function(HSR,method){
  if (method=="PRI"){ #Linked to xanthophyll cycle and leaf structure
    index <- (HSR[,"mean_531"]-HSR[,"mean_570"])/(HSR[,"mean_531"]+HSR[,"mean_570"])
  }else if (method=="mSR705"){
    index <- (HSR[,"mean_750"]-HSR[,"mean_445"])/(HSR[,"mean_705"]+HSR[,"mean_445"])
  }else if (method=="EVI"){
    index <- 2.5 * (HSR[,"mean_780"] - HSR[,"mean_675"])/(HSR[,"mean_782"] + 6 * HSR[,"mean_675"] - 7.5 * HSR[,"mean_445"] + 1)
  }else if (method=="PSRI"){
    index <- (HSR[,"mean_680"] - HSR[,"mean_500"])/HSR[,"mean_750"]
  }else if (method=="ARI"){
    index <- HSR[,"mean_800"] * (1/HSR[,"mean_550"] - 1/HSR[,"mean_700"])
  }else if (method=="SIPI"){
    index <- (HSR[,"mean_800"]-HSR[,"mean_445"])/(HSR[,"mean_800"]+HSR[,"mean_680"])
  }else if (method=="SGI"){
    ind <- which(colnames(HSR)%in%"mean_500")
    ind2 <- which(colnames(HSR)%in%"mean_600")
    index <- apply(HSR[,ind:ind2],1,mean)
  }else if (method=="NDVI"){ # red edge NDVI, captures chlorophyll and mesophyll structure
    index <- (HSR[,"mean_750"]-HSR[,"mean_705"])/(HSR[,"mean_750"]+HSR[,"mean_705"])
  }else if (method=="REM"){ #Reflectance at 800 nm (NIR), 700 nm (red-edge region)
    index <- HSR[,"mean_800"]/HSR[,"mean_700"]-1
  }else if (method=="NDMI-1"){
    index <- HSR[,"mean_860"]- (HSR[,"mean_1640"]-HSR[,"mean_2130"])/HSR[,"mean_860"]+(HSR[,"mean_1640"]-HSR[,"mean_2130"])
  }else if (method=="NDII"){
    index <- (HSR[,"mean_860"]-HSR[,"mean_1650"])/(HSR[,"mean_860"]+HSR[,"mean_1650"])
  }
}
