setwd("C:/Users/Rudan/Documents/GitHub/HyperspectralML")
files=c("SLA","ACi","fluorescence","CN_isotopes")
for (year in c(2021,2022,2023)){ #
  
  # load reference data
  ref<-read.csv(paste0("data/accessions/CAM_Maize_",year,"_ref.csv"))
  uniq_ref<-unique(ref[,colnames(ref)%in%c("Plot","Accession")])
  
  # load HSR data
  HSR <- read.csv(paste0("data/reflectance/", year, "_hyperspectral_reflectance_maize.csv"))
  HSR_plot<-aggregate(HSR, by=list(HSR$Plot), "mean")
  ind<-which(colnames(HSR_plot)%in%c("Group.1","Repeat"))
  HSR_plot <- HSR_plot[,-ind]
  HSR_plot0<-merge(uniq_ref,HSR_plot,by="Plot")
  
  HSR_acc<-aggregate(HSR_plot0, by=list(HSR_plot0$Accession), "mean")
  ind<-which(colnames(HSR_acc)%in%c("Accession","Plot"))
  HSR_acc <- HSR_acc[,-ind]
  colnames(HSR_acc)[colnames(HSR_acc)%in%"Group.1"]<-"Accession"
  
  HSR$ID<-paste0(HSR$Plot,"_",HSR$Repeat)
  
  
  ind<-which(colnames(HSR)%in%c("Plot","Repeat"))
  HSR <- HSR[,-ind]
  
  
  # load traits
  
  for (f in 1:length(files)){
    excelfile<-files[f]
    filen <- paste("data/traits/", year, "_", excelfile, "_maize.csv", sep = "")
    ydata <- read.csv(filen)
    
    if (excelfile == "CN_isotopes") {
      ydata <- ydata[complete.cases(ydata), ]
      ydata$CN <- as.numeric(ydata$CN)
      ind<-which(is.na(ydata),arr.ind = T)
      if (dim(ind)[1]!=0){
        ydata<-ydata[-ind[,2],]
      }
    }
    ydata<-ydata[complete.cases(ydata),]
    colnames(ydata)[c(1,2)]<-c("Plot","Repeat")
    ydata$Plot <- as.character(ydata$Plot)
    
    # calculate its average across plots
    temp_plot<-aggregate(ydata, by=list(ydata$Plot), "mean")
    ind<-which(colnames(temp_plot)%in%c("Plot","Repeat"))
    temp_plot <- temp_plot[,-ind]
    colnames(temp_plot)[colnames(temp_plot)%in%"Group.1"]<-"Plot"
    
    temp_plot0<-merge(uniq_ref,temp_plot,by="Plot")
    # calculate its average across accession
    temp_acc<-aggregate(temp_plot0, by=list(temp_plot0$Accession), "mean")
    ind<-which(colnames(temp_acc)%in%c("Accession","Plot"))
    temp_acc <- temp_acc[,-ind]
    colnames(temp_acc)[colnames(temp_acc)%in%"Group.1"]<-"Accession"

    
    ydata$ID<-paste0(ydata$Plot,"_",ydata$Repeat)
    ind<-which(colnames(ydata)%in%c("Plot","Repeat"))
    ydata <- ydata[,-ind]
    
    
    if (f==1){
      traits<-ydata
      traits_plot<-temp_plot
      traits_acc<-temp_acc
      
    }else{
      if (year==2021){
        traits<-merge(ydata,traits,by="ID", all=T)
      }else{
        if (excelfile!="CN_isotopes"){
          traits<-merge(ydata,traits,by="ID", all=T)
        }
      }
      
      traits_plot<-merge(temp_plot,traits_plot,by="Plot", all=T)
      traits_acc<-merge(temp_acc,traits_acc,by="Accession", all=T)
    }
  }

  traits_plot<-merge(uniq_ref,traits_plot,by="Plot")
  
  traitHSR<-merge(traits,HSR,by="ID",all=T)
  traitHSR$ID<-paste0(year,"_",traitHSR$ID)
  traitHSR<-merge(ref,traitHSR,by="ID")
  
  removeind<-which(is.na(traitHSR$mean_350))
  if (length(removeind)>0){
    traitHSR<-traitHSR[-removeind,]
  }
  
  
  traitHSR_plot<-merge(traits_plot,HSR_plot,by="Plot",all=T)
  removeind<-which(is.na(traitHSR_plot$mean_350))
  if (length(removeind)>0){
    traitHSR_plot<-traitHSR_plot[-removeind,]
  }
  
  traitHSR_acc<-merge(traits_acc,HSR_acc,by="Accession",all=T)
  removeind<-which(is.na(traitHSR_acc$mean_350))
  if (length(removeind)>0){
    traitHSR_acc<-traitHSR_acc[-removeind,]
  }
  
  ##############################
  
  write.csv(traitHSR,paste0("data/combined_data/traits_and_HSR",year,"_raw_data.csv"),row.names = F, quote=F)
  write.csv(traitHSR_plot,paste0("data/combined_data/traits_and_HSR",year,"_plot_averaged.csv"),row.names = F, quote=F)
  write.csv(traitHSR_acc,paste0("data/combined_data/traits_and_HSR",year,"_genotype_averaged.csv"),row.names = F, quote=F)
  
}


# some plots for CN isotopes from 2022 and 2023 come with already averaged data
