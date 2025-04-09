
rm(list = ls())

path0 <- "/home/mpimp-golm.mpg.de/xu2004/HyperspectralML/"
if (!dir.exists(path0)) {
  path0 <- "C:/Users/Rudan/Documents/GitHub/HyperspectralML/"
}


setwd(path0)

vars<-c("SLA","C",  "N",  "CN", "d13C","d15N","Vpmax","Vmax","SL", "a400","gsw","Vpmax.Vmax","iWUE" ,"NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate",    
        "NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res" , "NPQ_ind_linear",  
        "maxNPQ","endNPQ","endFvFm","initialFvFm")

cormat<-matrix(0,26,3)
ncount<-matrix(0,26,3)
rownames(cormat)<-vars


years1 <- c(2021,2021,2022)
years2 <- c(2022,2023,2023)

for (y in 1:3){
  traits1<- read.csv(paste0("data/combined_data/traits_and_HSR",years1[y],"_genotype_averaged.csv"))
  traits2<- read.csv(paste0("data/combined_data/traits_and_HSR",years2[y],"_genotype_averaged.csv"))
  
  ind <- which(colnames(traits1)%in%"mean_350")
  traits1 <- traits1[,seq(1,ind-1)]
  ind <- which(colnames(traits2)%in%"mean_350")
  traits2 <- traits2[,seq(1,ind-1)]
  

  
  for (v in 1:length(vars)){
    trait1_v <- data.frame(Accession=traits1[,"Accession"], value=traits1[,vars[v]])
    trait2_v <- data.frame(Accession=traits2[,"Accession"], value=traits2[,vars[v]])
    
    rownames(trait1_v) <- trait1_v$Accession
    rownames(trait2_v) <- trait2_v$Accession
    acc <- intersect(trait1_v$Accession,trait2_v$Accession)
    trait1_v <- trait1_v[acc,]
    trait2_v <- trait2_v[acc,]
    
    ind <- intersect(which(complete.cases(trait1_v)),which(complete.cases(trait2_v)))
    
    trait1_v <- trait1_v[ind,]
    trait2_v <- trait2_v[ind,]
    
    cormat[v,y]<-cor(trait1_v$value,trait2_v$value,method="pearson")
    ncount[v,y]<-nrow(trait1_v)
  }

}

meancor <- round(apply(cormat,1,mean),3)
cormat <- round(cormat,3)

sort(meancor)

