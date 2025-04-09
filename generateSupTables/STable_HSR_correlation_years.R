
rm(list = ls())

path0 <- "/home/mpimp-golm.mpg.de/xu2004/HyperspectralML/"
if (!dir.exists(path0)) {
  path0 <- "C:/Users/Rudan/Documents/GitHub/HyperspectralML/"
}


setwd(path0)


years1 <- c(2021,2021,2022)
years2 <- c(2022,2023,2023)

cormat<-matrix(0,3,length(seq(350,2500)))
colnames(cormat) <- paste0("mean_",seq(350,2500))

for (y in 1:3){
  HSR1<- read.csv(paste0("data/combined_data/traits_and_HSR",years1[y],"_genotype_averaged.csv"))
  HSR2<- read.csv(paste0("data/combined_data/traits_and_HSR",years2[y],"_genotype_averaged.csv"))
  
  rownames(HSR1) <- HSR1$Accession
  rownames(HSR2) <- HSR2$Accession
  
  ind <- which(colnames(HSR1)%in%"mean_350")
  HSR1 <- HSR1[,seq(ind,ncol(HSR1))]
  ind <- which(colnames(HSR2)%in%"mean_350")
  HSR2 <- HSR2[,seq(ind,ncol(HSR2))]

  acc <- intersect(rownames(HSR1),rownames(HSR2))
  
  HSR1 <- HSR1[acc,]
  HSR2 <- HSR2[acc,]
  
  cormat[y,] <- mapply(function(x, y) cor(x, y), HSR1,HSR2)
}


hstart<-which(colnames(cormat)%in%paste0("mean_",350))
hend<-which(colnames(cormat)%in%paste0("mean_",700))

# round(mean(cormat[,hstart:hend]),3)

round(mean(cormat[1,hstart:hend]),3)
round(mean(cormat[2,hstart:hend]),3)
round(mean(cormat[3,hstart:hend]),3)


hstart<-which(colnames(cormat)%in%paste0("mean_",700))
hend<-which(colnames(cormat)%in%paste0("mean_",1300))

round(mean(cormat[1,hstart:hend]),3)
round(mean(cormat[2,hstart:hend]),3)
round(mean(cormat[3,hstart:hend]),3)

hstart<-which(colnames(cormat)%in%paste0("mean_",1300))
hend<-which(colnames(cormat)%in%paste0("mean_",2500))

round(mean(cormat[,hstart:hend]),3)
round(mean(cormat[2,hstart:hend]),3)
round(mean(cormat[3,hstart:hend]),3)
