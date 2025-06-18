rm(list = ls())
list.of.packages <- c("pls","dplyr","reshape2","parallel","ggplot2","caret")
invisible(lapply(list.of.packages, library, character.only = TRUE))

path0 <- "/home/mpimp-golm.mpg.de/xu2004/HyperspectralML/"
if (!dir.exists(path0)) {
  path0 <- "C:/Users/Rudan/Documents/GitHub/HyperspectralML/"
}


setwd(path0)
for (file in c("R_Burnett/","myR/")){
  filen<-paste0(path0,"R/",file)
  myls <- list.files(path=filen,pattern="*.R")
  myls<-paste0(filen,myls)
  invisible(sapply(myls,FUN=source))
}



segments0 <- 1
nrep <- 1

repl<-10
segment<-3

agg <- "sampledHSR"


vars<-c("SLA","N","Vmax","a400","phiPSII_ind_res","initialFvFm") #
year1 <- c(2022,2023,2021,2023,2021,2022)
year2 <- c(2021,2021,2022,2022,2023,2023)

y <- 3
inVar <- "N"
if (inVar=="N" | inVar=="CN"){
  endd <- 2
}else{
  endd <- 3
}

datatype <- "genotype_averaged"
traitHSR1 <- read.csv(paste0("data/combined_data/traits_and_HSR",year1[y],"_",datatype,".csv"))
traitHSR1$Season <- year1[y]
traitHSR2 <- read.csv(paste0("data/combined_data/traits_and_HSR",year2[y],"_",datatype,".csv"))
traitHSR2$Season <- year2[y]
sametrait <- intersect(colnames(traitHSR1),colnames(traitHSR2))
traitHSR <- rbind(traitHSR1[,sametrait],traitHSR2[,sametrait])

if (inVar=="N" | inVar=="CN"){
  Start.wave <- 1500
  End.wave <- 2400
}else{
  Start.wave <- 400
  End.wave <- 2400
}

if (agg=="fullHSR"){
  wv <- seq(Start.wave,End.wave,1)
}else if(agg=="sampledHSR"){
  wv <- seq(Start.wave,End.wave,by=5)
}

spec <- as.matrix(traitHSR[, which(names(traitHSR) %in%  paste0("mean_",wv))])
traitHSR <- cbind(data.frame(season=traitHSR$Season,genotype=traitHSR$Accession,y=traitHSR[,inVar]),spec)
traitHSR <- traitHSR[!is.na(traitHSR$y),]


templ <- list(train = matrix(0,segments0,nrep), test = matrix(0,segments0,nrep))
temp <- matrix(0,segments0,nrep)
temp2 <- matrix(0,nrow(traitHSR),1)
res.rcv.svr <- list(R2 = templ, RMSEP = templ, RMSEPper = templ, bestnoc = temp, pred = temp2)


val.ind <- which(traitHSR$season==year2[y])
cal.ind <- setdiff(seq(1,nrow(traitHSR)),val.ind)
val.data <- traitHSR[val.ind, ]
cal.data <- traitHSR[cal.ind, ]
cal_spec <- as.matrix(cal.data[, which(names(cal.data) %in%  paste0("mean_",wv))])
cal.data <- data.frame(y=cal.data$y, Spectra=I(cal_spec))

val_spec <- as.matrix(val.data[, which(names(val.data) %in% paste0("mean_",wv))])
val.data <- data.frame(y=val.data$y, Spectra=I(val_spec))

datasets<-list(cal.data,val.data)


val_spec0 <- val_spec
cal_spec0 <- cal_spec

acc1 <- traitHSR$genotype[cal.ind]
acc2 <- traitHSR$genotype[val.ind]

acc <- intersect(acc1,acc2)
accdist <- setdiff(acc2,acc1)

ind1 <- match(acc,traitHSR$genotype[cal.ind])
ind2 <- match(acc,traitHSR$genotype[val.ind])

plot(traitHSR$y[cal.ind[ind1]],traitHSR$y[val.ind[ind2]])
cor(traitHSR$y[cal.ind[ind1]],traitHSR$y[val.ind[ind2]])


ind3 <- match(accdist,traitHSR$genotype[val.ind])


val_spec <- as.data.frame(val_spec0[ind2,])
colnames(val_spec) <- seq(1500,2400,by=5)
val_spec$genotype <- traitHSR$genotype[val.ind[ind2]]
df <- melt(val_spec,id.vars = "genotype")
df$season <- "Common genotypes (2022)"


cal_spec <- as.data.frame(cal_spec0[ind1,])
colnames(cal_spec) <- seq(1500,2400,by=5)
cal_spec$genotype <- traitHSR$genotype[cal.ind[ind1]]
df2 <- melt(cal_spec,id.vars = "genotype")
df2$season <- "Common genotypes (2021)"


val_spec <- as.data.frame(val_spec0[ind3,])
colnames(val_spec) <- seq(1500,2400,by=5)
val_spec$genotype <- traitHSR$genotype[val.ind[ind3]]
df3 <- melt(val_spec,id.vars = "genotype")
df3$season <- "Distinct genotypes (2022)"


dfall <- rbind(df,df2)
dfall <- rbind(dfall,df3)

dfall <- dfall[!is.na(dfall$variable),]
dfall$genotype <- factor(dfall$genotype,levels=c(acc,accdist))
color_palette <- colorRampPalette(c( "red","orange", "yellow", "green", "cyan", "blue", "purple"))(318)


ggplot(dfall, aes(x=variable,y=value,color=genotype)) +
  facet_grid(. ~ season, scales = "free")+
  geom_point(size=0.05)+
  theme(legend.position = "none",  axis.title=element_text(size=16,face="bold"), 
        strip.text = element_text(size=15,face="bold"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        ) +     # pure white plot background) +
  labs(y="Reflectance (%)",x="Wavelength (nm)") +
  scale_color_manual(values = color_palette)


ggsave("Figures/SF_sec4.4_HSR_2122.png",width = 14, height = 7)

# 
# ind0 <-  which(datasets[[1]]$Spectra[,"mean_1660"]>min(datasets[[2]]$Spectra[,"mean_1660"]))
# 
# model.svr <- train(as.formula("y~Spectra"),data = datasets[[1]],method="svmLinear", preProcess = c("center", "scale"),
#                    trControl=trainControl(method="none"))
# 
# 
# ind4a <- which(datasets[[2]]$Spectra[,"mean_1660"]<max(datasets[[1]]$Spectra[,"mean_1660"]))
# ind4b <- which(datasets[[2]]$Spectra[,"mean_1660"]>min(datasets[[1]]$Spectra[,"mean_1660"]))
# 
# t <- 2
# pred.svr<-predict(model.svr, newdata=datasets[[t]])
# observed <- datasets[[t]]$y
# plot(pred.svr[ind4], observed[ind4])
# 
# cor(pred.svr[ind4], observed[ind4])
# VIP <- function(object) {
#   SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
#   Wnorm2 <- colSums(object$loading.weights^2)
#   SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*") # Replace with matrix mult.
#   sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
# }
# 
# 
# vips <- VIP(model.mse)
# 
# check <- data.frame(values=c(vips[7,],datasets[[1]]$Spectra[1,]),x=c(wv,wv),type=c(rep("VIP",length(wv)),rep("HSR",length(wv))))
# 
# plot(wv,vips[7,],xlab="Wavelength (nm)",ylab="VIP")
# 
# plot(wv,datasets[[1]]$Spectra[1,],xlab="Wavelength (nm)",ylab="VIP")
# 
# 
# ggplot(check, aes(x=x,y=values,color=type)) +
#   geom_point(size=0.5)
