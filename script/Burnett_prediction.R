
# setwd("/work/xu2/HyperspectralML/")
setwd("C:/Users/Rudan/Documents/HyperspectralML_complete/")
list.of.packages <- c("pls","dplyr","plotrix","ggplot2","gridExtra")
invisible(lapply(list.of.packages, library, character.only = TRUE))

### Setup other functions and options
# not in
`%notin%` <- Negate(`%in%`)
# Script options
pls::pls.options(plsralg = "oscorespls")
pls::pls.options("plsralg")

# Default par options
opar <- par(no.readonly = T)

filen<-"C:/Users/Rudan/Documents/HyperspectralML_complete/TESTgroup-BNL-spectratrait-3182844/R/"
myls <- list.files(path=filen,pattern="*.R")
myls<-paste0(filen,myls)
sapply(myls,FUN=source)

source("script/Rfunctions_data.R")

aggmethods=c("none","plot","acc") 
files=c("ACi","SLA","fluorescence")#"CN_isotopes",
vars<-list()
# vars[[1]]<-c("N","CN")
vars[[1]]<-c("Vpmax","Vmax","gsw","a400")#,"gs_sat","ci_sat","iWUE_sat"
vars[[2]]<-c("SLA")
vars[[3]]=c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")


######################################################################
### Create plsr dataset
year=2023
for (agg in aggmethods){
  for (f in 1:3){
    plsr_data<-get_XYdata("maize", year,files[f],agg) 
    
    # What is the target variable?
    for (inVar in vars[[f]]){
      Start.wave <- 350
      End.wave <- 2500
      wv <- seq(Start.wave,End.wave,1)
      ############### Create cal/val datasets #####################
      
      method <- "dplyr" #base/dplyr
      # base R - a bit slow
      # dplyr - much faster
      split_data <- create_data_split(dataset=plsr_data, approach=method, 
                                      split_seed=7529075, prop=0.8)
      
      cal.plsr.data <- split_data$cal_data
      
      
      intersect(split_data$cal_data$Plant_Species,split_data$val_data$Plant_Species)
      cal.plsr.data <- split_data$cal_data
      val.plsr.data <- split_data$val_data
      
      
      rm(split_data)
      
      # Datasets:
      
      ## Step 8. Create calibration and validation PLSR datasets
      ### Format PLSR data for model fitting 
      cal_spec <- as.matrix(cal.plsr.data[, which(names(cal.plsr.data) %in% 
                                                    paste0("mean_",wv))])
      cal.plsr.data <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% 
                                                          paste0("mean_",wv))], 
                                  Spectra=I(cal_spec))
      
      val_spec <- as.matrix(val.plsr.data[, which(names(val.plsr.data) %in% 
                                                    paste0("mean_",wv))])
      val.plsr.data <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% 
                                                          paste0("mean_",wv))],
                                  Spectra=I(val_spec))
      
      
      
      ## Step 10. Use permutation to determine the optimal number of components
      
      ### Use permutation to determine the optimal number of components
      if(grepl("Windows", sessionInfo()$running)){
        pls.options(parallel = NULL)
      } else {
        pls.options(parallel = parallel::detectCores()-1)
      }
      
      method <- "firstMin" #pls, firstPlateau, firstMin
      random_seed <- 7529075
      seg <- 80
      maxComps <- 20
      iterations <- 50
      prop <- 0.70
      
      png(file=file.path(paste0("Burnett_result/",inVar,year,"_",agg,"_NoC_selection.png")),
          height=2800, width=4800,  res=340)
      if (method=="pls") {
        nComps <- find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                          method=method, 
                                          maxComps=maxComps, seg=seg, 
                                          random_seed=random_seed)
        print(paste0("*** Optimal number of components: ", nComps))
      } else {
        nComps <- find_optimal_components(dataset=cal.plsr.data,
                                          method=method, 
                                          maxComps=maxComps, iterations=iterations, 
                                          seg=seg, prop=prop, 
                                          random_seed=random_seed)
      }
      dev.off()
      
      ## Step 11. Fit final model
      
      ### Fit final model - using leave-one-out cross validation
      plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,
                       validation="CV",trace=FALSE,data=cal.plsr.data)
      
      fit <- plsr.out$fitted.values[,1,nComps]
      pls.options(parallel = NULL)
      
      # External validation fit stats
      text_loc <- c(max(RMSEP(plsr.out, newdata = val.plsr.data)$comps),
                    RMSEP(plsr.out, newdata = val.plsr.data)$val[1])
      
      png(file=file.path(paste0("Burnett_result/",inVar,year,"_",agg,"_Validation_RMSEP_R2_by_Component.png")),
                         height=2800, width=4800,  res=340)
      par(mfrow=c(1,2)) # B, L, T, R
      pls::RMSEP(plsr.out, newdata = val.plsr.data)
      
      plot(pls::RMSEP(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL RMSEP",
           xlab="Number of Components",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
      text(text_loc[1],text_loc[2],labels = "", cex=2)
      box(lwd=2.2)
      
      pls::R2(plsr.out, newdata = val.plsr.data)
      
      plot(pls::R2(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL R2",
           xlab="Number of Components",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)
      box(lwd=2.2)
      
      
      dev.off()
      ## Step 12. PLSR fit observed vs. predicted plot data
      #calibration
      cal.plsr.output <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin%
                                                            "Spectra")],
                                    PLSR_Predicted=fit,
                                    PLSR_CV_Predicted=as.vector(plsr.out$validation$pred[,,
                                                                                         nComps]))
      
      cal.plsr.output <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin%
                                                            "Spectra")],
                                    PLSR_Predicted=fit,
                                    PLSR_CV_Predicted=fit)
      
      cal.plsr.output <- cal.plsr.output %>%
        mutate(PLSR_CV_Residuals = PLSR_CV_Predicted-get(inVar))
      head(cal.plsr.output)
      
      
      cal.R2 <- round(pls::R2(plsr.out,intercept=F)[[1]][nComps],2)
      cal.RMSEP <- round(sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2)),2)
      
      val.plsr.output <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% 
                                                            "Spectra")],
                                    PLSR_Predicted=as.vector(predict(plsr.out, 
                                                                     newdata = val.plsr.data, 
                                                                     ncomp=nComps, 
                                                                     type="response")[,,1]))
      val.plsr.output <- val.plsr.output %>%
        mutate(PLSR_Residuals = PLSR_Predicted-get(inVar))
      head(val.plsr.output)
      
      
      val.R2 <- round(pls::R2(plsr.out,newdata=val.plsr.data,intercept=F)[[1]][nComps],2)
      val.RMSEP <- round(sqrt(mean(val.plsr.output$PLSR_Residuals^2)),2)
      
      rng_quant <- quantile(cal.plsr.output[,inVar], probs = c(0.001, 0.999))
      cal_scatter_plot <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Predicted, y=get(inVar))) + 
        theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                                linetype="dashed", size=1.5) + 
        xlim(rng_quant[1], rng_quant[2]) + 
        ylim(rng_quant[1], rng_quant[2]) +
        labs(x=paste0("Predicted ", paste(inVar), " (units)"),
             y=paste0("Observed ", paste(inVar), " (units)"),
             title=paste0("Calibration: ", paste0("Rsq = ", cal.R2), "; ", 
                          paste0("RMSEP = ", cal.RMSEP))) +
        theme(axis.text=element_text(size=18), legend.position="none",
              axis.title=element_text(size=20, face="bold"), 
              axis.text.x = element_text(angle = 0,vjust = 0.5),
              panel.border = element_rect(linetype = "solid", fill = NA, size=1.5)) + 
        annotate("text", x=rng_quant[1], y=rng_quant[2], label= "",size=10)
      
      cal_resid_histogram <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Residuals)) +
        geom_histogram(alpha=.5, position="identity") + 
        geom_vline(xintercept = 0, color="black", 
                   linetype="dashed", size=1) + theme_bw() + 
        theme(axis.text=element_text(size=18), legend.position="none",
              axis.title=element_text(size=20, face="bold"), 
              axis.text.x = element_text(angle = 0,vjust = 0.5),
              panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
      
      rng_quant <- quantile(val.plsr.output[,inVar], probs = c(0.001, 0.999))
      val_scatter_plot <- ggplot(val.plsr.output, aes(x=PLSR_Predicted, y=get(inVar))) + 
        theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                                linetype="dashed", size=1.5) + 
        xlim(rng_quant[1], rng_quant[2]) + 
        ylim(rng_quant[1], rng_quant[2]) +
        labs(x=paste0("Predicted ", paste(inVar), " (units)"),
             y=paste0("Observed ", paste(inVar), " (units)"),
             title=paste0("Validation: ", paste0("Rsq = ", val.R2), "; ", 
                          paste0("RMSEP = ", val.RMSEP))) +
        theme(axis.text=element_text(size=18), legend.position="none",
              axis.title=element_text(size=20, face="bold"), 
              axis.text.x = element_text(angle = 0,vjust = 0.5),
              panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
      
      val_resid_histogram <- ggplot(val.plsr.output, aes(x=PLSR_Residuals)) +
        geom_histogram(alpha=.5, position="identity") + 
        geom_vline(xintercept = 0, color="black", 
                   linetype="dashed", size=1) + theme_bw() + 
        theme(axis.text=element_text(size=18), legend.position="none",
              axis.title=element_text(size=20, face="bold"), 
              axis.text.x = element_text(angle = 0,vjust = 0.5),
              panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))
      
      png(file.path(paste0("Burnett_result/",inVar,year,"_",agg,"_scatterplot_cal_vs_val.png")), 
               height=2800, width=4800,  res=340)
      # plot cal/val side-by-side
      scatterplots <- grid.arrange(cal_scatter_plot, val_scatter_plot, cal_resid_histogram, 
                                   val_resid_histogram, nrow=2, ncol=2)
      
      dev.off()
      
      if(grepl("Windows", sessionInfo()$running)){
        pls.options(parallel =NULL)
      } else {
        pls.options(parallel = parallel::detectCores()-1)
      }
      
      jk.plsr.out <- pls::plsr(as.formula(paste(inVar,"~","Spectra")), scale=FALSE, 
                               center=TRUE, ncomp=nComps, validation="CV", trace=FALSE, 
                               jackknife=TRUE, 
                               data=cal.plsr.data)
      pls.options(parallel = NULL)
      
      Jackknife_coef <- f.coef.valid(plsr.out = jk.plsr.out, 
                                     data_plsr = cal.plsr.data, 
                                     ncomp = nComps, inVar=inVar)
      Jackknife_intercept <- Jackknife_coef[1,,,]
      Jackknife_coef <- Jackknife_coef[2:dim(Jackknife_coef)[1],,,]
      
      interval <- c(0.025,0.975)
      Jackknife_Pred <- val.plsr.data$Spectra %*% Jackknife_coef + 
        matrix(rep(Jackknife_intercept, length(val.plsr.data[,inVar])), byrow=TRUE, 
               ncol=length(Jackknife_intercept))
      
      Interval_Conf <- apply(X = Jackknife_Pred, MARGIN = 1, FUN = quantile, 
                             probs=c(interval[1], interval[2]))
      sd_mean <- apply(X = Jackknife_Pred, MARGIN = 1, FUN =sd)
      sd_res <- sd(val.plsr.output$PLSR_Residuals)
      sd_tot <- sqrt(sd_mean^2+sd_res^2)
      val.plsr.output$LCI <- Interval_Conf[1,]
      val.plsr.output$UCI <- Interval_Conf[2,]
      val.plsr.output$LPI <- val.plsr.output$PLSR_Predicted-1.96*sd_tot
      val.plsr.output$UPI <- val.plsr.output$PLSR_Predicted+1.96*sd_tot
      head(val.plsr.output)

      
      ### Permutation validation plot
      rmsep_percrmsep <- percent_rmse(plsr_dataset = val.plsr.output, 
                                      inVar = inVar, 
                                      residuals = val.plsr.output$PLSR_Residuals, 
                                      range="full")
      RMSEP <- rmsep_percrmsep$rmse
      perc_RMSEP <- rmsep_percrmsep$perc_rmse
      r2 <- round(pls::R2(plsr.out, newdata = val.plsr.data,intercept=F)$val[nComps],2)
      expr <- vector("expression", 3)
      expr[[1]] <- bquote(R^2==.(r2))
      expr[[2]] <- bquote(RMSEP==.(round(RMSEP,2)))
      expr[[3]] <- bquote("%RMSEP"==.(round(perc_RMSEP,2)))
      rng_vals <- c(min(val.plsr.output$LPI), max(val.plsr.output$UPI))
      
      png(file.path(paste0("Burnett_result/",inVar,year,"_",agg,"_PLSR_Validation_Scatterplot.png")), 
          height=2800, width=3200,  res=340)
      
      par(mfrow=c(1,1), mar=c(4.2,5.3,1,0.4), oma=c(0, 0.1, 0, 0.2))
      plotrix::plotCI(val.plsr.output$PLSR_Predicted,val.plsr.output[,inVar], 
                      li=val.plsr.output$LPI, ui=val.plsr.output$UPI, gap=0.009,sfrac=0.004, 
                      lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
                      err="x", pch=21, col="black", pt.bg=scales::alpha("grey70",0.7), scol="grey50",
                      cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
                      ylab=paste0("Observed ", paste(inVar), " (units)"),
                      cex.axis=1.5,cex.lab=1.8)
      abline(0,1,lty=2,lw=2)
      legend("topleft", legend=expr, bty="n", cex=1.5)
      legend("bottomright", legend="", bty="n", cex=2.2)
      box(lwd=2.2)
      
      dev.off()
      
      filename=paste0("Burnett_result/", inVar,year,"_",agg,".RData")
      save(plsr.out,jk.plsr.out,nComps,cal.plsr.output,val.plsr.output,val.plsr.output,cal.plsr.data,val.plsr.data,file=filename)
      
    }
  }
}


