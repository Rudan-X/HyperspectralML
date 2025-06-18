##' @param PRESS PRESS matrix (iterations x maxComps)
##' @param method Which approach to use to find optimal components. Options: pls, firstPlateau, firstMin
##' @param maxComps maximum number of components to consider
##' @return nComps the optimal number of PLSR components
##' @return pressDFres the dataframe to plot selection of NOC

selectNOC_MSE <- function(MSEP, select_strat="hastie",sdfact=1, repl=30,segments=3) {
  optcomp <- matrix(NA, nrow = 1, ncol = repl)
  for (i in 1:repl){
    MSEPm <- apply(MSEP[((i-1)*segments+1):(i*segments),], 2, mean)
    MSEPsd <- apply(MSEP[((i-1)*segments+1):(i*segments),], 2, sd)/sqrt(segments)
    if (select_strat == "diffnext") {
      fvec <- (diff(MSEPm) + sdfact * MSEPsd[-1]) < 0
      fvec <- c(TRUE, fvec)
      ind <- which.min(MSEPm)
      optcomp[,i] <- max((1:ind)[fvec[1:ind]])
    }else if (select_strat == "hastie") {
      ind <- which.min(MSEPm)
      fvec <- (MSEPm < (MSEPm[ind] + sdfact * MSEPsd[ind]))
      final_ind <- min((1:ind)[fvec[1:ind]])
      if (final_ind==Inf){
        optcomp[,i] <- ind
      }else{
        optcomp[,i] <- final_ind
      }
    }else if (select_strat == "relchange") {
      ind <- which.min(MSEPm)
      MSEPsel <- MSEPm[1:ind]
      relchange <- (MSEPsel - MSEPm[ind])/max(MSEPsel) > 0.001
      ind2 <- which.max((1:length(relchange))[relchange])
      MSEPm2 <- MSEPsel[1:ind2]
      MSEPsd2 <- MSEPsd[1:ind2]
      indm <- which.min(MSEPm2)
      fvec <- (MSEPm2 < (MSEPm2[indm] + sdfact * MSEPsd2[indm]))
      optcomp[,i] <- min((1:indm)[fvec[1:indm]])
    }
  }
  
  afinaldistr <- table(optcomp)/sum(table(optcomp))
  afinal <- as.numeric(names(which.max(afinaldistr)))
  afinaldistr<-melt(afinaldistr)
  
  
  bp<-ggplot(afinaldistr, aes(x = optcomp, y = value)) +
    geom_point() +  geom_line() + 
    geom_vline(xintercept = afinal,linetype = "dashed",  
               color = "blue")+
    labs(x = "Number of components",  
         y = "Relative frequency for optimal number") + 
    scale_x_continuous(breaks = seq(min(optcomp), max(optcomp), by = 1)) +
    theme_minimal() 
  print(bp)
  
  return(list(bestnoc=afinal,plot=bp))
}