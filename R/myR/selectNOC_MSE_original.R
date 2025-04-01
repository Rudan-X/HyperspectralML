##' @param PRESS PRESS matrix (iterations x maxComps)
##' @param method Which approach to use to find optimal components. Options: pls, firstPlateau, firstMin
##' @param maxComps maximum number of components to consider
##' @return nComps the optimal number of PLSR components
##' @return pressDFres the dataframe to plot selection of NOC

selectNOC_MSE_original <- function(MSEP, select_strat="hastie",sdfact=1,segments=3) {
  MSEPm <- apply(MSEP, 2, mean)
  MSEPsd <- apply(MSEP, 2, sd)/sqrt(segments)
  if (select_strat == "diffnext") {
    fvec <- (diff(MSEPm) + sdfact * MSEPsd[-1]) < 
      0
    fvec <- c(TRUE, fvec)
    ind <- which.min(MSEPm)
    optcomp <- max((1:ind)[fvec[1:ind]])
  }else if (select_strat == "hastie") {
    ind <- which.min(MSEPm)
    fvec <- (MSEPm < (MSEPm[ind] + sdfact * MSEPsd[ind]))
    optcomp <- min((1:ind)[fvec[1:ind]])
  }else if (select_strat == "relchange") {
    ind <- which.min(MSEPm)
    MSEPsel <- MSEPm[1:ind]
    relchange <- (MSEPsel - MSEPm[ind])/max(MSEPsel) > 
      0.001
    ind2 <- which.max((1:length(relchange))[relchange])
    MSEPm2 <- MSEPsel[1:ind2]
    MSEPsd2 <- MSEPsd[1:ind2]
    indm <- which.min(MSEPm2)
    fvec <- (MSEPm2 < (MSEPm2[indm] + sdfact * MSEPsd2[indm]))
    optcomp <- min((1:indm)[fvec[1:indm]])
  }
  
  return(list(bestnoc=optcomp))
}