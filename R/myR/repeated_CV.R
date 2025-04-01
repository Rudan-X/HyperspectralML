##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param iterations how many different permutations to run
##' @param sdfact The parsimony factor p controls the selection of the optimum
# number of components. A value of zero would choose the global
# minimum of MSE; sdfact=1 gives the one standard error rule; sdfact=2
# gives a two standard error rule and considers a 95% confidence
# interval of the minimum MSE, resulting in a model with a very
# small number of components. Mostly, sdfact=1 gives best results
##' @param segments the number of folds 
##' @param verbose Should the function report the current iteration status/progress to the terminal
##' or run silently? TRUE/FALSE. Default FALSE
##' @return output of optimal NoC


repeated_CV <- function(dataset=NULL, maxComps=20, method = "oscorespls", iterations=100, segments=3, 
                        selstrat = " ",sdfact=1,verbose=FALSE) {
  
  
  optcomp <- matrix(NA, nrow = 1, ncol = iterations)
  for (i in 1:iterations) {
    
    res <- plsr(as.formula("y~Spectra"), data=dataset, ncomp = maxComps,  
              method = method, segments = segments, validation = "CV")
    # MSEPj <- matrix(NA, nrow = segments, ncol = maxComps)
    MSEPj <- matrix(NA,nrow = segments , ncol = maxComps)
    for (j in 1:segments) {
      # for each tested segment:
      test_ind <- res$vali$seg[[j]]
      predj <- res$vali$pred[test_ind, , ]
      obsj <- dataset$y[test_ind]
      resj <- predj - obsj
      MSEPj[j, ] <- apply(resj^2, 2, mean)
    }
    MSEPm <- apply(MSEPj, 2, mean)
    MSEPsd <- apply(MSEPj, 2, sd)/sqrt(segments)
    if (selstrat == "diffnext") {
      fvec <- (diff(MSEPm) + sdfact * MSEPsd[-1]) < 
        0
      fvec <- c(TRUE, fvec)
      ind <- which.min(MSEPm)
      optcomp[, i] <- max((1:ind)[fvec[1:ind]])
    }else if (selstrat == "hastie") {
      ind <- which.min(MSEPm)
      fvec <- (MSEPm < (MSEPm[ind] + sdfact * MSEPsd[ind]))
      optcomp[, i] <- min((1:ind)[fvec[1:ind]])
    }else if (selstrat == "relchange") {
      ind <- which.min(MSEPm)
      MSEPsel <- MSEPm[1:ind]
      relchange <- (MSEPsel - MSEPm[ind])/max(MSEPsel) > 
        0.001
      ind2 <- which.max((1:length(relchange))[relchange])
      MSEPm2 <- MSEPsel[1:ind2]
      MSEPsd2 <- MSEPsd[1:ind2]
      indm <- which.min(MSEPm2)
      fvec <- (MSEPm2 < (MSEPm2[indm] + sdfact * MSEPsd2[indm]))
      optcomp[, i] <- min((1:indm)[fvec[1:indm]])
    }
  }
  afinaldistr <- table(optcomp)/sum(table(optcomp))
  afinal <- as.numeric(names(which.max(afinaldistr)))
  return(afinal)
}