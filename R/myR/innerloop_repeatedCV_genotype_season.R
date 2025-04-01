##' Inner loop structure: repeatedCV
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param method method for PLSR model
##' @param iterations how many different permutations to run
##' @param segments segments of data to split during each iteration
##' @return output a list containing the PRESS, MSEP and coef_array.


innerloop_repeatedCV_genotype_season <- function(dataset=NULL, maxComps=20, method = "oscorespls", iterations=30, segments=3,parallel=FALSE) {
  
  MSEP.out <- matrix(NA,nrow = segments*iterations , ncol = maxComps)
  press.out <- matrix(NA,nrow = segments*iterations , ncol = maxComps)
  optcomp <- matrix(NA, nrow = 1, ncol = iterations)
  for (i in 1:iterations) {
    
    gfolds <- genotype_split(dataset=dataset,segments)
    if (parallel){
      pls.options(parallel = makeCluster(segments))
      res <- plsr(as.formula("y~Spectra"), data=dataset, ncomp = maxComps,  
                  method = method, segments = gfolds, validation = "CV")
    }else{
      res <- plsr(as.formula("y~Spectra"), data=dataset, ncomp = maxComps,  
              method = method, segments = gfolds, validation = "CV")
    }
    # MSEPj <- matrix(NA, nrow = segments, ncol = maxComps)
    MSEPj <- matrix(NA,nrow = segments , ncol = maxComps)
    pressj <- matrix(NA,nrow = segments , ncol = maxComps)
    for (j in 1:segments) {
      # for each tested segment:
      test_ind <- res$vali$seg[[j]]
      predj <- res$vali$pred[test_ind, , ]
      obsj <- dataset$y[test_ind]
      sq_resid <- (predj - obsj)^2
      MSEPj[j, ] <- apply(sq_resid, 2, mean)
      pressj[j, ] <- apply(sq_resid, 2, sum)
    }
    
    MSEP.out[((i-1)*3+1):(i*3), ] <- MSEPj
    press.out[((i-1)*3+1):(i*3), ] <- pressj
    
    # optcomp[,i]<-selectNOC_MSE(MSEPj, selstrat="hastie",sdfact=1, maxComps=20)
  }
  return(list(PRESS=press.out, MSEP=MSEP.out))
}