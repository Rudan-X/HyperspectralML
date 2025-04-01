##' Inner loop structure: permutations
##' @param dataset input full PLSR dataset. Usually just the calibration dataset
##' @param maxComps maximum number of components to use for each PLSR fit
##' @param iterations how many different permutations to run
##' @param prop proportion of data to preserve for each permutation
##' @param verbose Should the function report the current iteration status/progress to the terminal
##' or run silently? TRUE/FALSE. Default FALSE
##' @return output a list containing the PRESS, MSEP and coef_array.
innerloop_permutations <- function(dataset=NULL, maxComps=20, method, iterations=20, prop=0.70,verbose=FALSE) {
  coefs <- array(0,dim=c((ncol(dataset$Spectra)+1),iterations,maxComps))
  press.out <- array(data=NA, dim=c(iterations,maxComps))
  MSEP.out <- array(data=NA, dim=c(iterations,maxComps))
  print("*** Running permutation test.  Please hang tight, this can take awhile ***")
  print("Options:")
  print(paste("Max Components:",maxComps, "Iterations:", iterations, 
              "Data Proportion (percent):", prop*100, sep=" "))
  
  if (verbose) {
    j <- 1 # <--- Numeric counter for progress bar
    pb <- utils::txtProgressBar(min = 0, max = iterations, 
                                char="*",width=70,style = 3)
  }
  
  for (i in seq_along(1:iterations)) {
    rows <- sample(1:nrow(dataset),floor(prop*nrow(dataset)))
    sub.data <- dataset[rows,]
    val.sub.data <- dataset[-rows,]
    plsr.out <- plsr(as.formula("y~Spectra"), scale=FALSE, center=TRUE, 
                     ncomp=maxComps,method=method, validation="none", data=sub.data)
    pred_val <- predict(plsr.out,newdata=val.sub.data)
    sq_resid <- (pred_val[,,]-val.sub.data[,"y"])^2
    press <- apply(X = sq_resid, MARGIN = 2, FUN = sum)
    press.out[i,] <- press
    MSEP.out[i,]<-apply(sq_resid, 2, mean)
    coefs[,i,] <- coef(plsr.out, intercept = TRUE, ncomp = 1:maxComps)
    rm(rows,sub.data,val.sub.data,plsr.out,pred_val,sq_resid,press)
    
    ### Display progress to console
    if (verbose) {
      setTxtProgressBar(pb, j)    # show progress bar
      j <- j+1                    # <--- increase counter by 1
      flush.console()             #<--- show output in real-time
    }
  }
  if (verbose) {
    close(pb)
  }
  
  # create a new list with PRESS and permuted coefficients x wavelength x component number
  print("*** Providing PRESS and coefficient array output ***")
  output <- list(PRESS=press.out, MSEP=MSEP.out,coef_array=coefs)
  return(output)
}