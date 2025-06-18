reshape_mat <- function(arr){
  reshaped_list <- lapply(1:dim(arr)[3], function(i) {
    slice <- arr[,,i]  # 5 x 20
    matrix(slice, ncol = 1)  # flatten to 100 x 1
  })
  
  # Combine each 100x1 column to get 100 x 40 matrix
  result <- do.call(cbind, reshaped_list)  # 100 x 40
  return(result)
}



calc_mean <- function(mat){
  temp <- apply(mat, c(2, 3), mean)
  colnames(temp) <- seq(1,ncol(temp))
  df <- melt(temp)
  df <- df[,2:3]
  colnames(df) <- c("NoC", "score")
  return(df)
}


get_nocR <- function(mat,noc){
  sel_R <- matrix(NA, dim(mat)[1], dim(mat)[2])
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      sel_R[i,j] <- mat[i,j,noc[1,j]]
    }
  }
  return(sel_R)
}
