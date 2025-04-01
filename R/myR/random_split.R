##' Generate k-folds
##' @param strat_id index to split
##' @param k k-fold
##' @return output a list containing the validation folds
##' 
random_split <- function(strat_id, k) {
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }	
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))	
  return(folds)
}