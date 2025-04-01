genotype_split <- function(dataset, k) {
  genotypes <- unique(dataset$genotype)
  genokfolds <- random_split(genotypes, k)
  kfolds <- list()
  for (k in 1:length(genokfolds)){
    kfolds[[k]] <- which(dataset$genotype%in%genotypes[genokfolds[[k]]])
  }
  return(kfolds)
}
