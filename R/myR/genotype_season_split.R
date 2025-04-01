genotype_season_split <- function(test_dataset,train_dataset, k) {
  genotypes <- unique(test_dataset$genotype)
  genokfolds <- random_split(genotypes, k)
  trainfolds <- list()
  testfolds <- list()
  for (j in 1:length(genokfolds)){
    testfolds[[j]] <- which(test_dataset$genotype%in%genotypes[genokfolds[[j]]])
    
    valind <- which(train_dataset$genotype%in%genotypes[genokfolds[[j]]])
    trainfolds[[j]] <- setdiff(seq(1,nrow(train_dataset)),valind)
  }
  return(list(trainfold=trainfolds,testfold=testfolds))
}
