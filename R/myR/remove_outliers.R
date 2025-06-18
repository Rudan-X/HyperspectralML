remove_outliers <- function(variable,trait_HSR, data_type){
  if (variable=="d13C" & data_type=="plot_averaged"){
    trait_HSR <- trait_HSR[trait_HSR$y>(-30),]
  }else if (variable=="d13C" & data_type=="genotype_averaged"){
    trait_HSR <- trait_HSR[trait_HSR$y>(-20),]
  }
  if (variable=="Vmax" & data_type=="raw_data"){ # only happens in 2023
    trait_HSR <- trait_HSR[trait_HSR$y<75,]
  }
  
  if (variable=="C"){
    trait_HSR <- trait_HSR[trait_HSR$y<60,]
  }
  
  if (variable=="endNPQ"){
    trait_HSR <- trait_HSR[trait_HSR$y<0.5,]
  }
  
  if (variable=="NPQ_rel_res"){
    trait_HSR <- trait_HSR[trait_HSR$y<0.5,]
  }
  
  if (variable=="NPQ_ind_rate"){
    trait_HSR <- trait_HSR[trait_HSR$y<0.04,]
  }
  return(trait_HSR)
}
