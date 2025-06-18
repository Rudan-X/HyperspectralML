vars<-c( "SLA","N", "CN", "d13C",  "Vpmax","Vmax", "a400", "gsw",  
         "NPQ_ind_amp", "NPQ_rel_amp", "NPQ_rel_rate","maxNPQ", 
         "phiPSII_ind_amp", "phiPSII_ind_rate", "phiPSII_ind_res", "endFvFm", "initialFvFm")


dtsi<-c( rep(3,3), 2,rep(1,4), rep(1,6),rep(2,3))
MLs <- c(rep("PLSR",3),rep("SVR",5),rep("PLSR",9))
dtlabels <- c("Raw data","Plot averaged","Genotype averaged")

vars2 <- long_var2_flat[match(vars,var1)] 

tab <- data.frame(Trait=vars2,Aggregation=dtlabels[dtsi],MLmethod=MLs)

library("openxlsx")
write.xlsx(
  x = list("Trait_specific" = tab),
  file = "results/supp_tables/Trait_aggregation_ML_combination.xlsx"
)
