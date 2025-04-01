var1<-c("SLA","C","N","CN","Area","C","d13C","d15N",
        "Vmax","Vpmax","Vpmax.Vmax","SL","a400","gsw","iWUE","A_sat","gs_sat","ci_sat","iWUE_sat","gs","sl",
        "NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res",
        "phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear",
        "maxNPQ","endNPQ","endFvFm","initialFvFm")
var2<-c("SLA","%C","%N","C/N","Leaf area","%C","d13C","d15N",
        "Vmax","Vpmax","Vpmax/Vmax","SL","A\n (sat-light)","Gsw\n (sat-light)","iWUE\n (sat-light)","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
        "ΦPSII\n induction\n amplitude","ΦPSII\n induction\n rate","ΦPSII\n recovery\n offset","NPQ\n induction\n slope",
        "MaxNPQ","Final\n NPQ","Final\n Fv/Fm","Fv/Fm")

varclass <- c(rep("Structural & biochemical",8),
              rep("Gas exchange",13),
              rep("Chlorophyll fluorescence",13))
