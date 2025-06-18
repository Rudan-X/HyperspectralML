var1<-c("SLA","C","N","CN","Area","C","d13C","d15N",
        "Vmax","Vpmax","Vpmax.Vmax","SL","a400","gsw","iWUE","A_sat","gs_sat","ci_sat","iWUE_sat","gs","sl",
        "NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res",
        "phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear",
        "maxNPQ","endNPQ","endFvFm","initialFvFm","VcmaxT_P25", "VpmaxT_P25", "JmaxT_P25")

var2<-c("SLA","%C","%N","C:N","Leaf area","%C","d13C","d15N",
        "Vmax","Vpmax","Vpmax/Vmax","SL","A\n","gsw\n","iWUE\n","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n offset",
        "ΦPSII\n recovery\n amplitude","ΦPSII\n recovery\n rate","ΦPSII\n recovery\n offset","NPQ\n induction\n slope",
        "MaxNPQ","Final\n NPQ","Final\n ΦPSII","Fv/Fm","VcmaxT_P25", "VpmaxT_P25", "JmaxT_P25")

short_var2<-c("SLA","%C","%N","C:N","Leaf area","%C","d13C","d15N",
        "Vmax","Vpmax","Vpmax/Vmax","SL","A","gsw","iWUE","A at\n sat light","Gs at\n sat light","Ci at\n sat light","iWUE at\n sat light",
        "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
        "NPQ\n ind.\n amp.","NPQ\n ind.\n rate","NPQ\n rel.\n amp.","NPQ\n rel.\n rate","NPQ\n rel.\n offset",
        "ΦPSII\n rec.\n amp.","ΦPSII\n rec.\n rate","ΦPSII\n rec.\n offset","NPQ\n ind.\n slope",
        "Max.\n NPQ","Final\n NPQ","Final\n ΦPSII","Fv/Fm","VcmaxT_P25", "VpmaxT_P25", "JmaxT_P25")


var2_flat<-c("SLA","%C","%N","C:N","Leaf area","%C","d13C","d15N",
             "Vmax","Vpmax","Vpmax/Vmax","SL","A","Gsw","iWUE","A at sat light","Gs at sat light","Ci at sat light","iWUE at sat light",
             "Gs at sat light sat CO2","Stomatal limitation",
             "NPQ ind. amp.","NPQ ind. rate","NPQ rel. amp.","NPQ rel. rate","NPQ rel. offset",
             "ΦPSII rec. amp.","ΦPSII rec. rate","ΦPSII rec. offset","NPQ ind. slope",
             "Max. NPQ","Final NPQ","Final ΦPSII","Fv/Fm","VcmaxT_P25", "VpmaxT_P25", "JmaxT_P25")

long_var2_flat<-c("SLA","%C","%N","C:N","Leaf area","%C","d13C","d15N",
             "Vmax","Vpmax","Vpmax/Vmax","SL","A","Gsw","iWUE","A at sat light","Gs at sat light","Ci at sat light","iWUE at sat light",
             "Gs at sat light sat CO2","Stomatal limitation",
             "NPQ induction amplitude","NPQ induction rate","NPQ relaxation amplitude","NPQ relaxation rate","NPQ relaxation offset",
             "ΦPSII recovery amplitude","ΦPSII recovery rate","ΦPSII recovery offset","NPQ induction slope",
             "Maximum NPQ","Final NPQ","Final ΦPSII","Fv/Fm","VcmaxT_P25", "VpmaxT_P25", "JmaxT_P25")

varclass <- c(rep("Structural & biochemical",8),
              rep("Gas exchange",13),
              rep("Chlorophyll fluorescence",13))
