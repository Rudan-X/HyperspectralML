plot_main_goodfit=function(methods_comb,species, year,dataset,legendp,datatype){
  methods=methods_comb[[1]]
  tunings=methods_comb[[2]]
  scenarios<-c("Seen genotype","Unseen genotype","Unseen year")
  score_met<-c("\n(R²)","\n(R²)","\n(Correlation)")
  methodslabels<-c("PLSR","SVM")
  trait_types<-c("A. Structural & biochemical traits", "B. Gas exchange measurements", "C. Chlorophyll fluorescence")
  # trait_types=c('SLA, C%, N% and isotopes','Fitted parameters using A-Ci response curve','Fitted parameters using chlorophyll fluorescence','Light saturated gas exchange rates')

  maize_vars<-list()
  maize_vars[[1]]<-c("SLA","N","CN")
  maize_vars[[2]]<-c("Vpmax","Vmax","Vpmax.Vmax","a400","A_sat")
  maize_vars[[3]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  maize_vars[[4]]<-c("SLA")
  maize_vars[[5]]<-c("A_sat")
  maize_vars[[6]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

  barley_vars<-list()
  barley_vars[[1]]<-c("SLA","C","N","CN" )
  barley_vars[[2]]<-c("Vcmax","Jmax","TPU","asat","a400")#,"A_sat","gs_sat","ci_sat","iWUE_sat"
  barley_vars[[3]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")
  barley_vars[[4]]<-c("SLA")
  barley_vars[[5]]<-c("")
  barley_vars[[6]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")

  var1<-c("SLA","C","N","CN","Vcmax","Jmax","TPU","asat","Vpmax","Vmax","Vpmax.Vmax","a400","A_sat","gs_sat","ci_sat","iWUE_sat","NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  var2<-c("SLA","%C","%N","C/N","Vcmax","Jmax","TPU","A at\n sat light\n sat CO2","Vpmax","Vmax","Vpmax/Vmax","A at\n sat light\n 400 CO2","A at\n saturating\n light","A at\n saturating\n light","Ci at\n saturating\n light","iWUE at\n saturating\n light","NPQ\n relaxation\n residual","phiPSII\n induction\n residual","Maximum\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")


  if (species=="maize"){
    variables<-maize_vars
    stopsce<-3
  }else{
    variables<-barley_vars
    stopsce<-2
  }
  for (trait_ind in 1:3){
     for (scena in 1:stopsce){
        if (scena==3){
          starti<-3
        }else{starti<-0}
      varnames<-variables[[starti+trait_ind]]
      nvar<-length(varnames)

      for (t in 1:length(methods)){
        method=methods[t]
        tuning=tunings[t]
        if (scena==3){
          R2mat2fold=matrix(0,1,nvar)
          R2mat=matrix(0,1,nvar)
        }else{
          R2mat2fold=matrix(0,10,nvar)
          R2mat=matrix(0,6,nvar)
        }
        for (i in 1:nvar){
          varname=varnames[i]

          if (scena==3){
            filen <- paste("results/scenario", scena, "/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
            data=read.csv(file = filen)
            R2mat[1,i]=round(cor(data$meas,data$pred),4)

            filen <- paste("results/scenario", scena, "_agg/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
            data=read.csv(file = filen)
            R2mat2fold[1,i]=round(cor(data$meas,data$pred),4)
          }else{
            for (f in 1:6){
              filen <- paste("results/scenario", scena, "/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_testfold", f, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat[f,i]=round(get_R2(data$meas,data$pred),4)
            }
            c=1
            for (rep in 1:5){
              for (p in 0:1){
                filen <- paste("results/scenario", scena, "_2fold/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_rep", rep,"_iter",p, ".txt", sep = "")
                data=read.csv(filen)
                R2mat2fold[c,i]=round(get_R2(data$meas,data$pred),4)
                c=c+1
              }
            }
          }
        }

        newvar<-match(varnames,var1)

        colnames(R2mat)<-var2[newvar]
        df=melt(R2mat)
        df$Method=methodslabels[t]

        colnames(R2mat2fold)<-var2[newvar]
        df2fold=melt(R2mat2fold)
        df2fold$Method=methodslabels[t]
        if (t==1){
          dfscore=df
          dfscore2fold=df2fold
        }else{
          dfscore=rbind(dfscore,df)
          dfscore2fold=rbind(dfscore2fold,df2fold)
        }
      }

      # provideDimnames(R2mat) # returns for each fold, method versus variables
      colnames(dfscore)=c("TestedFold","Variable","Score","Method")
      dfscore$Scenario=paste0(scenarios[scena],score_met[scena])
      colnames(dfscore2fold)=c("TestedFold","Variable","Score","Method")
      dfscore2fold$Scenario=paste0(scenarios[scena],score_met[scena])

      if (scena==1){
        totaldfscore=dfscore
        totaldfscore2fold=dfscore2fold
      }else{
        totaldfscore=rbind(totaldfscore,dfscore)
        totaldfscore2fold=rbind(totaldfscore2fold,dfscore2fold)
      }
    }
    totaldfscore$Traittype<-trait_types[trait_ind]
    totaldfscore2fold$Traittype<-trait_types[trait_ind]

    if (trait_ind==1){
      finaldf<-totaldfscore
      finaldf2fold<-totaldfscore2fold
    }else{
      finaldf<-rbind(finaldf,totaldfscore)
      finaldf2fold<-rbind(finaldf2fold,totaldfscore2fold)
    }

  }
  if (datatype=="Original"){
    datatoplot<-finaldf
  }else if (datatype=="Averaged"){
    datatoplot<-finaldf2fold
  }

  p <- ggplot(datatoplot, aes(x=Variable, y=Score,fill=Method)) +
    geom_boxplot()  + #+ scale_x_discrete(labels= varnames)
    theme(legend.position = legendp,axis.text.x=element_text(angle=0),text = element_text(size = 18),axis.title.y = element_blank(),axis.title.x = element_blank()) + #
    facet_grid(Scenario ~ Traittype,scales="free", switch = "y") +
    scale_y_continuous(position = "right" )+ #, limits = c(-1, 0.9)
    labs(y="Coefficient of determination")

# scale_fill_discrete(labels = paste0(methodslab,"_",tunings)) +
  return(p)
}



plot_scatter_goodfit=function(methods_comb,species,dataset,legendp,datatype){
  methods=methods_comb[[1]]
  tunings=methods_comb[[2]]
  scenarios<-c("Seen genotype","Unseen genotype","Unseen year")
  score_met<-c("\n(R²)","\n(R²)","\n(Correlation)")
  methodslabels<-c("PLSR","SVM")
  trait_types<-c("A. Structural & biochemical traits", "B. Gas exchange measurements", "C. Chlorophyll fluorescence")
  # trait_types=c('SLA, C%, N% and isotopes','Fitted parameters using A-Ci response curve','Fitted parameters using chlorophyll fluorescence','Light saturated gas exchange rates')

  maize_vars<-list()
  maize_vars[[1]]<-c("SLA","N","CN")
  maize_vars[[2]]<-c("Vpmax","Vmax","Vpmax.Vmax","a400","A_sat")
  maize_vars[[3]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  maize_vars[[4]]<-c("SLA")
  maize_vars[[5]]<-c("A_sat")
  maize_vars[[6]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

  barley_vars<-list()
  barley_vars[[1]]<-c("SLA","C","N","CN" )
  barley_vars[[2]]<-c("Vcmax","Jmax","TPU","asat","a400")#,"A_sat","gs_sat","ci_sat","iWUE_sat"
  barley_vars[[3]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")
  barley_vars[[4]]<-c("SLA")
  barley_vars[[5]]<-c("")
  barley_vars[[6]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")

  var1<-c("SLA","C","N","CN","Vcmax","Jmax","TPU","asat","Vpmax","Vmax","Vpmax.Vmax","a400","A_sat","gs_sat","ci_sat","iWUE_sat","NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  var2<-c("SLA","%C","%N","C/N","Vcmax","Jmax","TPU","A at\n sat light\n sat CO2","Vpmax","Vmax","Vpmax/Vmax","A at\n sat light\n 400 CO2","A at saturating light","A at\n saturating\n light","Ci at\n saturating\n light","iWUE at\n saturating\n light","NPQ\n relaxation\n residual","phiPSII induction residual","Maximum NPQ","Final Fv/Fm","Initial Fv/Fm")


  if (species=="maize"){
    variables<-maize_vars
    stopsce<-3
  }else{
    variables<-barley_vars
    stopsce<-2
  }

  for (year in c(2021,2022)){
    for (t in 1:length(methods)){
      for (trait_ind in 1:3){
        scena <-3
        starti<-3

        varnames<-variables[[starti+trait_ind]]
        nvar<-length(varnames)

        method=methods[t]
        tuning=tunings[t]
        for (i in 1:nvar){
          varname=varnames[i]

          if (datatype=="Original"){
            filen <- paste("results/scenario", scena, "/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
          }else if (datatype=="Averaged"){
            filen <- paste("results/scenario", scena, "_agg/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
          }

          data=read.csv(file = filen)

          newvar<-match(varname,var1)
          data$Variable<-var2[newvar]
          data$Method<-methodslabels[t]
          data$Year<-year

          if (trait_ind==1 & i==1){
            df<-data
          }else{
            df<-rbind(df,data)
          }
        }
      }
      if (t==1){
        dfscore<-df
      }else{
        dfscore<-rbind(dfscore,df)
      }
    }
    if (year==2021){
      dftot<-dfscore
    }else{dftot<-rbind(dftot,dfscore)}

  }

  p<-ggplot(dftot, aes(x=meas, y=pred,shape=Method,color=Method)) +
    geom_point()  + #+ scale_x_discrete(labels= varnames)
    theme(legend.position = legendp,axis.text.x=element_text(angle=0),text = element_text(size = 18),axis.title.y = element_blank(),axis.title.x = element_blank()) + #
    facet_wrap(Year~Variable ,scales="free",ncol = 6) +
    stat_cor(method = "pearson")+
#     labs(x="Predicted",y="Measred")+
    scale_y_continuous(position = "right" ) #, limits = c(-1, 0.9)


# scale_fill_discrete(labels = paste0(methodslab,"_",tunings)) +
  return(p)
}



plot_main_robustness=function(methods_comb,dataset){
  methods=methods_comb[[1]]
  tunings=methods_comb[[2]]
  methodslab<-c("PLSR","SVM")
  scenarios<-c("Seen genotype","Unseen genotype","Unseen year")
  score_met<-c("\n(R²)","\n(R²)","\n(Correlation)")
  trait_types<-c("A. Structural & biochemical traits", "B. Gas exchange measurements", "C. Chlorophyll fluorescence")

  maize_vars<-list()
  maize_vars[[1]]<-c("SLA","N","CN")
  maize_vars[[2]]<-c("Vpmax","Vmax","a400","A_sat")
  maize_vars[[3]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  maize_vars[[4]]<-c("SLA")
  maize_vars[[5]]<-c("A_sat")
  maize_vars[[6]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

  barley_vars<-list()
  barley_vars[[1]]<-c("SLA","C","N","CN" )
  barley_vars[[2]]<-c("Vcmax","Jmax","TPU","asat","a400")#,"A_sat","gs_sat","ci_sat","iWUE_sat"
  barley_vars[[3]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")
  barley_vars[[4]]<-c("SLA")
  barley_vars[[5]]<-c("")
  barley_vars[[6]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")

  var1<-c("SLA","C","N","CN","Vcmax","Jmax","TPU","asat","Vpmax","Vmax","Vpmax.Vmax","a400","A_sat","gs_sat","ci_sat","iWUE_sat","NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  var2<-c("SLA","%C","%N","C/N","Vcmax","Jmax","TPU","A at\n sat light\n sat CO2","Vpmax","Vmax","Vpmax/Vmax","A at\n sat light\n 400 CO2","A at\n saturating\n light","A at\n saturating\n light","Ci at\n saturating\n light","iWUE at\n saturating\n light","NPQ\n relaxation\n residual","phiPSII\n induction\n residual","Maximum\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")


  Spec<-c("A. Maize (2021)","B. Barley (2022)")
  s<-1
  for (species in c("maize","barley")){
    if (species=="maize"){
      variables<-maize_vars
      stopsce<-3
    }else{
      variables<-barley_vars
      stopsce<-2
    }
    if (species=="maize"){year<-2021}else{year<-2022}
    for (trait_ind in 1:2){
      for (scena in 1:2){
        if (scena==3){
          starti<-3
        }else{starti<-0}
        varnames<-variables[[starti+trait_ind]]
        nvar<-length(varnames)

        for (t in 1:length(methods)){
          method=methods[t]
          tuning=tunings[t]

          if (scena==3){
            R2mat2fold=matrix(0,1,nvar)
            R2mat6fold=matrix(0,1,nvar)
          }else{
            R2mat2fold=matrix(0,10,nvar)
            R2mat6fold=matrix(0,6,nvar)
          }
          for (i in 1:nvar){
            varname=varnames[i]

            if (scena==3){
              filen <- paste("results/scenario", scena, "_agg/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat2fold[1,i]=round(cor(data$meas,data$pred),4)

              filen <- paste("results/scenario", scena, "/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat6fold[1,i]=round(cor(data$meas,data$pred),4)
            }else{
              for (f in 1:6){
              filen <- paste("results/scenario", scena, "/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_testfold", f, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat6fold[f,i]=round(get_R2(data$meas,data$pred),4)
              }

              c=1
              for (rep in 1:5){
                for (p in 0:1){
                  filen <- paste("results/scenario", scena, "_2fold/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_rep", rep,"_iter",p, ".txt", sep = "")
                  data=read.csv(filen)
                  R2mat2fold[c,i]=round(get_R2(data$meas,data$pred),4)
                  c=c+1
                }
              }
            }
          }
          newvar<-match(varnames,var1)

          colnames(R2mat2fold)<-var2[newvar]
          colnames(R2mat6fold)<-var2[newvar]
          df2=melt(R2mat2fold)
          df6=melt(R2mat6fold)

          df2$Method_Data=paste0(methodslab[t],"_Averaged")
          df6$Method_Data=paste0(methodslab[t],"_Original")

          if (t==1){
            dfscore2=df2
            dfscore6=df6
          }else{
            dfscore2=rbind(dfscore2,df2)
            dfscore6=rbind(dfscore6,df6)
          }
        }

      # provideDimnames(R2mat) # returns for each fold, method versus variables
      colnames(dfscore2)=c("TestedFold","Variable","Score","Method_Data")
      colnames(dfscore6)=c("TestedFold","Variable","Score","Method_Data")

      dfscore2$Scenario=paste0(scenarios[scena], score_met[scena])
      dfscore6$Scenario=paste0(scenarios[scena], score_met[scena])

      if (scena==1){
        totaldfscore2=dfscore2
        totaldfscore6=dfscore6
      }else{
        totaldfscore2=rbind(totaldfscore2,dfscore2)
        totaldfscore6=rbind(totaldfscore6,dfscore6)
      }
    }
    totaldfscore2$Traittype<-trait_types[trait_ind]
    totaldfscore6$Traittype<-trait_types[trait_ind]
    if (trait_ind==1){
      finaldf6fold<-totaldfscore6
      finaldf2fold<-totaldfscore2
    }else{
      finaldf6fold<-rbind(finaldf6fold,totaldfscore6)
      finaldf2fold<-rbind(finaldf2fold,totaldfscore2)
    }
    }
    finaldf6fold$Species<-Spec[s]
    finaldf2fold$Species<-Spec[s]
    if (species=="maize"){
      finaldf6fold0<-finaldf6fold
      finaldf2fold0<-finaldf2fold
    }else{
      finaldf6fold0<-rbind(finaldf6fold0,finaldf6fold)
      finaldf2fold0<-rbind(finaldf2fold0,finaldf2fold)
    }
    s<-s+1
  }

  final<-rbind(finaldf6fold0,finaldf2fold0)

  p <- ggplot(final, aes(x=Variable, y=Score,fill=Method_Data)) +
    geom_boxplot()  + #+ scale_x_discrete(labels= varnames)
    theme(legend.position = c(.05, .1),text = element_text(size = 18),axis.title.y = element_blank(),axis.title.x = element_blank()) +
    facet_grid(Scenario ~ Species,scales="free", switch = "y") + #Scenario ~ Traittype
    scale_y_continuous(position = "right")


  return(p)
}

plot_supp_MLmethods<-function(methods_comb,dataset){
  methods=methods_comb[[1]]
  tunings=methods_comb[[2]]
  methodslab<-c("PLSR","SVM (Linear)","Elastic Net", "Random Forest","SVM (Polyn)")
  scenarios<-c("Seen genotype","Unseen genotype","Unseen year")
  score_met<-c("\n(R²)","\n(R²)","\n(Correlation)")
  trait_types<-c("A. Structural & biochemical traits", "B. Gas exchange measurements", "C. Chlorophyll fluorescence")

  maize_vars<-list()
  maize_vars[[1]]<-c("SLA","N","CN")
  maize_vars[[2]]<-c("Vpmax","Vmax") #,"A_sat","a400"
  maize_vars[[3]]<-c("phiPSII_ind_res","maxNPQ") #,"endFvFm","initialFvFm"
  maize_vars[[4]]<-c("SLA")
  maize_vars[[5]]<-c("A_sat")
  maize_vars[[6]]<-c("phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")

  barley_vars<-list()
  barley_vars[[1]]<-c("SLA","N","CN" )
  barley_vars[[2]]<-c("Vcmax","Jmax","TPU")#,"A_sat","gs_sat","ci_sat","iWUE_sat"
  barley_vars[[3]]<-c("phiPSII_ind_res","maxNPQ") #"NPQ_rel_res", ,"endFvFm"
  barley_vars[[4]]<-c("SLA")
  barley_vars[[5]]<-c("")
  barley_vars[[6]]<-c("NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm")

  var1<-c("SLA","C","N","CN","Vcmax","Jmax","TPU","asat","Vpmax","Vmax","Vpmax.Vmax","a400","A_sat","gs_sat","ci_sat","iWUE_sat","NPQ_rel_res","phiPSII_ind_res","maxNPQ","endFvFm","initialFvFm")
  var2<-c("SLA","%C","%N","C/N","Vcmax","Jmax","TPU","A at\n sat light\n sat CO2","Vpmax","Vmax","Vpmax/Vmax","A at\n sat light\n 400 CO2","A at\n saturating\n light","A at\n saturating\n light","Ci at\n saturating\n light","iWUE at\n saturating\n light","NPQ\n relaxation\n residual","phiPSII\n induction\n residual","Maximum\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")


  Spec<-c("A. Maize (2021)","B. Barley (2022)")
  s<-1
  for (species in c("maize","barley")){ #
    if (species=="maize"){
      variables<-maize_vars
      stopsce<-3
    }else{
      variables<-barley_vars
      stopsce<-2
    }
    if (species=="maize"){year<-2021}else{year<-2022}
    for (trait_ind in 1:3){
      for (scena in 1:1){
        if (scena==3){
          starti<-3
        }else{starti<-0}
        varnames<-variables[[starti+trait_ind]]
        nvar<-length(varnames)

        for (t in 1:5){ #length(methods)
          method=methods[t]
          tuning=tunings[t]

          if (scena==3){
            R2mat2fold=matrix(0,1,nvar)
            R2mat6fold=matrix(0,1,nvar)
          }else{
            R2mat2fold=matrix(0,10,nvar)
            R2mat6fold=matrix(0,6,nvar)
          }
          for (i in 1:nvar){
            varname=varnames[i]

            if (scena==3){
              filen <- paste("results/scenario", scena, "_agg/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
#              data=read.csv(file = filen)
#              R2mat2fold[1,i]=round(cor(data$meas,data$pred),4)

              filen <- paste("results/scenario", scena, "/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat6fold[1,i]=round(cor(data$meas,data$pred),4)
            }else{
              for (f in 1:6){
              filen <- paste("results/scenario", scena, "/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_testfold", f, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat6fold[f,i]=round(get_R2(data$meas,data$pred),4)
              }

              c=1
              for (rep in 1:5){
                for (p in 0:1){
                  filen <- paste("results/scenario", scena, "_2fold/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_rep", rep,"_iter",p, ".txt", sep = "")
#                 data=read.csv(filen)
#                 R2mat2fold[c,i]=round(get_R2(data$meas,data$pred),4)
                  c=c+1
                }
              }
            }
          }
          newvar<-match(varnames,var1)

          colnames(R2mat2fold)<-var2[newvar]
          colnames(R2mat6fold)<-var2[newvar]
          df2=melt(R2mat2fold)
          df6=melt(R2mat6fold)

          df2$Method_Data=methodslab[t] # paste0(methodslab[t],"_Averaged")
          df6$Method_Data=methodslab[t] # paste0(methodslab[t],"_Original")

          if (t==1){
            dfscore2=df2
            dfscore6=df6
          }else{
            dfscore2=rbind(dfscore2,df2)
            dfscore6=rbind(dfscore6,df6)
          }
        }

      # provideDimnames(R2mat) # returns for each fold, method versus variables
      colnames(dfscore2)=c("TestedFold","Variable","Score","Method_Data")
      colnames(dfscore6)=c("TestedFold","Variable","Score","Method_Data")

      dfscore2$Scenario=paste0(scenarios[scena], score_met[scena])
      dfscore6$Scenario=paste0(scenarios[scena], score_met[scena])

      if (scena==1){
        totaldfscore2=dfscore2
        totaldfscore6=dfscore6
      }else{
        totaldfscore2=rbind(totaldfscore2,dfscore2)
        totaldfscore6=rbind(totaldfscore6,dfscore6)
      }
    }
    totaldfscore2$Traittype<-trait_types[trait_ind]
    totaldfscore6$Traittype<-trait_types[trait_ind]
    if (trait_ind==1){
      finaldf6fold<-totaldfscore6
      finaldf2fold<-totaldfscore2
    }else{
      finaldf6fold<-rbind(finaldf6fold,totaldfscore6)
      finaldf2fold<-rbind(finaldf2fold,totaldfscore2)
    }
    }
    finaldf6fold$Species<-Spec[s]
    finaldf2fold$Species<-Spec[s]
    if (species=="maize"){
      finaldf6fold0<-finaldf6fold
      finaldf2fold0<-finaldf2fold
    }else{
      finaldf6fold0<-rbind(finaldf6fold0,finaldf6fold)
      finaldf2fold0<-rbind(finaldf2fold0,finaldf2fold)
    }
    s<-s+1
  }

  final<-rbind(finaldf6fold0,finaldf2fold0)

  p <- ggplot(finaldf6fold0, aes(x=Variable, y=Score,fill=Method_Data)) +
    geom_boxplot()  + #+ scale_x_discrete(labels= varnames)
    theme(legend.position = c(.9, .9),text = element_text(size = 18),axis.title.y = element_blank(),axis.title.x = element_blank()) +
    facet_grid(Species ~ Scenario,scales="free", switch = "y") + #Scenario ~ Traittype
    scale_y_continuous(position = "right")


  return(p)
}


plot_supp_badfit=function(methods_comb,species, year,dataset,legendp,datatype){
  methods=methods_comb[[1]]
  tunings=methods_comb[[2]]
  scenarios<-c("Seen genotype","Unseen genotype","Unseen year")
  score_met<-c("\n(R²)","\n(R²)","\n(Correlation)")
  methodslabels<-c("PLSR","SVM")
  if (year==2021){
    trait_types<-c("A. Structural &\n biochemical traits", "B. Gas exchange measurements", "C. Chlorophyll fluorescence")
  }else{
    trait_types<-c("A. Struct.\n traits", "B. Gas exchange\n measurements", "C. Chlorophyll fluorescence")
  }

  maize_vars<-list()
  maize_vars[[1]]<-c("C","d13C","d15N")
  maize_vars[[2]]<-c("SL","a400","gsw","iWUE","gs_sat","ci_sat","iWUE_sat")
  maize_vars[[3]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","NPQ_ind_linear","endNPQ")
  maize_vars[[4]]<-c("SLA")
  maize_vars[[5]]<-c("A_sat","gs_sat","ci_sat","iWUE_sat")
  maize_vars[[6]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res","phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear","maxNPQ","endNPQ","endFvFm","initialFvFm")

  barley_vars<-list()
  barley_vars[[1]]<-c("Area","d15N","d13C")
  barley_vars[[2]]<-c("gs","iWUE","sl","A_sat","gs_sat","ci_sat","iWUE_sat")
  barley_vars[[3]]<-c("NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","phiPSII_ind_amp","phiPSII_ind_rate","NPQ_ind_linear","endNPQ")

  var1<-c("SLA","Area","C","d13C","d15N","SL","a400","gsw","iWUE","A_sat","gs_sat","ci_sat","iWUE_sat",
          "gs","sl",
          "NPQ_ind_amp","NPQ_ind_rate","NPQ_rel_amp","NPQ_rel_rate","NPQ_rel_res",
          "phiPSII_ind_amp","phiPSII_ind_rate","phiPSII_ind_res","NPQ_ind_linear",
          "maxNPQ","endNPQ","endFvFm","initialFvFm")
  var2<-c("SLA","Leaf area","%C","d13C","d15N","Stomatal\n limitation","A at\n sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","iWUE at\n sat light\n 400 CO2","A at\n saturating\n light","Gs at\n saturating\n light","Ci at\n saturating\n light","iWUE at\n saturating\n light",
          "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
          "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
          "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n induction\n residual","NPQ\n induction\n slope",
          "Maximum\n NPQ","Final\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")
#  var2<-c("SLA","Leaf area","%C","d13C","d15N","Stomatal\n limitation","A at\n sat light\n 400 CO2","Gsw\n sat light\n 400 CO2","A/Gsw(iWUE)\n at light\n 400 CO2","A at\n sat\n light","Gs at\n sat\n light","Ci at\n sat\n light","iWUE at\n sat\n light",
#          "Gs at\n sat light\n sat CO2","Stomatal\n limitation",
#          "NPQ\n induction\n amplitude","NPQ\n induction\n rate","NPQ\n relaxation\n amplitude","NPQ\n relaxation\n rate","NPQ\n relaxation\n residual",
#          "phiPSII\n induction\n amplitude","phiPSII\n induction\n rate","phiPSII\n induction\n residual","NPQ\n induction\n slope",
#          "Maximum\n NPQ","Final\n NPQ","Final\n Fv/Fm","Initial\n Fv/Fm")

  if (species=="maize"){
    variables<-maize_vars
    stopsce<-3
  }else{
    variables<-barley_vars
    stopsce<-2
  }
  for (trait_ind in 1:3){
     for (scena in 1:2){
        if (species=="maize" & year==2022){
          starti<-3
        }else{starti<-0}
      varnames<-variables[[starti+trait_ind]]
      nvar<-length(varnames)

      for (t in 1:length(methods)){
        method=methods[t]
        tuning=tunings[t]
        if (scena==3){
          R2mat2fold=matrix(0,1,nvar)
          R2mat=matrix(0,1,nvar)
        }else{
          R2mat2fold=matrix(0,10,nvar)
          R2mat=matrix(0,6,nvar)
        }
        for (i in 1:nvar){
          varname=varnames[i]

          if (scena==3){
            filen <- paste("results/scenario", scena, "/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
            data=read.csv(file = filen)
            R2mat[1,i]=round(cor(data$meas,data$pred),4)

            filen <- paste("results/scenario", scena, "_agg/whole_data_",dataset, "_", species, "_", method,"_", tuning, "_", varname,"_testingyear",year, ".txt", sep = "")
            data=read.csv(file = filen)
            R2mat2fold[1,i]=round(cor(data$meas,data$pred),4)
          }else{
            for (f in 1:6){
              filen <- paste("results/scenario", scena, "/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_testfold", f, ".txt", sep = "")
              data=read.csv(file = filen)
              R2mat[f,i]=round(get_R2(data$meas,data$pred),4)
            }
            c=1
            for (rep in 1:5){
              for (p in 0:1){
                filen <- paste("results/scenario", scena, "_2fold/",dataset, "_", species, "_", method,"_", tuning, "_", varname, year, "_rep", rep,"_iter",p, ".txt", sep = "")
                data=read.csv(filen)
                R2mat2fold[c,i]=round(get_R2(data$meas,data$pred),4)
                c=c+1
              }
            }
          }
        }

        newvar<-match(varnames,var1)

        colnames(R2mat)<-var2[newvar]
        df=melt(R2mat)
        df$Method=methodslabels[t]

        colnames(R2mat2fold)<-var2[newvar]
        df2fold=melt(R2mat2fold)
        df2fold$Method=methodslabels[t]
        if (t==1){
          dfscore=df
          dfscore2fold=df2fold
        }else{
          dfscore=rbind(dfscore,df)
          dfscore2fold=rbind(dfscore2fold,df2fold)
        }
      }

      # provideDimnames(R2mat) # returns for each fold, method versus variables
      colnames(dfscore)=c("TestedFold","Variable","Score","Method")
      dfscore$Scenario=paste0(scenarios[scena],score_met[scena])
      colnames(dfscore2fold)=c("TestedFold","Variable","Score","Method")
      dfscore2fold$Scenario=paste0(scenarios[scena],score_met[scena])

      if (scena==1){
        totaldfscore=dfscore
        totaldfscore2fold=dfscore2fold
      }else{
        totaldfscore=rbind(totaldfscore,dfscore)
        totaldfscore2fold=rbind(totaldfscore2fold,dfscore2fold)
      }
    }
    totaldfscore$Traittype<-trait_types[trait_ind]
    totaldfscore2fold$Traittype<-trait_types[trait_ind]

    if (trait_ind==1){
      finaldf<-totaldfscore
      finaldf2fold<-totaldfscore2fold
    }else{
      finaldf<-rbind(finaldf,totaldfscore)
      finaldf2fold<-rbind(finaldf2fold,totaldfscore2fold)
    }

  }
  if (datatype=="Original"){
    datatoplot<-finaldf
  }else if (datatype=="Averaged"){
    datatoplot<-finaldf2fold
  }

  if (min(datatoplot$Score) < (-1)){
    p <- ggplot(datatoplot, aes(x=Variable, y=Score,fill=Method)) +
    geom_boxplot()  + #+ scale_x_discrete(labels= varnames)
    theme(legend.position = legendp,axis.text.x=element_text(angle=0),text = element_text(size = 18),axis.title.y = element_blank(),axis.title.x = element_blank()) + #
    facet_grid(Scenario ~ Traittype,scales="free", space="free_x",switch = "y") +
    scale_y_continuous(position = "right", limits=c(-1,0.9))+
    labs(y="Coefficient of determination")+
    theme(axis.text.x=element_text(angle=0,size=12))
  }else{
    p <- ggplot(datatoplot, aes(x=Variable, y=Score,fill=Method)) +
    geom_boxplot()  + #+ scale_x_discrete(labels= varnames)
    theme(legend.position = legendp,axis.text.x=element_text(angle=0),text = element_text(size = 18),axis.title.y = element_blank(),axis.title.x = element_blank()) + #
    facet_grid(Scenario ~ Traittype,scales="free", space="free_x",switch = "y") +
    scale_y_continuous(position = "right")+
    labs(y="Coefficient of determination")+
    theme(axis.text.x=element_text(angle=0,size=12))
  }




# scale_fill_discrete(labels = paste0(methodslab,"_",tunings)) +
  return(p)
}

