


#' @title Evaluation of different spline variants.
#' 
#' @description The model selection method evaluates which spline models achieve the best quality among all tested metabolites. 
#' @param selectedLDO LDO containing all selected metabolites to be used for the model selection.
#' @param potentialSplines Vector of all possible knots to be used for the spline modeling.
#' @param nknots Vector of number of spline knots to be used.
#' @param splinetype spline type default is "linear".
#' @param qualityMeasure Vector of quality measures to be used. Possible options are "AIC", "BIC", and "logLik".
#' @return square of the in put
#' @export Use this function 
#' @examples \dontrun{
#' } 
#'   data(ldoExample)
#'   potentialSplines <- c(6,8,10,12,14,16)
#'   
#'   nknots=c(0,1, 2, 3, 4)
#'   splinetype="linear"
#'   qualityMeasure=c("AIC", "BIC", "logLik")
#'   components <- ldos@selectedPeaks[,"bf"]
#'   components <- names(components)[components]
#'   selectedLDO <- selectComponents(ldo, components)
#'   lobraModelSelectionObject<- lobraModelSelection(selectedLDO, potentialSplines, nknots, splinetype, qualityMeasure)
#'   # save(modelSelectionObject,file="data/modelSelectionObjectExample.RData")
#'   
#'   
lobraModelSelection<-function(selectedLDO, potentialSplines=c(), nknots=c(0,1, 2), splinetype="linear", qualityMeasure=c("AIC", "BIC", "logLik")){
  
  library(nlme)	
  
  modelList= list()
  quality= matrix(0,0,0)
  ctrl <- lmeControl(opt='optim');
  AICTable<-matrix(0,nrow=0,ncol=length(selectedLDO@peaknames));
  colnames(AICTable)<-selectedLDO@peaknames;
  BICTable<-AICTable;
  logLikTable<-AICTable;
  modelList<-list();
  breaks<-list();
  
  currentP <-selectedLDO@peaknames[1]
  if(0 %in% nknots || length(potentialSplines)<1){
    AICS<-c();
    BICS<-c();
    logLikS<-c();
    mlist<-list();
    for(currentP in selectedLDO@peaknames){
      print(currentP);
      mat<-selectedLDO@dataMatrices[[currentP]]
      model <- lme(fixed = value ~ class + time + class:time, random = ~ 1+time | id, data = mat, method = "ML", na.action = na.exclude, control=ctrl);
      mlist[[currentP]]<-model;
      AICS<- c(AICS, summary(model)$AIC);
      BICS<-c(BICS, summary(model)$BIC);
      logLikS<-c(logLikS, summary(model)$logLik);
    }
    
    AICTable<-rbind(AICTable, AICS);
    BICTable<-rbind(BICTable, BICS);
    logLikTable<-rbind(logLikTable, logLikS);
    rownames(AICTable)<-rownames(BICTable)<-rownames(logLikTable)<-"No Spline"
    breaks[["No Spline"]] <-c();
    modelList[["No Spline"]]<-    mlist
  }

  
  ################################# Single knot spline models ##############################################################################################
  if(1 %in% nknots && length(potentialSplines)>=1){
    n<-rownames(AICTable);
    t<-8
    for(t in potentialSplines){
      print(paste("time: ",t));
      AICS<-c();
      BICS<-c();
      logLikS<-c();
      mlist<-list();
      
      for(currentP in selectedLDO@peaknames){
        # print(currentP);
        mat<-selectedLDO@dataMatrices[[currentP]]
        inter.knot <- c(t)
        names(inter.knot) <- c("TimeV")
        
        
        E2<- outer(mat$time, inter.knot,"-")
        ls.mat <- E2*(E2>0)
        mat$timeV<- ls.mat[,1]
        
        model.spline <- lme(fixed = value ~ class + time + timeV + class:time + class:timeV , random = ~ 1+time | id, data = mat, method = "ML", na.action = na.exclude, control=ctrl);
        mlist[[currentP]]<-model.spline;
        
        
        AICS<- c(AICS, summary(model.spline)$AIC);
        BICS<-c(BICS, summary(model.spline)$BIC);
        logLikS<-c(logLikS, summary(model.spline)$logLik);
      }
      
      AICTable<-rbind(AICTable, AICS);
      BICTable<-rbind(BICTable, BICS);
      logLikTable<-rbind(logLikTable, logLikS);
      breaks[[paste("T-",t, sep="")]] <-c(t);
      modelList[[paste("T-",t, sep="")]]<-    mlist
      
    }
    n<- c(n, paste("T-",potentialSplines, sep=""))
    rownames(AICTable)<-rownames(BICTable)<-rownames(logLikTable)<-n
    
  }
  
  
  ################################# Two know spline #################################################################################################
  nknots<- nknots[nknots>1]
  if(length(nknots)>0 && length(potentialSplines)>= min(nknots)){
    
    pS<-powerSet(1:length(potentialSplines))
    pSlengths<-lapply(pS, length)
    k<-2
    for(k in nknots){
      n<-rownames(AICTable);
      print(paste("PowerSets of size:", k));
      idsSets<- pS[pSlengths==k];
      ids<-idsSets[[4]]
      for(ids in idsSets){
        inter.knot <- potentialSplines[ids]
        names(inter.knot) <- paste("time",inter.knot, sep = "")
        tcomb<-paste("time",inter.knot, collapse = "-", sep = "")
        breaks[[tcomb]] <-c(inter.knot);

        AICS<-c();
        BICS<-c();
        logLikS<-c();
        mlist<-list();
        
        
        currentP <- selectedLDO@peaknames[[4]]
        for(currentP in selectedLDO@peaknames){
          # print(currentP);
          mat<-selectedLDO@dataMatrices[[currentP]]
          head(mat,10)
          E2<- outer(mat$time, inter.knot,"-")
          ls.mat <- E2*(E2>0)
          colnames(ls.mat)<- names(inter.knot)
          
          for(t in names(inter.knot)){
            mat[t]<- as.vector(ls.mat[,t])
          }
          head(mat,10)
          
          fixedf <- as.formula(paste(c("value ~ class + time + class:time", names(inter.knot) , paste("class:", names(inter.knot), sep="")) ,collapse =" + "))
          model.spline <- lme(fixed = fixedf , random = ~ 1+time | id, data = mat, method = "ML", na.action = na.exclude, control=ctrl);
          mlist[[currentP]]<-model.spline;
          
          AICS<- c(AICS, summary(model.spline)$AIC);
          BICS<-c(BICS, summary(model.spline)$BIC);
          logLikS<-c(logLikS, summary(model.spline)$logLik);
        }
        
        AICTable<-rbind(AICTable, AICS);
        BICTable<-rbind(BICTable, BICS);
        logLikTable<-rbind(logLikTable, logLikS);
        modelList[[tcomb]]<-    mlist
        n<- c(n, tcomb)
      }
      rownames(AICTable)<-rownames(BICTable)<-rownames(logLikTable)<-n
    }
  }
  quality=list()
  if("AIC" %in% qualityMeasure){
    quality[["AIC"]]<- AICTable;
  }
  if("BIC" %in% qualityMeasure){
    quality[["BIC"]]<- BICTable;
  }
  if("logLik" %in% qualityMeasure){
    quality[["logLik"]]<- logLikTable;
  }
  
  modelSelectionObject<- new("LDOmodelselection", ldo=ldo, 
                             potentialSlines=potentialSplines, 
                             splinetype=splinetype, 
                             qualityMeasure=qualityMeasure, 
                             modelList=modelList,
                             quality=quality,
                             breaks=breaks)
  return(modelSelectionObject)
}


#' @title Plotting background screening results.
#' 
#' @description Plotting the results of the background screening of backbround or confounding components.
#'  
#' @param modelSelectionObject Object of type LDOmodelselection that was created during the model evaluation. @seealso lobraModelSelection
#' @param qualityMeasure Quality measure to be visualized. 
#' @export Use this function 
#' @examples \dontrun{
#' } 
#'   
#' 
plot.modelSelectionEvaluation<-function(modelSelectionObject, qualityMeasure){
  q<-qualityMeasure[1]
  for(q in qualityMeasure){
    data<-modelSelectionObject@quality[[q]]
    
    par( mar = c(5, 12, 4, 6)+0.1)
    vals<-sort(apply(data, 1, FUN=median))
    boxplot(t(data)[,rev(names(vals))], horizontal=TRUE, col="lightblue", ylab='Spline Type', xlab=q, las=2,  main=paste("", q, " spline model comparison."))
    axis(4,at=1:length(rownames(data)),adj=1,labels=rev(round(vals)), las=2)
  }
}


#' @title Extract the optimal spline model parameters from the ModelSelection Object.
#' 
#'  @description The method calculates which spline model and parameters worked best over all different components.
#' @param modelSelectionObject
#' @param qualityMeasure Quality measure to be used to select the optimal spline.
#' @return modelSelectionObject that contains the optimal model.
#' @export Use this function 
#' @examples \dontrun{
#' } 
#'   data(modelSelectionObjectExample)
#'   
#'   optimal<-getOptimalSpline(modelSelectionObject, qualityMeasure="AIC")
#' 
getOptimalSpline<-function(modelSelectionObject, qualityMeasure="AIC"){

  
  meds<-apply(modelSelectionObject@quality[[qualityMeasure]], 1, median)
  which(meds==min(meds))
  optimalID<-names(which(meds==min(meds)))
  optimalModel<- modelSelectionObject@modelList[[optimalID]]
  length(optimalModel)
  optimalSpline<- modelSelectionObject@breaks[[optimalID]]
  
  
  optimalModelObject<- new("LDOmodelselection", ldo=modelSelectionObject@ldo, 
                             potentialSlines=optimalSpline, 
                             splinetype=modelSelectionObject@splinetype, 
                             qualityMeasure=qualityMeasure, 
                             modelList=optimalModel,
                             quality=list(modelSelectionObject@quality[[qualityMeasure]][optimalID,] ),
                             breaks=list(optimalID=optimalSpline))
  
  return(optimalModelObject)
}


#' @title Creating the power set of a set..
#' 
#' @param set Set of numbers or strings.
#' @return Returns powerset of the given set.
powerSet <- function(set) { 
  n <- length(set)
  masks <- 2^(1:n-1)
  lapply( 1:2^n-1, function(u) set[ bitwAnd(u, masks) != 0 ] )
}
