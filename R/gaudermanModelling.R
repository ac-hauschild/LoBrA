#' Create the Gouderman Data and Model
#' 
#' The screening method extracts peaks that are likely to 
#' @param s A numeric value to 
#' @return square of the in put
#' @export Use this function 
#' @examples \dontrun{
#' } 
#'   data(modelSelectionObjectExample)
#'   
#'   optimal<-getOptimalSpline(lobraModelSelectionObject, qualityMeasure="AIC")
#'   breaks<- optimal@breaks$optimalID
#'   center<- 14
#'   timeperiod <- 4;
#'   gaudermanLDOexample4TP <- createGoudermanData(modelSelectionObject@ldo, breaks, center, timeperiod)
#'   save(gaudermanLDOexample4TP,file="data/mygaudermanLDOExample4.RData")
#'   
#'   range <- c(9,13)
#'   center<- 9;
#'   gaudermanLDOexample3TP <- createGoudermanData(modelSelectionObject@ldo, breaks=range, center, range=range);
#'   save(gaudermanLDOexample3TP,file="data/mygaudermanLDOExample3.RData")
#'   
createGoudermanData<-function(selectedLDO, breaks, center, timeperiod=NA, range=NA){
  
  
  ## Checking timeperiod and range values.
  times<-selectedLDO@times
  bs<-c(min(times), breaks, max(times))
  if(!is.na(timeperiod)){
    range<-c(bs[timeperiod], bs[timeperiod+1])
  }else if( length(!is.na(range))>1 ){
    if(length(range)!=2){
      stop("Range must be a vector with two entries!")
    }
    range<- sort(range)
    i<- which(bs==range[1])
    j<- which(bs==range[2])
    # i<- which(bs==4)
    if(length(i)==1 & length(j)==1){
      timeperiod<-i;
    }else{
      stop("Each value in range must included once and only once in the breaks vector!")
    }
  }else{
    stop("Function requires either timeperiod or range to be defined!")
  }
  
  if(!(center>=range[1] & center<=range[2])){
    stop("The center variable must be within the defined range!")
  }
  
  
  
  ## Read required parameters.
  sampleIds<-selectedLDO@ids
  peaknames<-selectedLDO@peaknames
  classes<- selectedLDO@labels
  

  ##Caluculate Gauderman range and the new beginning and end of this range. If one of them is equal to the center, it will be set to 0, the remaining on will be either -1 or 1.
  gaudermanRange<- abs(range[1] - range[2]);

  ## Determine the new breaks for the spline.
  k<-c()
  for(i in 1:length(breaks)){
    b<-breaks[i]
    b<- (b-center)/gaudermanRange
    k<-c(k, b)
  }
  names(k)<-1:length(k)
  
  # Create new Gouderman data frames
  dataFrames<-list();
  newtimes<-c()
  # selectedLDO@times
  p<- peaknames[1]
  for(p in peaknames){
    peakmatrix<-selectedLDO@dataMatrices[[p]]
    
    myDataFrame<- getGaudermanDataFrame(peakmatrix, sampleIds, classes, center, timeperiod, gaudermanRange, k)
    dataFrames[[p]]<-myDataFrame;
    newtimes<-c(  newtimes, myDataFrame[,"time"])
    
  }
  newTimeVars<-setdiff(colnames(myDataFrame), colnames(peakmatrix))
  
  newtimes<- as.matrix(unique(newtimes))
  
  myGaudermanLDO<- new("GaudermanLDO",
                       name=paste0("Gauderman-",selectedLDO@name), 
                       dataFrames=dataFrames, 
                       peaknames=peaknames, 
                       k= k,
                       times=newtimes, 
                       newTimeVars=newTimeVars,
                       ids=sampleIds, 
                       labels=classes)
  return(myGaudermanLDO)
}





######## Create Peak Matrices for LMER Model with parameterized Times	#####################################
#############################################################################################################
# k1=-200; k2=0; center = 0
# k1 = 100; k2=300; center = 300


getGeneralizedGaudermanDataFrame=function(peakmatrix, sampleIds, classes, center, timeperiod, gaudermanRange, k){
  testrange<-1:15
  # timeperiod<-4
  
  ## Scale the old times to the new Gauderman range, standardizing the range between rangeStart and rangeEnd to 1.
  oldMins<- as.numeric(peakmatrix[,'time']);
  # head(oldMins);
  newMins = (oldMins - center)/gaudermanRange;
  # head(mins);
  
  
  ## Define the new Gouderman time variables for (1) the time points before k1 (2) the time points between k1 and k2 (k1>=t<=k2) which is the time of interest and (3) the time points larger than k2.
  
  numberS<-length(k)+1
  
  timematrix<-matrix(rep(newMins, numberS), ncol = numberS, nrow = length(newMins))
  
  k0<- c(min(newMins), k, max(newMins))
  names(k0)<-1:length(k0)
  
  i<-4
  for(i in 1:numberS){
    t<- timematrix[,i]
    timematrix[t<=k0[i],i]<- k0[i]
    timematrix[t>k0[i+1],i]<- k0[i+1]
  }

  ## Create a data.frame with all previous fields.
  peakFrame <- data.frame(class=factor(peakmatrix[,"class"], levels=classes), 
                          id = factor(peakmatrix[,"id"], levels=sampleIds), 
                          time=as.numeric(mins), 
                          value=as.numeric(peakmatrix[,"value"]))

  ## Remove k[i] for each column of the time-matrix that is not equal to the specified time period according to the generalized Gauderman algorithm.
  ## Add each column to the data.frame.
  i=1;
  j=1;
  while(i <=numberS){
    # print(paste("i=",i," j=",j))
    t<-as.numeric(timematrix[,i])
    n<-paste0("timep",i)
    if(i==timeperiod){
      # print("break")
      peakFrame[,n]<-t
      i=i+1
      next;
    }
    peakFrame[,n] <- t-as.numeric(k[j])
    j=j+1;
    i=i+1;
  }
  
  ## Return the new data frame for this peak.
  return(peakFrame);
}









#' @title Create the Gouderman Model for each of the metabolites
#' 
#' @description Fitting the Gouderman Model with using Gouderman-Data Arangement.
#' @param mygaudermanLDO GaudermanLDO data object, created by the generalized gauderman algorithm (GGA).
#' @return GaudermanModelEvaluation Results of the evaluation of the Fitted linear mixed effect models for the defined time periods.
#' @export Use this function 
#' @examples \dontrun{
#' } 
#'   data(mygaudermanLDOExample3)
#'   evaluationresult1<- modelGoudermanLongitudinal(gaudermanLDOexample3TP)
#'   
#'   data(mygaudermanLDOExample4)
#'   evaluationresult1<- modelGoudermanLongitudinal(gaudermanLDOexample4TP)
#'   

modelGoudermanLongitudinal<-function(mygaudermanLDO, splinetype="linear"){

  ## Get required parameter from mygaudermanLDO. 
  peaknames<- mygaudermanLDO@peaknames
  inter.knot2<- mygaudermanLDO@newTimeVars
  k<- mygaudermanLDO@k;
  classes<- mygaudermanLDO@labels
  
  
  ## Arrays and Lists to gather the models, p-values that go into the GaudermanModelEvaluation.
  ctrl <- lmeControl(opt='optim');
  pValues<- c();
  modelparameter<- c();
  modellist<- list();
  allIntercepts <- c()
  allSlopes<-c()
  
  ## Run over all selected peaks and fit the specified linear mixed effect model.
  # length(peaknames)
  # currentP <- peaknames[50]
  for(currentP in peaknames){
    # print(currentP);
    ## Get the generalized gauderman data frame.
    mat<-mygaudermanLDO@dataFrames[[currentP]]
    head(mat, 10)

    ## Fit the linear mixed effect model.
    fixedf <- as.formula(paste(c("value ~ class ", inter.knot2 , paste("class:", inter.knot2, sep="")) ,collapse =" + "))
    model.spline <- lme(fixed = fixedf , random = ~ 1+time | id, data = mat, method = "ML", na.action = na.exclude, control=ctrl);
    modellist[[currentP]]<-model.spline;
    
    ## Aquiring pvalues for the Class intercept comparison and the slops for all modeled time periods.
    classname<-paste0("class",classes[2])
    pvaluesOfinterest<- c( classname, paste0(classname, ":", inter.knot2))
    pv<- summary(model.spline)$tTable[pvaluesOfinterest,"p-value"];
    pValues<- rbind(pValues, pv);
    
    ## Aquiring the model parameter.
    modelparameter<- rbind(modelparameter, fixef(model.spline));

  }
  
  correctedpValues<- c();
  mycol<-colnames(pValues)[1]
  for(mycol in colnames(pValues)){
    corrected<- p.adjust(pValues[,mycol], method = "fdr");
    correctedpValues<- cbind(correctedpValues,corrected);
  }
  rownames(pValues) <- peaknames;
  rownames(correctedpValues) <- peaknames;
  rownames(modelparameter) <- peaknames;
  
  
  
  gModelEvaluation<- new("GaudermanModelEvaluation",
                       name=paste0("Evaluation-",mygaudermanLDO@name), 
                       gaudermanLDO=mygaudermanLDO, 
                       models=modellist, 
                       pvalues=pValues,
                       correctedpvalues=correctedpValues,
                       modelparameter=modelparameter)
  
  return(gModelEvaluation) 
}
