#' @title Screening of backbround or confounding components
#' 
#' @description Background noise signals originating from experimental settings or random events can hugely influence the signal pattern of the breath. Background data enables the detailed evaluation and differentiation of the compounds originating primarily from the background or confounding factors as compared to those from the sample itself. The method assumes that all compounds of interest show a larger variation in the sample as compared to the background noise.
#' @param ldo Longitudinal Data Object
#' @param method list of tests to perform, standard values: "bf", "levene" or "bartlett"). bf relates to Brown-Forsythe Levene-type procedure, levene uses classical Levene's procedure and bartlett applies Bartlett's test.
#' @param alpha A numeric value to defining the cutoff to select peaks.
#' @param criteria indicators which criterias to use for screening decision.
#' 
#' @return square of the in put
#' @export Use this function 
#' @return LDO data object
#' @examples \dontrun{
#' } 
#' 
#' data(ldoExample)
#' method= c("bf", "levene", "bartlett")
#' alpha =0.05
#' criteria=c(1,1)
#' ldos<-screening(ldo, method, alpha, criteria)
#' #save(ldos,file="data/ldosExample.RData")
#' 
#' 
screening<-function(ldo, method=c("bf", "levene", "bartlett"), alpha =0.05, criteria=c(1,1)){
  require(lawstat)
  allalphapvals<-c();
  allrespvals<-c();
  
  sAlphasVar<- c();
  bgAlphasVar<- c();
  sResVar<- c();
  bgResVar<- c();
  experimentResiduals<-list();
  experimentIntercept <-list();
  
  currentP<- ldo@peaknames[71]
  for(currentP in ldo@peaknames){
    print(currentP)
    par(mfrow=c(1,2));
    sPeaks<- c()
    ventilatorPeaks <- c()
    rmeans<- c()
    rvars<-c()
    vmeans<- c()
    vvars<-c()
    
    smodel<- lm(value~id, ldo@dataMatrices[[currentP]])
    sAlphas<- smodel$coefficients[-1]
    sRes<- smodel$residuals
    sResVar<- c(sResVar, var(sAlphas));
    sAlphasVar<- c(sAlphasVar, var(sRes));
    
    bmodel<- lm(value~id, ldo@backgroundMatrices[[currentP]])
    bgAlphas<- bmodel$coefficients[-1]
    bgRes<- bmodel$residuals
    bgAlphasVar<- c(bgAlphasVar, var(bgAlphas));
    bgResVar<- c(bgResVar, var(bgRes));
    
    alpy<- c(sAlphas, bgAlphas)
    alpgroup<-factor(c(rep("S", length(sAlphas)), rep("B", length(bgAlphas))))
    resy<- c(sRes, bgRes)
    resgroup<-factor(c(rep("S", length(sRes)), rep("B", length(bgRes))))
    alphapvals<-c();
    respvals<-c();
    
    for(m in method){
      alphapvals<-c( alphapvals, getPvalue(alpy, alpgroup, m));
      respvals<-c( respvals, getPvalue(resy, resgroup, m));
    }
    
    allalphapvals<-rbind(allalphapvals, alphapvals);
    allrespvals<-rbind(allrespvals,respvals);
    experimentIntercept[[currentP]] <- list(sAlphas, bgAlphas);
    experimentResiduals[[currentP]] <- list(sRes, bgRes);
  }
  rownames(allalphapvals)<-ldo@peaknames;
  rownames(allrespvals)<-ldo@peaknames;
  colnames(allalphapvals)<-method;
  colnames(allrespvals)<-method;
  
  
  selectedPeaks<-c(); 
  for(m in method){
    pv1<-allalphapvals[,m];
    pv2<-allrespvals[,m];
    cbind(pv1,pv2);
    tested <- cbind(pv1<alpha, pv2<alpha, bgAlphasVar<sAlphasVar, bgResVar<sResVar);
    if(sum(criteria)>1){
      accepted <- (tested[,1] & tested[,3]) | (tested[,2] & tested[,4]);
    }else if(criteria[1]==1){
      accepted <- (tested[,1]) | (tested[,2]);
    }else if(criteria[2]==1){
      accepted <- (tested[,3]) | (tested[,4]);
    }
    selectedPeaks<-cbind(selectedPeaks, accepted); 
    # Accepted <- rbind(Accepted, apply(Accepted, 2, FUN=sum));
    
  }
  colnames(selectedPeaks)<-method;
    
  
  ldoscreen<-new("LDOscreening",ldo=ldo, experimentIntercept=experimentIntercept, experimentResiduals=experimentResiduals,  interceptPvalues=allalphapvals, residualPvalues=allrespvals, selectedPeaks=selectedPeaks)

  
  
  return(ldoscreen);
}




#' @title Plotting the screening results.
#' 
#' @param ldoscreen LDO screening result
#' @param plotAll Select all components to be plottet. Default plots only the selected peaks using the correctionmethod.
#' @param correctionmethod Version of correction method to be used to select the peaks. Valid values are "bf", "levene", and "bartlett".
#' @param decs decimal numbers of p-values to be plotted.
#' @export Use this function 
#' @return LDO data object
#' @examples \dontrun{
#' } 
#' 
#' data(ldosExample)
#' correctionmethod= "levene"
#' 
#' filename<-paste("screeningresults.pdf",sep="") 
#' pdf(filename, width=12, height=6)
#' plotScreening(ldoscreen)
#' dev.off();
#' 
plotScreening=function(ldoscreen, plotAll=FALSE, correctionmethod="levene", decs=3){
  ldoscreen@ldo@peaknames
  peaknames<-rownames(ldoscreen@selectedPeaks)
  if(!plotAll){
    peaknames<- names(which(ldoscreen@selectedPeaks[,correctionmethod] ))
  }
  
  p<-peaknames[30]
  for(p in peaknames){
    par(mfrow=c(1,2))
    intensities <-ldoscreen@experimentIntercept[[p]]
    pv<-ldoscreen@interceptPvalues[p,]
    names(intensities)<- c("Sample", "Background")
    boxplot(intensities, main = paste("Sample Intersept Comparison", p), col="lightblue")
    legend("topright", c("P-values:", paste("Brown-Forsythe:", round(pv["bf"], decs)), paste("Levene Test" ,round(pv["levene"], decs)), paste("Bartlett's Test" ,round(pv["bartlett"], decs))))
    
    intensities <-ldoscreen@experimentIntercept[[p]]
    pv<-ldoscreen@interceptPvalues[p,]
    names(intensities)<- c("Rat", "Ventilator")
    boxplot(intensities, main = paste("Rat Residual Comparison", p), col="lightblue")
    legend("topright", c("P-values:", paste("Brown-Forsythe:", round(pv["bf"], decs)), paste("Levene Test" ,round(pv["levene"], decs)), paste("Bartlett's Test" ,round(pv["bartlett"], decs))))
    
  }

  

  
}


getPvalue=function(y, group, test){
  if(test=="bf"){
    pv<- levene.test(y, group)$p.value
  }else if(test=="levene"){
    pv<- levene.test(y, group, location="mean")$p.value
  }else if(test=="bartlett"){
    pv<- bartlett.test(x=y, g=group)$p.value
  }else{
    stop("Error: test needs to be either of bf, levene or bartlett"); 
  }
  
  return(pv)
}
