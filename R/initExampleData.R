#

#' @title Create example data set for LoBrA
#' 
#' @description Real signals and background noise originating from experimental settings or random events 
#' @param component vector numbers of background and informative components to be created.
#' @param samples number of experiments per class
#' @param classes number of classes
#' @param bg number of background measurements 
#' @param timepoints number of sample measurements for each spline
#' @param myfile filename of the pdf file created.
#' @export Use this function 
#' @return finalMatrix matrix of example data...
#' @examples \dontrun{
#' } 
#'   components = c(70,10,10,10)
#'   samples = 10
#'   classes = 2;
#'   bg = 3; 
#'   timepoints = rep(5,3)
#'   myfile = "D:/Dropbox/work/projects/LoBrA/results/sample"
#'   p=TRUE;
#'   longDataExample <- createExampleData(components, samples, classes, bg, timepoints)
#'   dim(longDataExample)
#'   save(longDataExample,file="data/longDataExample.RData")
#'   
#' 
createExampleData<-function(components=c(70,10,10,10), samples=10, classes=2, bg=3, timepoints=rep(5,3), myfile=NA) {
  p<-!is.na(myfile)
  
  finalMatrix<-c()
  
  labels<-rep(1:classes,each=samples)
  s=2*samples
  # labels<-rep(1:classes,each=5)
  
  if(p){
    pdf(file = paste(myfile,"bg.pdf",  sep = ""), width=8, height=6)
    bgComp<-createBGComponents(components=components[1], samples=s, labels=labels, timepoints=15, bg=3, mean=5, sd=3, experimentSD=4, randomnoise=0.1, plotting=TRUE)
    dev.off();
  }else{
    bgComp<-createBGComponents(components=components[1], samples=s, labels=labels, timepoints=15, bg=3, mean=5, sd=3, experimentSD=4, randomnoise=0.1, plotting=FALSE)
  }
  finalMatrix<-bgComp
  dim(bgComp)
  head(bgComp)
  se=1
  if(p){
    pdf(file = paste(myfile,"IC.pdf",  sep = ""), width=8, height=6)
  }
  for(se in 1:3){
    iC<- createInformativeComponents(components=components[1+se], samples=s, labels=labels, timepoints=c(5,5,5), bg=3, mean=5, sd=3, segment=se, myl=labels[1], slopeSD=2, randomnoise=0.5, plotting=p)
    finalMatrix<-cbind(finalMatrix, iC[,-c(1:4)])
  }
  if(p){
    dev.off();
  }
  
  finalMatrix[,"id"]<-paste("Exp",finalMatrix[,"id"], sep="")
  
  return(finalMatrix)
}


#' Create example informative peaks
#' 
#' Background noise signals originating from experimental settings or random events 
#' @param component number of background components to be created 
#' @param samples number of experiments 
#' @param label label of each experiment
#' @param timepoints number of sample measurements 
#' @param bg number of background measurements 
#' @param mean mean value of noise components
#' @param sd standard deviation value of noise for this component
#' @param experiementSD standard deviation value of each experiment for this component
#' @param randomnoise random variation changing at each time point
#' @return matrix of background components
#' 
createInformativeComponents<-function(components=10, samples=10, labels=NA, timepoints=c(5,5,5), bg=3, mean=5, sd=3, segment=1, myl=labels[1], slopeSD=2, randomnoise=0.5, plotting=FALSE){
  alltimes <- bg + sum(timepoints);
  types<- c(rep("b",bg), rep("s",sum(timepoints)))
  s<- length(labels)
  finalMatrix <- data.frame(rep(1:s,each=alltimes), rep(1:alltimes,s), rep(types, s), rep(labels, each=alltimes))
  
  initialvalues<-rnorm(components,mean,sd);
  shiftvalues<-rnorm(components,0,1);
  m<-2
  for(m in 1:components){
    backgroundMatrix<- createBGData(s, bg, mean=initialvalues[m], sd=1, randomnoise)
    background<-backgroundMatrix[,1]
    background<-background +  rnorm(s, 0, sd=randomnoise)
    
    
    ts=1
    for(ts in 1:length(timepoints)){
      if(ts==segment){
        componentSlope<-rnorm(1, 0, sd=slopeSD);
        for(t in 1:timepoints[ts]){
          background<- background + rnorm(s, 0, sd=randomnoise)
          ids<-labels==myl
          background[ids]<- background[ids] + componentSlope
          backgroundMatrix<-cbind(backgroundMatrix, background)
        }
      }else{
        for(t in 1:timepoints[ts]){
          background<- background + rnorm(s, 0, sd=randomnoise)
          backgroundMatrix<-cbind(backgroundMatrix, background)
        }
      }
    }
    colnames(backgroundMatrix)<- 1:dim(backgroundMatrix)[2]
    if(plotting){
      plotTimeSeries(backgroundMatrix, main=paste("Component",m,", Segment",segment), labels=labels);
    }
    
    oneColumn<-c()
    for(s in 1:s){
      oneColumn<-c(oneColumn, backgroundMatrix[s,])
    }
    
    
    finalMatrix<- cbind(finalMatrix, oneColumn)
  }
  colnames(finalMatrix)<- c("id", "time", "type", "class", paste("component-", 1:components, "-", segment, sep=""))
  
  return(finalMatrix);
}

#' Create example background noise peaks
#' 
#' Background noise signals originating from experimental settings or random events 
#' @param component number of background components to be created 
#' @param samples number of experiments 
#' @param label label of each experiment
#' @param timepoints number of sample measurements 
#' @param bg number of background measurements 
#' @param mean mean value of noise components
#' @param sd standard deviation value of noise for this component
#' @param experiementSD standard deviation value of each experiment for this component
#' @param randomnoise random variation changing at each time point
#' @return matrix of background components
#' 
createBGComponents<-function(components=50, samples=10, labels=NA, timepoints=15, bg=3, mean=5, sd=3, experimentSD=2, randomnoise=0.1, plotting=FALSE){
  alltimes <- bg + timepoints;
  types<- c(rep("b",bg), rep("s",timepoints))
  finalMatrix <- data.frame(rep(1:samples,each=alltimes), rep(1:alltimes,samples), rep(types, samples), rep(labels, each=alltimes))
  
  initialvalues<-rnorm(components,mean,sd);
  # shiftvalues<-rnorm(components,0,1);
  m<-1
  for(m in 1:components){
    backgroundMatrix<- createBGData(samples, bg, mean=initialvalues[m], sd=experimentSD, randomnoise*2)
    background<-backgroundMatrix[,1]
    background<-background +  rnorm(samples, 0, sd=randomnoise)
    
    for(t in 1:timepoints){
      background<- background + rnorm(samples, 0, sd=randomnoise)
      backgroundMatrix<-cbind(backgroundMatrix, background)
    }
    colnames(backgroundMatrix)<- 1:dim(backgroundMatrix)[2]
    if(plotting){
      plotTimeSeries(backgroundMatrix, main=paste("BG Component",m), labels=labels);
    }
    
    oneColumn<-c()
    for(s in 1:samples){
      oneColumn<-c(oneColumn, backgroundMatrix[s,])
    }
    
    
    finalMatrix<- cbind(finalMatrix, oneColumn)
  }
  colnames(finalMatrix)<- c("id", "time", "type", "class", paste("bgcomponent-",1:components, sep=""))
  
  return(finalMatrix);
}


#' Create example background measurements
#' 
#' Background noise signals originating from experimental settings or random events 
#' @param samples number of experiments 
#' @param bg number of background measurements 
#' @param mean mean value of noise for this component
#' @param sd standard deviation value of noise for this component
#' @param randomnoise random variation changing at each time point
#' @return matrix of background measurements
#' 
createBGData<-function(samples=10, bg=3, mean=0, sd=1, randomnoise=0.1){
  background<-rnorm(samples,mean,sd);
  backgroundMatrix<-background
  for(i in 2:bg){
    background<-background+rnorm(samples, 0, randomnoise)
    backgroundMatrix<-cbind(backgroundMatrix, background)
  }
  colnames(backgroundMatrix)<-1:bg
  rownames(backgroundMatrix)<-1:samples
  
  #plotTimeSeries(backgroundMatrix)
  return(backgroundMatrix);
}


#' Plotting Function for a longitudinal data matrix
#' 
#' @param myMatrix matrix to be plotted
#' @param labels class labels of samples
#' @param ylab Label of y achis
#' @param xlab Label of x achis
#' @param legend of plot
#' @param col vector of colors for plot
#' 
plotTimeSeries=function(myMatrix, main="", labels=NA, ylab="Expression", xlab="Time Point", legend="", col=1:dim(myMatrix)[1]){
  ylim<-c(min(myMatrix),max(myMatrix));
  times<-as.numeric(colnames(myMatrix))
  xlim<-c(min(times),max(times));
  mycol<-1:dim(myMatrix)[1]
  if(sum(!is.na(labels))>0){
    ul<- unique(labels)
    mycol<-rep("", length(labels))
    for(l in ul){
      mycol[labels==l]<- getColor(l, sum(labels==l))
    }
  }
  plot(x=xlim, y=c(ylim[1], ylim[1]), type="l", col="white", xlim=xlim, ylim=ylim, lwd=2, main=main, ylab=ylab, xlab=xlab)
  for(l in 1:dim(myMatrix)[1]){
    # print(rat$name)
    lines(x=times, y=myMatrix[l,], col=mycol[l], lty=l, lwd=2)
  }
  
  legend("topright", legend=1:dim(myMatrix)[1], col=mycol, lty=1:dim(myMatrix)[1], lwd=2, cex=0.5, bg="white")

  }

#' Get colors for the plotting function.
#' 
#' @param label class labels of the samples
#' @param size size of the color vector to be created
#' @return col vector of colors created
#' 
getColor=function(label, size){
  size<-10
  if(label==1){
    col<- colorRampPalette(c("khaki", "yellow", "darkgoldenrod3"))(size);
  }else   if(label==2){
    col<- colorRampPalette(c("cadetblue1", "darkblue"))(size);
  }else   if(label==3){
    col<- colorRampPalette(c("lightgreen", "green","darkgreen"))(size);
  }else {
    col<- colorRampPalette(c("pink", "red", "darkred"))(size);
  }
  return(col)
}
