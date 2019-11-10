rm(list=ls())
require(sp)
require(raster)
require(maptools)
require(rgdal)
require(dismo)
require(rJava)
require(arm)
require(plyr)
require(dplyr, warn.conflicts = FALSE)
require(sf)
require(ENMeval)
require(randomForest)
require(caret)
require(pdp)
require(ggplot2)
require(cowplot)
require(splines)
require(rfUtilities)

#To deal with functions with the same name in tidyr
.rs.unloadPackage("tidyr")
#To deal with java issues
.jinit()

wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
#wd <- "~/Desktop/Coastlight/SDM"
setwd(wd)

##This part is to be run remotely
#Set species types for observation data.
speciesList <- c("Grunion","Plover")
for(species in speciesList){
  # Read in grunion or plover presence and pseudo-absence points.
  if(species=="Grunion"){
    obs.data <- read.csv(file="RandomGrunionPointsWGS84.csv")
    abs.data <- read.csv(file="RandomGrunionAbsencesWGS84.csv")
    # Read in environmental map layers.
    env.files <- list.files(pattern="10mAligned.tif$",full.names=TRUE)
    env.files <- env.files[env.files != "./DEM10mAligned.tif" & env.files != "./Saltwater10mAligned.tif"]
  }
  if(species=="Plover"){
    obs.data <- read.csv(file="RandomPloverPointsWGS84.csv")
    abs.data <- read.csv(file="RandomPloverAbsencesWGS84.csv")
    # Read in environmental map layers.
    env.files <- list.files(pattern="10mAligned.tif$",full.names=TRUE)
  }
  
  # Drop unused columns
  obs.data <- obs.data[, c("xcoord", "ycoord")]
  abs.data <- abs.data[, c("xcoord", "ycoord")]
  
  # Read in environmental map layers.
  env.data <- stack(c(env.files))
  # Initialize data containing environmental layer values at presence and pseudo-absence locations.
  presvals <- obs.data
  absvals <- abs.data
  for(env.file in env.files){
    # Get layer names for data frame.
    env.filename <- gsub("^./","",gsub(".tif","",env.file))
    # Extract environmental map layer values at presence points.
    tmp1 <- as.data.frame(extract(raster(env.file),obs.data))
    colnames(tmp1) <- env.filename
    presvals <- cbind(presvals,tmp1)
    # Extract environmental map layer values at presence points.
    tmp2 <- as.data.frame(extract(raster(env.file),abs.data))
    colnames(tmp2) <- env.filename
    absvals <- cbind(absvals,tmp2)
  }
  
  # Standardize missing data
  presvals[is.na(presvals)] <- NA
  presvals <- na.omit(presvals)
  absvals[is.na(absvals)] <- NA
  absvals <- na.omit(absvals)
  
  #Intialize summary statistics variables
  RFImportanceTotal <- data.frame()
  RFEvaluationTotal <- data.frame()
  RFp1Total <- data.frame()
  RFp2Total <- data.frame()
  RFp3Total <- data.frame()
  RFp4Total <- data.frame()
  RFp5Total <- data.frame()
  RFp6Total <- data.frame()
  RFp7Total <- data.frame()
  RFp8Total <- data.frame()
  for(i in 1:3){
    #Subsample presence and pseudo-absence points for training and testing sets for SDMs.
    sample_Num <- 100
    presSubset <- presvals
    presSubset$pa <- 1
    presSubset <- presSubset[,c(ncol(presSubset),1:ncol(presSubset)-1)]
    presSubset <- presSubset[sample(nrow(presSubset),sample_Num),]
    absSubset <- absvals
    absSubset$pa <- 0
    absSubset <- absSubset[sample(nrow(absSubset),10*sample_Num),]
    absSubset <- absSubset[,c(ncol(absSubset),1:ncol(absSubset)-1)]
    
    #Construct a training and testing set for the presence data.
    group <- kfold(presSubset,5)
    pres_train <- presSubset[group!=1,]
    pres_test <- presSubset[group==1,]
    
    #Construct a training and testing set for the pseudo-absence data.
    group <- kfold(absSubset,5)
    backgr_train <- absSubset[group!=1,]
    backgr_test <- absSubset[group==1,]
    
    #Construct presence / pseudo-absence training sets.
    envtrain <- rbind(pres_train,backgr_train)
    envtrain$SoCalBeachType10mAligned <- factor(envtrain$SoCalBeachType10mAligned,levels=1:6)
    testpres <- pres_test
    testbackgr <- backgr_test
    testpres$SoCalBeachType10mAligned <- factor(testpres$SoCalBeachType10mAligned,levels=1:6)
    testbackgr$SoCalBeachType10mAligned <- factor(testbackgr$SoCalBeachType10mAligned,levels=1:6)
    
    #Parsimonious random forest model
    rf.regress <- rf.modelSel(envtrain[,4:ncol(envtrain)],envtrain[,c(1)], imp.scale="mir", parsimony=0.03,final.model=TRUE,seed=1)
    RFImportance <- importance(rf.regress$rf.final)
    RFImportance <- data.frame(names=row.names(RFImportance),RFImportance)
    RFImportanceTotal <- rbind(RFImportanceTotal,RFImportance)
    erf <- suppressWarnings(evaluate(testpres,testbackgr,rf.regress$rf.final))
    RFEvaluation <-  data.frame(matrix(nrow=1,ncol=5))
    colnames(RFEvaluation) <- c("AUC","cor","kappa","Q","TSS")
    RFEvaluation$AUC <- erf@auc
    RFEvaluation$cor <- erf@cor
    RFEvaluation$kappa <- max(erf@kappa)
    # Calculate Yule's Q.
    tmp <- erf@OR
    tmp[!is.finite(tmp)] <- NA 
    RFEvaluation$Q <- (mean(tmp,na.rm=T)-1)/(mean(tmp,na.rm=T)+1)
    RFEvaluation$TSS <- mean(erf@TPR,na.rm=T)+mean(erf@TNR,na.rm=T)-1
    RFEvaluationTotal <- rbind(RFEvaluationTotal,RFEvaluation)
    #Store partial response data for each environmental factor in the random forest model.
    if(species=="Plover"){
      RFp1 <- partial(rf.regress$rf.final,pred.var = "DEM10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp1Total <- rbind(RFp1Total,RFp1)
      RFp3 <- partial(rf.regress$rf.final,pred.var = "Saltwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp3Total <- rbind(RFp3Total,RFp3)
    }
    RFp2 <- partial(rf.regress$rf.final,pred.var = "Freshwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp2Total <- rbind(RFp2Total,RFp2)
    RFp4 <- partial(rf.regress$rf.final,pred.var = "LogSI10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp4Total <- rbind(RFp4Total,RFp4)
    RFp5 <- partial(rf.regress$rf.final,pred.var = "Slope10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp5Total <- rbind(RFp5Total,RFp5)
    RFp6 <- partial(rf.regress$rf.final,pred.var = "SoCalBeachType10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp6Total <- rbind(RFp6Total,RFp6)
    RFp7 <- partial(rf.regress$rf.final,pred.var = "SoCalBeachWidth10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp7Total <- rbind(RFp7Total,RFp7)
    RFp8 <- partial(rf.regress$rf.final,pred.var = "SVF10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp8Total <- rbind(RFp8Total,RFp8)
  }
  #Summary statistics on the frequency and importance of environmental parameters in random forest model.
  tmp <- as.data.frame(table(RFImportanceTotal$names))
  colnames(tmp) <- c("Variable","Freq")
  RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
  colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
  RFImportanceTotal <- left_join(tmp,RFImportanceTotal)
  #To save aggregated data frame.
  write.table(RFImportanceTotal,paste(species,"ParsimonyRFImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(RFEvaluationTotal)
  tmpSD <- apply(RFEvaluationTotal,2,sd)
  RFEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(RFEvaluationTotal,paste(species,"ParsimonyRFEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
  #To save all of the raw points used in the partial response plots for the random forest model.
  if(species=="Plover"){
    colnames(RFp1Total)[which(names(RFp1Total) == "yhat")] <- paste(colnames(RFp1Total)[2],colnames(RFp1Total)[1],sep="")
    colnames(RFp3Total)[which(names(RFp3Total) == "yhat")] <- paste(colnames(RFp3Total)[2],colnames(RFp3Total)[1],sep="")
  }
  colnames(RFp2Total)[which(names(RFp2Total) == "yhat")] <- paste(colnames(RFp2Total)[2],colnames(RFp2Total)[1],sep="")
  colnames(RFp4Total)[which(names(RFp4Total) == "yhat")] <- paste(colnames(RFp4Total)[2],colnames(RFp4Total)[1],sep="")
  colnames(RFp5Total)[which(names(RFp5Total) == "yhat")] <- paste(colnames(RFp5Total)[2],colnames(RFp5Total)[1],sep="")
  colnames(RFp6Total)[which(names(RFp6Total) == "yhat")] <- paste(colnames(RFp6Total)[2],colnames(RFp6Total)[1],sep="")
  colnames(RFp7Total)[which(names(RFp7Total) == "yhat")] <- paste(colnames(RFp7Total)[2],colnames(RFp7Total)[1],sep="")
  colnames(RFp8Total)[which(names(RFp8Total) == "yhat")] <- paste(colnames(RFp8Total)[2],colnames(RFp8Total)[1],sep="")
  if(species=="Plover"){
    RFTotal <- bind_cols(RFp1Total,RFp2Total,RFp3Total,RFp4Total,RFp5Total,RFp7Total,RFp8Total)
  }
  if(species=="Grunion"){
    RFTotal <- bind_cols(RFp2Total,RFp4Total,RFp5Total,RFp7Total,RFp8Total)
  }
  write.table(RFTotal,paste(species,"RFParsimonyPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  write.table(RFp6Total,paste(species,"RFParsimonyBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  #Create 2d histograms with best-fit splines for the partial response curves.
  RFTotal <- read.table(paste(species,"RFParsimonyPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  RFp6Total <- read.table(paste(species,"RFParsimonyBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  if(species=="Plover"){
    RFp1Plot <- ggplot(RFTotal, aes(x=DEM10mAligned, y=yhatDEM10mAligned) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatDEM10mAligned, fill=yhatDEM10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp3Plot <- ggplot(RFTotal, aes(x=Saltwater10mAligned, y=yhatSaltwater10mAligned) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSaltwater10mAligned, fill=yhatSaltwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  }
  RFp2Plot <- ggplot(RFTotal, aes(x=Freshwater10mAligned, y=yhatFreshwater10mAligned) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatFreshwater10mAligned, fill=yhatFreshwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  RFp4Plot <- ggplot(RFTotal, aes(x=LogSI10mAligned, y=yhatLogSI10mAligned) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatLogSI10mAligned, fill=yhatLogSI10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  RFp5Plot <- ggplot(RFTotal, aes(x=Slope10mAligned, y=yhatSlope10mAligned) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSlope10mAligned, fill=yhatSlope10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  RFp6Plot <- ggplot(RFp6Total, aes(x=as.factor(SoCalBeachType10mAligned),y=yhatSoCalBeachType10mAligned))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
  RFp7Plot <- ggplot(RFTotal, aes(x=SoCalBeachWidth10mAligned, y=yhatSoCalBeachWidth10mAligned) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSoCalBeachWidth10mAligned, fill=yhatSoCalBeachWidth10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  RFp8Plot <- ggplot(RFTotal, aes(x=SVF10mAligned, y=yhatSVF10mAligned) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSVF10mAligned, fill=yhatSVF10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  if(species=="Plover"){
    RFPlots <- plot_grid(RFp1Plot,RFp2Plot,RFp3Plot,RFp4Plot,RFp5Plot,RFp6Plot,RFp7Plot,RFp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Grunion"){
    RFPlots <- plot_grid(RFp2Plot,RFp4Plot,RFp5Plot,RFp6Plot,RFp7Plot,RFp8Plot,ncol=2,labels="AUTO")
  }
  png(paste(species,"RFParsimonyPlots.png",sep=""),width=2*800,height=400*length(env.files))
  RFPlots
  print(RFPlots)
  dev.off()
}
