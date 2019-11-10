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

#To deal with functions with the same name in tidyr
.rs.unloadPackage("tidyr")
#To deal with java issues
.jinit()

#wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
wd <- "~/Desktop/Coastlight/SDM"
setwd(wd)

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
  GLMImportanceTotal <- data.frame()
  GLMEvaluationTotal <- data.frame()
  GLMp1Total <- data.frame()
  GLMp2Total <- data.frame()
  GLMp3Total <- data.frame()
  GLMp4Total <- data.frame()
  GLMp5Total <- data.frame()
  GLMp6Total <- data.frame()
  GLMp7Total <- data.frame()
  GLMp8Total <- data.frame()
  XMImportanceTotal <- data.frame()
  XMEvaluationTotal <- data.frame()
  XMp1Total <- data.frame()
  XMp2Total <- data.frame()
  XMp3Total <- data.frame()
  XMp4Total <- data.frame()
  XMp5Total <- data.frame()
  XMp6Total <- data.frame()
  XMp7Total <- data.frame()
  XMp8Total <- data.frame()
  set.seed(1)
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
    
    #Random forest model
    rf1 <- suppressWarnings(tuneRF(envtrain[,4:ncol(envtrain)],envtrain[,c(1)],stepFactor=1,plot=FALSE,doBest=TRUE))
    RFImportance <- importance(rf1)
    RFImportance <- data.frame(names=row.names(RFImportance),RFImportance)
    RFImportanceTotal <- rbind(RFImportanceTotal,RFImportance)
    erf <- suppressWarnings(evaluate(testpres,testbackgr,rf1))
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
    #Store partial response curve for random forest model.
    if(species=="Grunion"){
      #RFp1 <- partial(rf1,pred.var = "DEM10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      #RFp1Total <- rbind(RFp1Total,RFp1)
      RFp2 <- partial(rf1,pred.var = "Freshwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp2Total <- rbind(RFp2Total,RFp2)
      #RFp3 <- partial(rf1,pred.var = "Saltwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      #RFp3Total <- rbind(RFp3Total,RFp3)
      RFp4 <- partial(rf1,pred.var = "LogSI10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp4Total <- rbind(RFp4Total,RFp4)
      RFp5 <- partial(rf1,pred.var = "Slope10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp5Total <- rbind(RFp5Total,RFp5)
      RFp6 <- partial(rf1,pred.var = "SoCalBeachType10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp6Total <- rbind(RFp6Total,RFp6)
      RFp7 <- partial(rf1,pred.var = "SoCalBeachWidth10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp7Total <- rbind(RFp7Total,RFp7)
      RFp8 <- partial(rf1,pred.var = "SVF10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp8Total <- rbind(RFp8Total,RFp8)
    }
    if(species=="Plover"){
      RFp1 <- partial(rf1,pred.var = "DEM10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp1Total <- rbind(RFp1Total,RFp1)
      RFp2 <- partial(rf1,pred.var = "Freshwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp2Total <- rbind(RFp2Total,RFp2)
      RFp3 <- partial(rf1,pred.var = "Saltwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp3Total <- rbind(RFp3Total,RFp3)
      RFp4 <- partial(rf1,pred.var = "LogSI10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp4Total <- rbind(RFp4Total,RFp4)
      RFp5 <- partial(rf1,pred.var = "Slope10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp5Total <- rbind(RFp5Total,RFp5)
      RFp6 <- partial(rf1,pred.var = "SoCalBeachType10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp6Total <- rbind(RFp6Total,RFp6)
      RFp7 <- partial(rf1,pred.var = "SoCalBeachWidth10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp7Total <- rbind(RFp7Total,RFp7)
      RFp8 <- partial(rf1,pred.var = "SVF10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp8Total <- rbind(RFp8Total,RFp8)
    }
    
    #Generalized linear model
    m1 <- glm(pa ~ .-1, data=envtrain[,-c(2,3)],family = binomial(link = "logit"))
    #Determine the variable relative importance using the absolute value of the z-statistic.
    tmp <- summary(m1)
    GLMImportance <- as.data.frame(tmp$coefficients)
    GLMImportance <- abs(GLMImportance)
    GLMImportance <- data.frame(names=row.names(GLMImportance),GLMImportance)
    GLMImportance <- GLMImportance[,c(1,4)]
    GLMImportanceTotal <- rbind(GLMImportanceTotal,GLMImportance)
    em1 <- suppressWarnings(evaluate(testpres,testbackgr,m1))
    GLMEvaluation <- data.frame(matrix(nrow=1,ncol=5))
    colnames(GLMEvaluation) <- c("AUC","cor","kappa","Q","TSS")
    GLMEvaluation$AUC <- em1@auc
    GLMEvaluation$cor <- em1@cor
    GLMEvaluation$kappa <- max(em1@kappa)
    # Calculate Yule's Q.
    tmp <- em1@OR
    tmp[!is.finite(tmp)] <- NA 
    GLMEvaluation$Q <- (mean(tmp,na.rm=T)-1)/(mean(tmp,na.rm=T)+1)
    GLMEvaluation$TSS <- mean(em1@TPR,na.rm=T)+mean(em1@TNR,na.rm=T)-1
    GLMEvaluationTotal <- rbind(GLMEvaluationTotal,GLMEvaluation)
    #Store partial response curve for generalized linear model.
    #Note that raw outputs are logit units.  To convert to observation probabilities:
    # prob = exp(logit)/(1+exp(logit))
    if(species=="Grunion"){
      #GLMp1 <- partial(m1,pred.var = "DEM10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      #GLMp1Total <- rbind(GLMp1Total,GLMp1)
      GLMp2 <- partial(m1,pred.var = "Freshwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp2Total <- rbind(GLMp2Total,GLMp2)
      #GLMp3 <- partial(m1,pred.var = "Saltwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      #GLMp3Total <- rbind(GLMp3Total,GLMp3)
      GLMp4 <- partial(m1,pred.var = "LogSI10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp4Total <- rbind(GLMp4Total,GLMp4)
      GLMp5 <- partial(m1,pred.var = "Slope10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp5Total <- rbind(GLMp5Total,GLMp5)
      GLMp6 <- partial(m1,pred.var = "SoCalBeachType10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp6Total <- rbind(GLMp6Total,GLMp6)
      GLMp7 <- partial(m1,pred.var = "SoCalBeachWidth10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp7Total <- rbind(GLMp7Total,GLMp7)
      GLMp8 <- partial(m1,pred.var = "SVF10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp8Total <- rbind(GLMp8Total,GLMp8)
    }
    if(species=="Plover"){
      GLMp1 <- partial(m1,pred.var = "DEM10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp1Total <- rbind(GLMp1Total,GLMp1)
      GLMp2 <- partial(m1,pred.var = "Freshwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp2Total <- rbind(GLMp2Total,GLMp2)
      GLMp3 <- partial(m1,pred.var = "Saltwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp3Total <- rbind(GLMp3Total,GLMp3)
      GLMp4 <- partial(m1,pred.var = "LogSI10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp4Total <- rbind(GLMp4Total,GLMp4)
      GLMp5 <- partial(m1,pred.var = "Slope10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp5Total <- rbind(GLMp5Total,GLMp5)
      GLMp6 <- partial(m1,pred.var = "SoCalBeachType10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp6Total <- rbind(GLMp6Total,GLMp6)
      GLMp7 <- partial(m1,pred.var = "SoCalBeachWidth10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp7Total <- rbind(GLMp7Total,GLMp7)
      GLMp8 <- partial(m1,pred.var = "SVF10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      GLMp8Total <- rbind(GLMp8Total,GLMp8)
    }
    
    #MaxEnt model
    maxent()
    xm <- maxent(x=env.data,p=presSubset[,c(2,3)],a=absSubset[,c(2,3)],factors='SoCalBeachType10mAligned')
    XMImportance <- var.importance(xm)
    XMImportance <- XMImportance[,c(1,3)]
    XMImportanceTotal <- rbind(XMImportanceTotal,XMImportance)
    exm <- suppressWarnings(evaluate(testpres,testbackgr,xm))
    XMEvaluation <- data.frame(matrix(nrow=1,ncol=5))
    colnames(XMEvaluation) <- c("AUC","cor","kappa","Q","TSS")
    XMEvaluation$AUC <- exm@auc
    XMEvaluation$cor <- exm@cor
    XMEvaluation$kappa <- max(exm@kappa)
    # Calculate Yule's Q.
    tmp <- exm@OR
    tmp[!is.finite(tmp)] <- NA 
    XMEvaluation$Q <- (mean(tmp,na.rm=T)-1)/(mean(tmp,na.rm=T)+1)
    XMEvaluation$TSS <- mean(exm@TPR,na.rm=T)+mean(exm@TNR,na.rm=T)-1
    XMEvaluationTotal <- rbind(XMEvaluationTotal,XMEvaluation)
    #Store the variable response functions for Maxent model.
    if(species=="Grunion"){
      #XMp1 <- response(xm,var="DEM10mAligned",expand=0,range="pa")
      #XMp1Total <- rbind(XMp1Total,XMp1)
      XMp2 <- response(xm,var="Freshwater10mAligned",expand=0,range="pa")
      XMp2Total <- rbind(XMp2Total,XMp2)
      #XMp3 <- response(xm,var="Saltwater10mAligned",expand=0,range="pa")
      #XMp3Total <- rbind(XMp3Total,XMp3)
      XMp4 <- response(xm,var="LogSI10mAligned",expand=0,range="pa")
      XMp4Total <- rbind(XMp4Total,XMp4)
      XMp5 <- response(xm,var="Slope10mAligned",expand=0,range="pa")
      XMp5Total <- rbind(XMp5Total,XMp5)
      XMp6 <- response(xm,var="SoCalBeachType10mAligned",expand=0,range="pa")
      XMp6 <- XMp6[!duplicated(XMp6),]
      XMp6Total <- rbind(XMp6Total,XMp6)
      XMp7 <- response(xm,var="SoCalBeachWidth10mAligned",expand=0,range="pa")
      XMp7Total <- rbind(XMp7Total,XMp7)
      XMp8 <- response(xm,var="SVF10mAligned",expand=0,range="pa")
      XMp8Total <- rbind(XMp8Total,XMp8)
    }
    if(species=="Plover"){
      XMp1 <- response(xm,var="DEM10mAligned",expand=0,range="pa")
      XMp1Total <- rbind(XMp1Total,XMp1)
      XMp2 <- response(xm,var="Freshwater10mAligned",expand=0,range="pa")
      XMp2Total <- rbind(XMp2Total,XMp2)
      XMp3 <- response(xm,var="Saltwater10mAligned",expand=0,range="pa")
      XMp3Total <- rbind(XMp3Total,XMp3)
      XMp4 <- response(xm,var="LogSI10mAligned",expand=0,range="pa")
      XMp4Total <- rbind(XMp4Total,XMp4)
      XMp5 <- response(xm,var="Slope10mAligned",expand=0,range="pa")
      XMp5Total <- rbind(XMp5Total,XMp5)
      XMp6 <- response(xm,var="SoCalBeachType10mAligned",expand=0,range="pa")
      XMp6 <- XMp6[!duplicated(XMp6),]
      XMp6Total <- rbind(XMp6Total,XMp6)
      XMp7 <- response(xm,var="SoCalBeachWidth10mAligned",expand=0,range="pa")
      XMp7Total <- rbind(XMp7Total,XMp7)
      XMp8 <- response(xm,var="SVF10mAligned",expand=0,range="pa")
      XMp8Total <- rbind(XMp8Total,XMp8)
    }
  }
  #Random forest summary statistics
  RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
  #RFImportanceTotal <- left_join(tmpMean,tmpSD,by=c("names"))
  colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
  #To save aggregated data frame.
  write.table(RFImportanceTotal,paste(species,"RFImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(RFEvaluationTotal)
  tmpSD <- apply(RFEvaluationTotal,2,sd)
  RFEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(RFEvaluationTotal,paste(species,"RFEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
  
  #GLM summary statistics
  GLMImportanceTotal <- ddply(GLMImportanceTotal, .(names), summarize,  MeanImportance=mean(z.value), SDImportance=sd(z.value))
  #GLMImportanceTotal <- left_join(tmpMean,tmpSD,by=c("names"))
  colnames(GLMImportanceTotal) <- c("Variable","MeanImportance","SDImportance")
  #To save aggregated data frame.
  write.table(GLMImportanceTotal,paste(species,"GLMImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(GLMEvaluationTotal)
  tmpSD <- apply(GLMEvaluationTotal,2,sd)
  GLMEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(GLMEvaluationTotal,paste(species,"GLMEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
  
  #Maxent summary statistics
  XMImportanceTotal <- ddply(XMImportanceTotal, .(variable), summarize,  MeanImportance=mean(permutation.importance), SDImportance=sd(permutation.importance))
  #XMImportanceTotal <- left_join(tmpMean,tmpSD,by=c("variable"))
  colnames(XMImportanceTotal) <- c("Variable","MeanImportance","SDImportance")
  #To save aggregated data frame.
  write.table(XMImportanceTotal,paste(species,"XMImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(XMEvaluationTotal)
  tmpSD <- apply(XMEvaluationTotal,2,sd)
  XMEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(XMEvaluationTotal,paste(species,"XMEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
  #To save all of the raw points used in the partial response plots for the generalized linear model.
  if(species=="Grunion"){
    #colnames(GLMp1Total)[which(names(GLMp1Total) == "yhat")] <- paste(colnames(GLMp1Total)[2],colnames(GLMp1Total)[1],sep="")
    colnames(GLMp2Total)[which(names(GLMp2Total) == "yhat")] <- paste(colnames(GLMp2Total)[2],colnames(GLMp2Total)[1],sep="")
    #colnames(GLMp3Total)[which(names(GLMp3Total) == "yhat")] <- paste(colnames(GLMp3Total)[2],colnames(GLMp3Total)[1],sep="")
    colnames(GLMp4Total)[which(names(GLMp4Total) == "yhat")] <- paste(colnames(GLMp4Total)[2],colnames(GLMp4Total)[1],sep="")
    colnames(GLMp5Total)[which(names(GLMp5Total) == "yhat")] <- paste(colnames(GLMp5Total)[2],colnames(GLMp5Total)[1],sep="")
    colnames(GLMp6Total)[which(names(GLMp6Total) == "yhat")] <- paste(colnames(GLMp6Total)[2],colnames(GLMp6Total)[1],sep="")
    colnames(GLMp7Total)[which(names(GLMp7Total) == "yhat")] <- paste(colnames(GLMp7Total)[2],colnames(GLMp7Total)[1],sep="")
    colnames(GLMp8Total)[which(names(GLMp8Total) == "yhat")] <- paste(colnames(GLMp8Total)[2],colnames(GLMp8Total)[1],sep="")
    GLMTotal <- bind_cols(GLMp2Total,GLMp4Total,GLMp5Total,GLMp7Total,GLMp8Total)
    write.table(GLMTotal,paste(species,"GLMPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(GLMp6Total,paste(species,"GLMBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    #Create 2d histograms with best-fit splines for the partial response curves.
    GLMTotal <- read.table(paste(species,"GLMPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    GLMp6Total <- read.table(paste(species,"GLMBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    #GLMp1Plot <- ggplot(GLMTotal, aes(x=DEM10mAligned, y=exp(yhatDEM10mAligned)/(1+exp(yhatDEM10mAligned))) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatDEM10mAligned)/(1+exp(yhatDEM10mAligned)), fill=exp(yhatDEM10mAligned)/(1+exp(yhatDEM10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp2Plot <- ggplot(GLMTotal, aes(x=Freshwater10mAligned, y=exp(yhatFreshwater10mAligned)/(1+exp(yhatFreshwater10mAligned))) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatFreshwater10mAligned)/(1+exp(yhatFreshwater10mAligned)), fill=exp(yhatFreshwater10mAligned)/(1+exp(yhatFreshwater10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    #GLMp3Plot <- ggplot(GLMTotal, aes(x=Saltwater10mAligned, y=exp(yhatSaltwater10mAligned)/(1+exp(yhatSaltwater10mAligned))) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSaltwater10mAligned)/(1+exp(yhatSaltwater10mAligned)), fill=exp(yhatSaltwater10mAligned)/(1+exp(yhatSaltwater10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp4Plot <- ggplot(GLMTotal, aes(x=LogSI10mAligned, y=exp(yhatLogSI10mAligned)/(1+exp(yhatLogSI10mAligned))) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatLogSI10mAligned)/(1+exp(yhatLogSI10mAligned)), fill=exp(yhatLogSI10mAligned)/(1+exp(yhatLogSI10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp5Plot <- ggplot(GLMTotal, aes(x=Slope10mAligned, y=exp(yhatSlope10mAligned)/(1+exp(yhatSlope10mAligned))) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSlope10mAligned)/(1+exp(yhatSlope10mAligned)), fill=exp(yhatSlope10mAligned)/(1+exp(yhatSlope10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp6Plot <- ggplot(GLMp6Total, aes(x=as.factor(SoCalBeachType10mAligned),y=exp(yhatSoCalBeachType10mAligned)/(1+exp(yhatSoCalBeachType10mAligned))))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
    GLMp7Plot <- ggplot(GLMTotal, aes(x=SoCalBeachWidth10mAligned, y=exp(yhatSoCalBeachWidth10mAligned)/(1+exp(yhatSoCalBeachWidth10mAligned))) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSoCalBeachWidth10mAligned)/(1+exp(yhatSoCalBeachWidth10mAligned)), fill=exp(yhatSoCalBeachWidth10mAligned)/(1+exp(yhatSoCalBeachWidth10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp8Plot <- ggplot(GLMTotal, aes(x=SVF10mAligned, y=exp(yhatSVF10mAligned)/(1+exp(yhatSVF10mAligned))) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSVF10mAligned)/(1+exp(yhatSVF10mAligned)), fill=exp(yhatSVF10mAligned)/(1+exp(yhatSVF10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMPlots <- plot_grid(GLMp2Plot,GLMp4Plot,GLMp5Plot,GLMp6Plot,GLMp7Plot,GLMp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Plover"){
    colnames(GLMp1Total)[which(names(GLMp1Total) == "yhat")] <- paste(colnames(GLMp1Total)[2],colnames(GLMp1Total)[1],sep="")
    colnames(GLMp2Total)[which(names(GLMp2Total) == "yhat")] <- paste(colnames(GLMp2Total)[2],colnames(GLMp2Total)[1],sep="")
    colnames(GLMp3Total)[which(names(GLMp3Total) == "yhat")] <- paste(colnames(GLMp3Total)[2],colnames(GLMp3Total)[1],sep="")
    colnames(GLMp4Total)[which(names(GLMp4Total) == "yhat")] <- paste(colnames(GLMp4Total)[2],colnames(GLMp4Total)[1],sep="")
    colnames(GLMp5Total)[which(names(GLMp5Total) == "yhat")] <- paste(colnames(GLMp5Total)[2],colnames(GLMp5Total)[1],sep="")
    colnames(GLMp6Total)[which(names(GLMp6Total) == "yhat")] <- paste(colnames(GLMp6Total)[2],colnames(GLMp6Total)[1],sep="")
    colnames(GLMp7Total)[which(names(GLMp7Total) == "yhat")] <- paste(colnames(GLMp7Total)[2],colnames(GLMp7Total)[1],sep="")
    colnames(GLMp8Total)[which(names(GLMp8Total) == "yhat")] <- paste(colnames(GLMp8Total)[2],colnames(GLMp8Total)[1],sep="")
    GLMTotal <- bind_cols(GLMp1Total,GLMp2Total,GLMp3Total,GLMp4Total,GLMp5Total,GLMp7Total,GLMp8Total)
    write.table(GLMTotal,paste(species,"GLMPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(GLMp6Total,paste(species,"GLMBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    #Create 2d histograms with best-fit splines for the partial response curves.
    GLMTotal <- read.table(paste(species,"GLMPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    GLMp6Total <- read.table(paste(species,"GLMBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    GLMp1Plot <- ggplot(GLMTotal, aes(x=DEM10mAligned, y=exp(yhatDEM10mAligned)/(1+exp(yhatDEM10mAligned))) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatDEM10mAligned)/(1+exp(yhatDEM10mAligned)), fill=exp(yhatDEM10mAligned)/(1+exp(yhatDEM10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp2Plot <- ggplot(GLMTotal, aes(x=Freshwater10mAligned, y=exp(yhatFreshwater10mAligned)/(1+exp(yhatFreshwater10mAligned))) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatFreshwater10mAligned)/(1+exp(yhatFreshwater10mAligned)), fill=exp(yhatFreshwater10mAligned)/(1+exp(yhatFreshwater10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp3Plot <- ggplot(GLMTotal, aes(x=Saltwater10mAligned, y=exp(yhatSaltwater10mAligned)/(1+exp(yhatSaltwater10mAligned))) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSaltwater10mAligned)/(1+exp(yhatSaltwater10mAligned)), fill=exp(yhatSaltwater10mAligned)/(1+exp(yhatSaltwater10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp4Plot <- ggplot(GLMTotal, aes(x=LogSI10mAligned, y=exp(yhatLogSI10mAligned)/(1+exp(yhatLogSI10mAligned))) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatLogSI10mAligned)/(1+exp(yhatLogSI10mAligned)), fill=exp(yhatLogSI10mAligned)/(1+exp(yhatLogSI10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp5Plot <- ggplot(GLMTotal, aes(x=Slope10mAligned, y=exp(yhatSlope10mAligned)/(1+exp(yhatSlope10mAligned))) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSlope10mAligned)/(1+exp(yhatSlope10mAligned)), fill=exp(yhatSlope10mAligned)/(1+exp(yhatSlope10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp6Plot <- ggplot(GLMp6Total, aes(x=as.factor(SoCalBeachType10mAligned),y=exp(yhatSoCalBeachType10mAligned)/(1+exp(yhatSoCalBeachType10mAligned))))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
    GLMp7Plot <- ggplot(GLMTotal, aes(x=SoCalBeachWidth10mAligned, y=exp(yhatSoCalBeachWidth10mAligned)/(1+exp(yhatSoCalBeachWidth10mAligned))) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSoCalBeachWidth10mAligned)/(1+exp(yhatSoCalBeachWidth10mAligned)), fill=exp(yhatSoCalBeachWidth10mAligned)/(1+exp(yhatSoCalBeachWidth10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMp8Plot <- ggplot(GLMTotal, aes(x=SVF10mAligned, y=exp(yhatSVF10mAligned)/(1+exp(yhatSVF10mAligned))) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = exp(yhatSVF10mAligned)/(1+exp(yhatSVF10mAligned)), fill=exp(yhatSVF10mAligned)/(1+exp(yhatSVF10mAligned))),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(GLMTotal))+theme_bw(base_size=25)
    GLMPlots <- plot_grid(GLMp1Plot,GLMp2Plot,GLMp3Plot,GLMp4Plot,GLMp5Plot,GLMp6Plot,GLMp7Plot,GLMp8Plot,ncol=2,labels="AUTO")
  }
  #To save all of the raw points used in the partial response plots for the generalized linear model.
  if(species=="Grunion"){
    #colnames(XMp1Total) <- c("Elevation","yElevation")
    colnames(XMp2Total) <- c("DistanceToFreshwater","yDistanceToFreshwater")
    #colnames(XMp3Total) <- c("DistanceToSaltwater","yDistanceToSaltwater")
    colnames(XMp4Total) <- c("LogSI","yLogSI")
    colnames(XMp5Total) <- c("Slope","ySlope")
    colnames(XMp6Total) <- c("BeachType","yBeachType")
    colnames(XMp7Total) <- c("BeachWidth","yBeachWidth")
    colnames(XMp8Total) <- c("SVF","ySVF")
    XMTotal <- bind_cols(XMp2Total,XMp4Total,XMp5Total,XMp7Total,XMp8Total)
    print(colnames(XMTotal))
    write.table(XMTotal,paste(species,"XMPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(XMp6Total,paste(species,"XMBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    #Create 2d histograms with best-fit splines for the partial response curves.
    XMTotal <- read.table(paste(species,"XMPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    XMp6Total <- read.table(paste(species,"XMBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    #XMp1Plot <- ggplot(XMTotal, aes(x=Elevation, y=yElevation) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yElevation, fill=yElevation),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp2Plot <- ggplot(XMTotal, aes(x=DistanceToFreshwater, y=yDistanceToFreshwater) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yDistanceToFreshwater, fill=yDistanceToFreshwater),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    #XMp3Plot <- ggplot(XMTotal, aes(x=DistanceToSaltwater, y=yDistanceToSaltwater) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yDistanceToSaltwater, fill=yDistanceToSaltwater),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp4Plot <- ggplot(XMTotal, aes(x=LogSI, y=yLogSI) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yLogSI, fill=yLogSI),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp5Plot <- ggplot(XMTotal, aes(x=Slope, y=ySlope) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = ySlope, fill=ySlope),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp6Plot <- ggplot(XMp6Total, aes(x=as.factor(BeachType),y=yBeachType))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
    XMp7Plot <- ggplot(XMTotal, aes(x=BeachWidth, y=yBeachWidth) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yBeachWidth, fill=yBeachWidth),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp8Plot <- ggplot(XMTotal, aes(x=SVF, y=ySVF) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = ySVF, fill=ySVF),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMPlots <- plot_grid(XMp2Plot,XMp4Plot,XMp5Plot,XMp6Plot,XMp7Plot,XMp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Plover"){
    colnames(XMp1Total) <- c("Elevation","yElevation")
    colnames(XMp2Total) <- c("DistanceToFreshwater","yDistanceToFreshwater")
    colnames(XMp3Total) <- c("DistanceToSaltwater","yDistanceToSaltwater")
    colnames(XMp4Total) <- c("LogSI","yLogSI")
    colnames(XMp5Total) <- c("Slope","ySlope")
    colnames(XMp6Total) <- c("BeachType","yBeachType")
    colnames(XMp7Total) <- c("BeachWidth","yBeachWidth")
    colnames(XMp8Total) <- c("SVF","ySVF")
    XMTotal <- bind_cols(XMp1Total,XMp2Total,XMp3Total,XMp4Total,XMp5Total,XMp7Total,XMp8Total)
    print(colnames(XMTotal))
    write.table(XMTotal,paste(species,"XMPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(XMp6Total,paste(species,"XMBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    #Create 2d histograms with best-fit splines for the partial response curves.
    XMTotal <- read.table(paste(species,"XMPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    XMp6Total <- read.table(paste(species,"XMBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    XMp1Plot <- ggplot(XMTotal, aes(x=Elevation, y=yElevation) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yElevation, fill=yElevation),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp2Plot <- ggplot(XMTotal, aes(x=DistanceToFreshwater, y=yDistanceToFreshwater) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yDistanceToFreshwater, fill=yDistanceToFreshwater),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp3Plot <- ggplot(XMTotal, aes(x=DistanceToSaltwater, y=yDistanceToSaltwater) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yDistanceToSaltwater, fill=yDistanceToSaltwater),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp4Plot <- ggplot(XMTotal, aes(x=LogSI, y=yLogSI) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yLogSI, fill=yLogSI),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp5Plot <- ggplot(XMTotal, aes(x=Slope, y=ySlope) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = ySlope, fill=ySlope),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp6Plot <- ggplot(XMp6Total, aes(x=as.factor(BeachType),y=yBeachType))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
    XMp7Plot <- ggplot(XMTotal, aes(x=BeachWidth, y=yBeachWidth) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yBeachWidth, fill=yBeachWidth),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMp8Plot <- ggplot(XMTotal, aes(x=SVF, y=ySVF) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = ySVF, fill=ySVF),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(XMTotal))+theme_bw(base_size=25)
    XMPlots <- plot_grid(XMp1Plot,XMp2Plot,XMp3Plot,XMp4Plot,XMp5Plot,XMp6Plot,XMp7Plot,XMp8Plot,ncol=2,labels="AUTO")
  }
  #To save all of the raw points used in the partial response plots for the random forest model.
  if(species=="Grunion"){
    #colnames(RFp1Total)[which(names(RFp1Total) == "yhat")] <- paste(colnames(RFp1Total)[2],colnames(RFp1Total)[1],sep="")
    colnames(RFp2Total)[which(names(RFp2Total) == "yhat")] <- paste(colnames(RFp2Total)[2],colnames(RFp2Total)[1],sep="")
    #colnames(RFp3Total)[which(names(RFp3Total) == "yhat")] <- paste(colnames(RFp3Total)[2],colnames(RFp3Total)[1],sep="")
    colnames(RFp4Total)[which(names(RFp4Total) == "yhat")] <- paste(colnames(RFp4Total)[2],colnames(RFp4Total)[1],sep="")
    colnames(RFp5Total)[which(names(RFp5Total) == "yhat")] <- paste(colnames(RFp5Total)[2],colnames(RFp5Total)[1],sep="")
    colnames(RFp6Total)[which(names(RFp6Total) == "yhat")] <- paste(colnames(RFp6Total)[2],colnames(RFp6Total)[1],sep="")
    colnames(RFp7Total)[which(names(RFp7Total) == "yhat")] <- paste(colnames(RFp7Total)[2],colnames(RFp7Total)[1],sep="")
    colnames(RFp8Total)[which(names(RFp8Total) == "yhat")] <- paste(colnames(RFp8Total)[2],colnames(RFp8Total)[1],sep="")
    RFTotal <- bind_cols(RFp2Total,RFp4Total,RFp5Total,RFp7Total,RFp8Total)
    write.table(RFTotal,paste(species,"RFPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(RFp6Total,paste(species,"RFBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    #Create 2d histograms with best-fit splines for the partial response curves.
    RFTotal <- read.table(paste(species,"RFPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    RFp6Total <- read.table(paste(species,"RFBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    #RFp1Plot <- ggplot(RFTotal, aes(x=DEM10mAligned, y=yhatDEM10mAligned) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatDEM10mAligned, fill=yhatDEM10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp2Plot <- ggplot(RFTotal, aes(x=Freshwater10mAligned, y=yhatFreshwater10mAligned) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatFreshwater10mAligned, fill=yhatFreshwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    #RFp3Plot <- ggplot(RFTotal, aes(x=Saltwater10mAligned, y=yhatSaltwater10mAligned) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSaltwater10mAligned, fill=yhatSaltwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp4Plot <- ggplot(RFTotal, aes(x=LogSI10mAligned, y=yhatLogSI10mAligned) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatLogSI10mAligned, fill=yhatLogSI10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp5Plot <- ggplot(RFTotal, aes(x=Slope10mAligned, y=yhatSlope10mAligned) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSlope10mAligned, fill=yhatSlope10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp6Plot <- ggplot(RFp6Total, aes(x=as.factor(SoCalBeachType10mAligned),y=yhatSoCalBeachType10mAligned))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
    RFp7Plot <- ggplot(RFTotal, aes(x=SoCalBeachWidth10mAligned, y=yhatSoCalBeachWidth10mAligned) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSoCalBeachWidth10mAligned, fill=yhatSoCalBeachWidth10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp8Plot <- ggplot(RFTotal, aes(x=SVF10mAligned, y=yhatSVF10mAligned) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSVF10mAligned, fill=yhatSVF10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFPlots <- plot_grid(RFp2Plot,RFp4Plot,RFp5Plot,RFp6Plot,RFp7Plot,RFp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Plover"){
    colnames(RFp1Total)[which(names(RFp1Total) == "yhat")] <- paste(colnames(RFp1Total)[2],colnames(RFp1Total)[1],sep="")
    colnames(RFp2Total)[which(names(RFp2Total) == "yhat")] <- paste(colnames(RFp2Total)[2],colnames(RFp2Total)[1],sep="")
    colnames(RFp3Total)[which(names(RFp3Total) == "yhat")] <- paste(colnames(RFp3Total)[2],colnames(RFp3Total)[1],sep="")
    colnames(RFp4Total)[which(names(RFp4Total) == "yhat")] <- paste(colnames(RFp4Total)[2],colnames(RFp4Total)[1],sep="")
    colnames(RFp5Total)[which(names(RFp5Total) == "yhat")] <- paste(colnames(RFp5Total)[2],colnames(RFp5Total)[1],sep="")
    colnames(RFp6Total)[which(names(RFp6Total) == "yhat")] <- paste(colnames(RFp6Total)[2],colnames(RFp6Total)[1],sep="")
    colnames(RFp7Total)[which(names(RFp7Total) == "yhat")] <- paste(colnames(RFp7Total)[2],colnames(RFp7Total)[1],sep="")
    colnames(RFp8Total)[which(names(RFp8Total) == "yhat")] <- paste(colnames(RFp8Total)[2],colnames(RFp8Total)[1],sep="")
    RFTotal <- bind_cols(RFp1Total,RFp2Total,RFp3Total,RFp4Total,RFp5Total,RFp7Total,RFp8Total)
    write.table(RFTotal,paste(species,"RFPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(RFp6Total,paste(species,"RFBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    #Create 2d histograms with best-fit splines for the partial response curves.
    RFTotal <- read.table(paste(species,"RFPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    RFp6Total <- read.table(paste(species,"RFBeachCategoryResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
    RFp1Plot <- ggplot(RFTotal, aes(x=DEM10mAligned, y=yhatDEM10mAligned) )+xlab("Elevation (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatDEM10mAligned, fill=yhatDEM10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp2Plot <- ggplot(RFTotal, aes(x=Freshwater10mAligned, y=yhatFreshwater10mAligned) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatFreshwater10mAligned, fill=yhatFreshwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp3Plot <- ggplot(RFTotal, aes(x=Saltwater10mAligned, y=yhatSaltwater10mAligned) )+xlab("Distance to Saltwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSaltwater10mAligned, fill=yhatSaltwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp4Plot <- ggplot(RFTotal, aes(x=LogSI10mAligned, y=yhatLogSI10mAligned) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatLogSI10mAligned, fill=yhatLogSI10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp5Plot <- ggplot(RFTotal, aes(x=Slope10mAligned, y=yhatSlope10mAligned) )+xlab("Slope (%)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSlope10mAligned, fill=yhatSlope10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp6Plot <- ggplot(RFp6Total, aes(x=as.factor(SoCalBeachType10mAligned),y=yhatSoCalBeachType10mAligned))+geom_boxplot(notch=FALSE)+xlab("Beach category")+ylab("Detection\nProbability")+theme_bw(base_size=25)
    RFp7Plot <- ggplot(RFTotal, aes(x=SoCalBeachWidth10mAligned, y=yhatSoCalBeachWidth10mAligned) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSoCalBeachWidth10mAligned, fill=yhatSoCalBeachWidth10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFp8Plot <- ggplot(RFTotal, aes(x=SVF10mAligned, y=yhatSVF10mAligned) )+xlab("SVF")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSVF10mAligned, fill=yhatSVF10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
    RFPlots <- plot_grid(RFp1Plot,RFp2Plot,RFp3Plot,RFp4Plot,RFp5Plot,RFp6Plot,RFp7Plot,RFp8Plot,ncol=2,labels="AUTO")
  }
  png(paste(species,"GLMPlots.png",sep=""),width=2*800,height=400*length(env.files))
  GLMPlots
  print(GLMPlots)
  dev.off()
  png(paste(species,"XMPlots.png",sep=""),width=2*800,height=400*length(env.files))
  XMPlots
  print(XMPlots)
  dev.off()
  png(paste(species,"RFPlots.png",sep=""),width=2*800,height=400*length(env.files))
  RFPlots
  print(RFPlots)
  dev.off()
}
