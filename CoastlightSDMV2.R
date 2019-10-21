rm(list=ls())
require(sp)
require(raster)
require(maptools)
require(rgdal)
require(dismo)
require(rJava)
require(arm)
require(dplyr)
require(plyr)
require(sf)
require(ENMeval)
require(randomForest)
require(caret)
require(pdp)
require(ggplot2)

#wd <- "~/Desktop/Coastlight/SDM"
wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
setwd(wd)

#Set species type for observation data.
species <- "Grunion"
#species <- "Plover"

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
  RFEvaluation$Q <- (max(tmp,na.rm=T)-1)/(max(tmp,na.rm=T)+1)
  RFEvaluation$TSS <- max(erf@TPR,na.rm=T)+max(erf@TNR,na.rm=T)-1
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
  m1 <- glm(pa ~ ., data=envtrain[,-c(2,3)],family = binomial(link = "logit"))
  GLMImportance <- varImp(m1,scale=TRUE)
  GLMImportance <- data.frame(names=row.names(GLMImportance),GLMImportance)
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
  GLMEvaluation$Q <- (max(tmp,na.rm=T)-1)/(max(tmp,na.rm=T)+1)
  GLMEvaluation$TSS <- max(em1@TPR,na.rm=T)+max(em1@TNR,na.rm=T)-1
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
  XMEvaluation$Q <- (max(tmp,na.rm=T)-1)/(max(tmp,na.rm=T)+1)
  XMEvaluation$TSS <- max(exm@TPR,na.rm=T)+max(exm@TNR,na.rm=T)-1
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
tmpMean <- RFImportanceTotal %>% group_by(names) %>% summarise_all(funs(mean(IncNodePurity)))
tmpSD <- RFImportanceTotal %>% group_by(names) %>% summarise_all(funs(sd(IncNodePurity)))
RFImportanceTotal <- left_join(tmpMean,tmpSD,by=c("names"))
colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
RFImportanceTotal
tmpMean <- colMeans(RFEvaluationTotal)
tmpSD <- apply(RFEvaluationTotal,2,sd)
RFEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
RFEvaluationTotal
#GLM summary statistics
tmpMean <- GLMImportanceTotal %>% group_by(names) %>% summarise_all(funs(mean(Overall)))
tmpSD <- GLMImportanceTotal %>% group_by(names) %>% summarise_all(funs(sd(Overall)))
GLMImportanceTotal <- left_join(tmpMean,tmpSD,by=c("names"))
colnames(GLMImportanceTotal) <- c("Variable","MeanImportance","SDImportance")
tmpMean <- colMeans(GLMEvaluationTotal)
tmpSD <- apply(GLMEvaluationTotal,2,sd)
GLMEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
#Maxent summary statistics
tmpMean <- XMImportanceTotal %>% group_by(variable) %>% summarise_all(funs(mean(permutation.importance)))
tmpSD <- XMImportanceTotal %>% group_by(variable) %>% summarise_all(funs(sd(permutation.importance)))
XMImportanceTotal <- left_join(tmpMean,tmpSD,by=c("variable"))
colnames(XMImportanceTotal) <- c("Variable","MeanImportance","SDImportance")
tmpMean <- colMeans(XMEvaluationTotal)
tmpSD <- apply(XMEvaluationTotal,2,sd)
XMEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))

