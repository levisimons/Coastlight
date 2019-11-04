rm(list=ls())
require(sp)
require(raster)
require(maptools)
require(rgdal)
require(dismo)
require(rJava)
require(arm)
require(plyr)
require(dplyr)
require(sf)
require(ENMeval)
require(randomForest)
require(caret)
require(pdp)
require(ggplot2)
require(cowplot)
require(splines)

.jinit()

#wd <- "~/Desktop/Coastlight/SDM"
wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
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
  tmpMean <- RFImportanceTotal %>% group_by(names) %>% summarise_all(funs(mean(IncNodePurity)))
  tmpSD <- RFImportanceTotal %>% group_by(names) %>% summarise_all(funs(sd(IncNodePurity)))
  RFImportanceTotal <- left_join(tmpMean,tmpSD,by=c("names"))
  colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
  #To save aggregated data frame.
  write.table(RFImportanceTotal,paste(species,"RFImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(RFEvaluationTotal)
  tmpSD <- apply(RFEvaluationTotal,2,sd)
  RFEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(RFEvaluationTotal,paste(species,"RFEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
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
  }
  
  #GLM summary statistics
  tmpMean <- GLMImportanceTotal %>% group_by(names) %>% summarise_all(funs(mean(z.value)))
  tmpSD <- GLMImportanceTotal %>% group_by(names) %>% summarise_all(funs(sd(z.value)))
  GLMImportanceTotal <- left_join(tmpMean,tmpSD,by=c("names"))
  colnames(GLMImportanceTotal) <- c("Variable","MeanImportance","SDImportance")
  #To save aggregated data frame.
  write.table(GLMImportanceTotal,paste(species,"GLMImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(GLMEvaluationTotal)
  tmpSD <- apply(GLMEvaluationTotal,2,sd)
  GLMEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(GLMEvaluationTotal,paste(species,"GLMEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
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
  }
  
  #Maxent summary statistics
  tmpMean <- XMImportanceTotal %>% group_by(variable) %>% summarise_all(funs(mean(permutation.importance)))
  tmpSD <- XMImportanceTotal %>% group_by(variable) %>% summarise_all(funs(sd(permutation.importance)))
  XMImportanceTotal <- left_join(tmpMean,tmpSD,by=c("variable"))
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
    #colnames(XMp1Total)[which(names(XMp1Total) == "yhat")] <- paste(colnames(XMp1Total)[2],colnames(XMp1Total)[1],sep="")
    colnames(XMp2Total)[which(names(XMp2Total) == "yhat")] <- paste(colnames(XMp2Total)[2],colnames(XMp2Total)[1],sep="")
    #colnames(XMp3Total)[which(names(XMp3Total) == "yhat")] <- paste(colnames(XMp3Total)[2],colnames(XMp3Total)[1],sep="")
    colnames(XMp4Total)[which(names(XMp4Total) == "yhat")] <- paste(colnames(XMp4Total)[2],colnames(XMp4Total)[1],sep="")
    colnames(XMp5Total)[which(names(XMp5Total) == "yhat")] <- paste(colnames(XMp5Total)[2],colnames(XMp5Total)[1],sep="")
    colnames(XMp6Total)[which(names(XMp6Total) == "yhat")] <- paste(colnames(XMp6Total)[2],colnames(XMp6Total)[1],sep="")
    colnames(XMp7Total)[which(names(XMp7Total) == "yhat")] <- paste(colnames(XMp7Total)[2],colnames(XMp7Total)[1],sep="")
    colnames(XMp8Total)[which(names(XMp8Total) == "yhat")] <- paste(colnames(XMp8Total)[2],colnames(XMp8Total)[1],sep="")
    GLMTotal <- bind_cols(XMp2Total,XMp4Total,XMp5Total,XMp7Total,XMp8Total)
    write.table(XMTotal,paste(species,"XMPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(XMp6Total,paste(species,"XMBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  }
  if(species=="Plover"){
    colnames(XMp1Total)[which(names(XMp1Total) == "yhat")] <- paste(colnames(XMp1Total)[2],colnames(XMp1Total)[1],sep="")
    colnames(XMp2Total)[which(names(XMp2Total) == "yhat")] <- paste(colnames(XMp2Total)[2],colnames(XMp2Total)[1],sep="")
    colnames(XMp3Total)[which(names(XMp3Total) == "yhat")] <- paste(colnames(XMp3Total)[2],colnames(XMp3Total)[1],sep="")
    colnames(XMp4Total)[which(names(XMp4Total) == "yhat")] <- paste(colnames(XMp4Total)[2],colnames(XMp4Total)[1],sep="")
    colnames(XMp5Total)[which(names(XMp5Total) == "yhat")] <- paste(colnames(XMp5Total)[2],colnames(XMp5Total)[1],sep="")
    colnames(XMp6Total)[which(names(XMp6Total) == "yhat")] <- paste(colnames(XMp6Total)[2],colnames(XMp6Total)[1],sep="")
    colnames(XMp7Total)[which(names(XMp7Total) == "yhat")] <- paste(colnames(XMp7Total)[2],colnames(XMp7Total)[1],sep="")
    colnames(XMp8Total)[which(names(XMp8Total) == "yhat")] <- paste(colnames(XMp8Total)[2],colnames(XMp8Total)[1],sep="")
    GLMTotal <- bind_cols(XMp1Total,XMp2Total,XMp3Total,XMp4Total,XMp5Total,XMp7Total,XMp8Total)
    write.table(XMTotal,paste(species,"XMPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    write.table(XMp6Total,paste(species,"XMBeachCategoryResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  }
  
  #Partial dependence plots for the random forest model
  if(species=="Plover"){
    RFp1Plot <- ggplot(RFp1Total, aes(x=DEM10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp1Plot <- RFp1Plot+xlab("Elevation (m)")+ylab("Detection\nProbability")
    RFp2Plot <- ggplot(RFp2Total, aes(x=Freshwater10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp2Plot <- RFp2Plot+xlab("Distance to freshwater (m)")+ylab("Detection\nProbability")
    RFp3Plot <- ggplot(RFp3Total, aes(x=Saltwater10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp3Plot <- RFp3Plot+xlab("Distance to saltwater (m)")+ylab("Detection\nProbability")
    RFp4Plot <- ggplot(RFp4Total, aes(x=LogSI10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp4Plot <- RFp4Plot+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")
    RFp5Plot <- ggplot(RFp5Total, aes(x=Slope10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp5Plot <- RFp5Plot+xlab("Slope (%)")+ylab("Detection\nProbability")
    RFp6Plot <- ggplot(RFp6Total, aes(x=SoCalBeachType10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_boxplot(notch=FALSE)
    RFp6Plot <- RFp6Plot+xlab("Beach category")+ylab("Detection\nProbability")
    RFp7Plot <- ggplot(RFp7Total, aes(x=SoCalBeachWidth10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp7Plot <- RFp7Plot+xlab("Beach width (m)")+ylab("Detection\nProbability")
    RFp8Plot <- ggplot(RFp8Total, aes(x=SVF10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp8Plot <- RFp8Plot+xlab("SVF")+ylab("Detection\nProbability")
    RFPlots <- plot_grid(RFp1Plot,RFp2Plot,RFp3Plot,RFp4Plot,RFp5Plot,RFp6Plot,RFp7Plot,RFp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Grunion"){
    #RFp1Plot <- ggplot(RFp1Total, aes(x=DEM10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    #RFp1Plot <- RFp1Plot+xlab("Elevation (m)")+ylab("Detection\nProbability")
    RFp2Plot <- ggplot(RFp2Total, aes(x=Freshwater10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp2Plot <- RFp2Plot+xlab("Distance to freshwater (m)")+ylab("Detection\nProbability")
    #RFp3Plot <- ggplot(RFp3Total, aes(x=Saltwater10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    #RFp3Plot <- RFp3Plot+xlab("Distance to saltwater (m)")+ylab("Detection\nProbability")
    RFp4Plot <- ggplot(RFp4Total, aes(x=LogSI10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp4Plot <- RFp4Plot+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")
    RFp5Plot <- ggplot(RFp5Total, aes(x=Slope10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp5Plot <- RFp5Plot+xlab("Slope (%)")+ylab("Detection\nProbability")
    RFp6Plot <- ggplot(RFp6Total, aes(x=SoCalBeachType10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_boxplot(notch=FALSE)
    RFp6Plot <- RFp6Plot+xlab("Beach category")+ylab("Detection\nProbability")
    RFp7Plot <- ggplot(RFp7Total, aes(x=SoCalBeachWidth10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp7Plot <- RFp7Plot+xlab("Beach width (m)")+ylab("Detection\nProbability")
    RFp8Plot <- ggplot(RFp8Total, aes(x=SVF10mAligned,y=yhat))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    RFp8Plot <- RFp8Plot+xlab("SVF")+ylab("Detection\nProbability")
    RFPlots <- plot_grid(RFp2Plot,RFp4Plot,RFp5Plot,RFp6Plot,RFp7Plot,RFp8Plot,ncol=2,labels="AUTO")
  }
  png(paste(species,"RFPlots.png",sep=""),width=2*400,height=200*length(env.files))
  RFPlots
  dev.off()
  
  #Partial dependence plots for the GLM model
  if(species=="Plover"){
    GLMp1Plot <- ggplot(GLMp1Total, aes(x=DEM10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp1Plot <- GLMp1Plot+xlab("Elevation (m)")+ylab("Detection\nProbability")
    GLMp2Plot <- ggplot(GLMp2Total, aes(x=Freshwater10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp2Plot <- GLMp2Plot+xlab("Distance to freshwater (m)")+ylab("Detection\nProbability")
    GLMp3Plot <- ggplot(GLMp3Total, aes(x=Saltwater10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp3Plot <- GLMp3Plot+xlab("Distance to saltwater (m)")+ylab("Detection\nProbability")
    GLMp4Plot <- ggplot(GLMp4Total, aes(x=LogSI10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp4Plot <- GLMp4Plot+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")
    GLMp5Plot <- ggplot(GLMp5Total, aes(x=Slope10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp5Plot <- GLMp5Plot+xlab("Slope (%)")+ylab("Detection\nProbability")
    GLMp6Plot <- ggplot(GLMp6Total, aes(x=SoCalBeachType10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_boxplot(notch=FALSE)
    GLMp6Plot <- GLMp6Plot+xlab("Beach category")+ylab("Detection\nProbability")
    GLMp7Plot <- ggplot(GLMp7Total, aes(x=SoCalBeachWidth10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp7Plot <- GLMp7Plot+xlab("Beach width (m)")+ylab("Detection\nProbability")
    GLMp8Plot <- ggplot(GLMp8Total, aes(x=SVF10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp8Plot <- GLMp8Plot+xlab("SVF")+ylab("Detection\nProbability")
    GLMPlots <- plot_grid(GLMp1Plot,GLMp2Plot,GLMp3Plot,GLMp4Plot,GLMp5Plot,GLMp6Plot,GLMp7Plot,GLMp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Grunion"){
    #GLMp1Plot <- ggplot(GLMp1Total, aes(x=DEM10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    #GLMp1Plot <- GLMp1Plot+xlab("Elevation (m)")+ylab("Detection\nProbability")
    GLMp2Plot <- ggplot(GLMp2Total, aes(x=Freshwater10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp2Plot <- GLMp2Plot+xlab("Distance to freshwater (m)")+ylab("Detection\nProbability")
    #GLMp3Plot <- ggplot(GLMp3Total, aes(x=Saltwater10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    #GLMp3Plot <- GLMp3Plot+xlab("Distance to saltwater (m)")+ylab("Detection\nProbability")
    GLMp4Plot <- ggplot(GLMp4Total, aes(x=LogSI10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp4Plot <- GLMp4Plot+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")
    GLMp5Plot <- ggplot(GLMp5Total, aes(x=Slope10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp5Plot <- GLMp5Plot+xlab("Slope (%)")+ylab("Detection\nProbability")
    GLMp6Plot <- ggplot(GLMp6Total, aes(x=SoCalBeachType10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_boxplot(notch=FALSE)
    GLMp6Plot <- GLMp6Plot+xlab("Beach category")+ylab("Detection\nProbability")
    GLMp7Plot <- ggplot(GLMp7Total, aes(x=SoCalBeachWidth10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp7Plot <- GLMp7Plot+xlab("Beach width (m)")+ylab("Detection\nProbability")
    GLMp8Plot <- ggplot(GLMp8Total, aes(x=SVF10mAligned,y=exp(yhat)/(1+exp(yhat))))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=yhat))
    GLMp8Plot <- GLMp8Plot+xlab("SVF")+ylab("Detection\nProbability")
    GLMPlots <- plot_grid(GLMp2Plot,GLMp4Plot,GLMp5Plot,GLMp6Plot,GLMp7Plot,GLMp8Plot,ncol=2,labels="AUTO")
  }
  png(paste(species,"GLMPlots.png",sep=""),width=2*400,height=200*length(env.files))
  GLMPlots
  dev.off()
  
  #Partial dependence plots for the Maxent model
  if(species=="Plover"){
    XMp1Plot <- ggplot(XMp1Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp1Plot <- XMp1Plot+xlab("Elevation (m)")+ylab("Detection\nProbability")
    XMp2Plot <- ggplot(XMp2Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp2Plot <- XMp2Plot+xlab("Distance to freshwater (m)")+ylab("Detection\nProbability")
    XMp3Plot <- ggplot(XMp3Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp3Plot <- XMp3Plot+xlab("Distance to saltwater (m)")+ylab("Detection\nProbability")
    XMp4Plot <- ggplot(XMp4Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp4Plot <- XMp4Plot+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")
    XMp5Plot <- ggplot(XMp5Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp5Plot <- XMp5Plot+xlab("Slope (%)")+ylab("Detection\nProbability")
    XMp6Plot <- ggplot(XMp6Total, aes(x=as.factor(V1),y=p))+theme(text = element_text(size=25))+geom_boxplot(notch=FALSE)
    XMp6Plot <- XMp6Plot+xlab("Beach category")+ylab("Detection\nProbability")
    XMp7Plot <- ggplot(XMp7Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp7Plot <- XMp7Plot+xlab("Beach width (m)")+ylab("Detection\nProbability")
    XMp8Plot <- ggplot(XMp8Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp8Plot <- XMp8Plot+xlab("SVF")+ylab("Detection\nProbability")
    XMPlots <- plot_grid(XMp1Plot,XMp2Plot,XMp3Plot,XMp4Plot,XMp5Plot,XMp6Plot,XMp7Plot,XMp8Plot,ncol=2,labels="AUTO")
  }
  if(species=="Grunion"){
    #XMp1Plot <- ggplot(XMp1Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    #XMp1Plot <- XMp1Plot+xlab("Elevation (m)")+ylab("Detection\nProbability")
    XMp2Plot <- ggplot(XMp2Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp2Plot <- XMp2Plot+xlab("Distance to freshwater (m)")+ylab("Detection\nProbability")
    #XMp3Plot <- ggplot(XMp3Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    #XMp3Plot <- XMp3Plot+xlab("Distance to saltwater (m)")+ylab("Detection\nProbability")
    XMp4Plot <- ggplot(XMp4Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp4Plot <- XMp4Plot+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")
    XMp5Plot <- ggplot(XMp5Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp5Plot <- XMp5Plot+xlab("Slope (%)")+ylab("Detection\nProbability")
    XMp6Plot <- ggplot(XMp6Total, aes(x=as.factor(V1),y=p))+theme(text = element_text(size=25))+geom_boxplot(notch=FALSE)
    XMp6Plot <- XMp6Plot+xlab("Beach category")+ylab("Detection\nProbability")
    XMp7Plot <- ggplot(XMp7Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp7Plot <- XMp7Plot+xlab("Beach width (m)")+ylab("Detection\nProbability")
    XMp8Plot <- ggplot(XMp8Total, aes(x=V1,y=p))+theme(text = element_text(size=25))+geom_smooth(method=lm, formula = y ~ bs(x, degree = 2), aes(fill=p))
    XMp8Plot <- XMp8Plot+xlab("SVF")+ylab("Detection\nProbability")
    XMPlots <- plot_grid(XMp2Plot,XMp4Plot,XMp5Plot,XMp6Plot,XMp7Plot,XMp8Plot,ncol=2,labels="AUTO")
  }
  png(paste(species,"XMPlots.png",sep=""),width=2*400,height=200*length(env.files))
  XMPlots
  dev.off()
}
