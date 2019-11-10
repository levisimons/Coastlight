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
#.rs.unloadPackage("tidyr")
#To deal with java issues
.jinit()

wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
#wd <- "~/Desktop/Coastlight/SDM"
setwd(wd)

##This part is to be run first to determine variables to keep for parsimonious random forest models
#Set species types for observation data.
#Set number of model iterations
modelNum <- 3
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
  for(i in 1:modelNums){
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
  }
  #Summary statistics on the frequency and importance of environmental parameters in random forest model.
  tmp <- as.data.frame(table(RFImportanceTotal$names))
  colnames(tmp) <- c("Variable","Freq")
  RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
  colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
  RFImportanceTotal <- left_join(tmp,RFImportanceTotal)
  #To save aggregated data frame.
  write.table(RFImportanceTotal,paste(species,"ParsimonyRFImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
}

#Using the environmental variables found to be included in parsimonious
#random forest models rerun the model iterations.
speciesList <- c("Grunion","Plover")
for(species in speciesList){
  #Select the variables to include in the parsimonious random forest model.
  RFImportanceTotal <- read.table(paste(species,"ParsimonyRFImportance.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8")
  RFVars <- RFImportanceTotal[RFImportanceTotal$Freq==100,c("Variable")]
  # Read in grunion or plover presence and pseudo-absence points.
  if(species=="Grunion"){
    obs.data <- read.csv(file="RandomGrunionPointsWGS84.csv")
    abs.data <- read.csv(file="RandomGrunionAbsencesWGS84.csv")
  }
  if(species=="Plover"){
    obs.data <- read.csv(file="RandomPloverPointsWGS84.csv")
    abs.data <- read.csv(file="RandomPloverAbsencesWGS84.csv")
  }
  # Read in environmental map layers.
  env.files <- list.files(pattern="10mAligned.tif$",full.names=TRUE)
  env.files <- env.files[grep(paste(RFVars,collapse="|"),env.files)]
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
  for(i in 1:modelNum){
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
    #envtrain$SoCalBeachType10mAligned <- factor(envtrain$SoCalBeachType10mAligned,levels=1:6)
    testpres <- pres_test
    testbackgr <- backgr_test
    #testpres$SoCalBeachType10mAligned <- factor(testpres$SoCalBeachType10mAligned,levels=1:6)
    #testbackgr$SoCalBeachType10mAligned <- factor(testbackgr$SoCalBeachType10mAligned,levels=1:6)
    
    #Parsimonious random forest model
    rf.regress <- suppressWarnings(randomForest(envtrain[,4:ncol(envtrain)],envtrain[,c(1)], imp.scale="mir", parsimony=0.03,final.model=TRUE,seed=1))
    RFImportance <- importance(rf.regress)
    RFImportance <- data.frame(names=row.names(RFImportance),RFImportance)
    RFImportanceTotal <- rbind(RFImportanceTotal,RFImportance)
    erf <- suppressWarnings(evaluate(testpres,testbackgr,rf.regress))
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
      RFp1 <- partial(rf.regress,pred.var = "Freshwater10mAligned",train=envtrain[,c(4:ncol(envtrain))])
      RFp1Total <- rbind(RFp1Total,RFp1)
    }
    RFp2 <- partial(rf.regress,pred.var = "LogSI10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp2Total <- rbind(RFp2Total,RFp2)
    RFp3 <- partial(rf.regress,pred.var = "SoCalBeachWidth10mAligned",train=envtrain[,c(4:ncol(envtrain))])
    RFp3Total <- rbind(RFp3Total,RFp3)
  }
  #Summary statistics on the frequency and importance of environmental parameters in random forest model.
  tmp <- as.data.frame(table(RFImportanceTotal$names))
  colnames(tmp) <- c("Variable","Freq")
  RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
  colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
  RFImportanceTotal <- left_join(tmp,RFImportanceTotal)
  #To save aggregated data frame.
  write.table(RFImportanceTotal,paste(species,"ParsimonyFinalRFImportance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  tmpMean <- colMeans(RFEvaluationTotal)
  tmpSD <- apply(RFEvaluationTotal,2,sd)
  RFEvaluationTotal <- as.data.frame(rbind(tmpMean,tmpSD))
  #To save aggregated data frame.
  write.table(RFEvaluationTotal,paste(species,"ParsimonyFinalRFEvaluation.txt",sep=""),quote=FALSE,sep="\t",row.names = TRUE)
  #To save all of the raw points used in the partial response plots for the random forest model.
  if(species=="Plover"){
    colnames(RFp1Total)[which(names(RFp1Total) == "yhat")] <- paste(colnames(RFp1Total)[2],colnames(RFp1Total)[1],sep="")
  }
  colnames(RFp2Total)[which(names(RFp2Total) == "yhat")] <- paste(colnames(RFp2Total)[2],colnames(RFp2Total)[1],sep="")
  colnames(RFp3Total)[which(names(RFp3Total) == "yhat")] <- paste(colnames(RFp3Total)[2],colnames(RFp3Total)[1],sep="")
  if(species=="Plover"){
    RFTotal <- bind_cols(RFp1Total,RFp2Total,RFp3Total)
  }
  if(species=="Grunion"){
    RFTotal <- bind_cols(RFp2Total,RFp3Total)
  }
  write.table(RFTotal,paste(species,"RFParsimonyFinalPartialResponse.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  #Create 2d histograms with best-fit splines for the partial response curves.
  RFTotal <- read.table(paste(species,"RFParsimonyFinalPartialResponse.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  if(species=="Plover"){
    RFp1Plot <- ggplot(RFTotal, aes(x=Freshwater10mAligned, y=yhatFreshwater10mAligned) )+xlab("Distance to Freshwater (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatFreshwater10mAligned, fill=yhatFreshwater10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  }
  RFp2Plot <- ggplot(RFTotal, aes(x=LogSI10mAligned, y=yhatLogSI10mAligned) )+xlab("Log(SI) log(mlx)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatLogSI10mAligned, fill=yhatLogSI10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  RFp3Plot <- ggplot(RFTotal, aes(x=SoCalBeachWidth10mAligned, y=yhatSoCalBeachWidth10mAligned) )+xlab("Beach width (m)")+ylab("Detection\nProbability")+geom_bin2d(bins = 50)+scale_fill_continuous(type = "viridis")+stat_smooth(aes(y = yhatSoCalBeachWidth10mAligned, fill=yhatSoCalBeachWidth10mAligned),method="auto",formula=y~x,color="violet",fill="red",n=0.1*nrow(RFTotal))+theme_bw(base_size=25)
  if(species=="Plover"){
    RFPlots <- plot_grid(RFp1Plot,RFp2Plot,RFp3Plot,ncol=2,labels="AUTO")
  }
  if(species=="Grunion"){
    RFPlots <- plot_grid(RFp2Plot,RFp3Plot,ncol=2,labels="AUTO")
  }
  png(paste(species,"RFParsimonyFinalPlots.png",sep=""),width=2*800,height=400*length(env.files))
  RFPlots
  print(RFPlots)
  dev.off()
}
