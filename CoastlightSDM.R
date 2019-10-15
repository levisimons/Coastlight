rm(list=ls())
require(sp)
require(raster)
require(maptools)
require(rgdal)
require(dismo)
require(rJava)
require(arm)
require(dplyr)
require(sf)
require(ENMeval)
require(randomForest)
require(caret)
require(pdp)

#wd <- "~/Desktop/Coastlight/SDM"
wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
setwd(wd)

#Set species type for observation data.
species <- "Grunion"
#species <- "Plover"

#Set output file.  Send all screen output to this file.
outputFile <- paste("SDMTest",species,".txt",sep="")
sink(outputFile)

# Read in grunion or plover observation points.
if(species=="Grunion"){
  obs.data <- read.csv(file="RandomGrunionPointsWGS84.csv")
  # Read in environmental map layers.
  env.files <- list.files(pattern="Aligned.tif$",full.names=TRUE)
  env.files <- env.files[env.files != "./DEM2mRasterAligned.tif" & env.files != "./DistanceToSaltwaterClipAligned.tif"]
}
if(species=="Plover"){
  obs.data <- read.csv(file="RandomPloverPointsWGS84.csv")
  # Read in environmental map layers.
  env.files <- list.files(pattern="Aligned.tif$",full.names=TRUE)
}

# Drop unused column
obs.data <- obs.data[, c("xcoord", "ycoord")]
set.seed(0)
#obs.data <- sample_n(obs.data,25)

# Read in environmental map layers.
#env.files <- list.files(pattern="Aligned.tif$",full.names=TRUE)
# Make map stack of environmental layers.
env.data <- stack(c(env.files))
# Initialize data containing environmental layer values at presence data locations.
presvals <- obs.data
for(env.file in env.files){
  # Extract environmental map layer values at observation points.
  tmp <- as.data.frame(extract(raster(env.file),obs.data))
  # Get layer names for data frame.
  env.filename <- gsub("^./","",gsub(".tif","",env.file))
  colnames(tmp) <- env.filename
  presvals <- cbind(presvals,tmp)
}
# Standardize missing data
presvals[is.na(presvals)] <- NA
presvals <- na.omit(presvals)

#Record the locations of where full data is available.
obs.data <- presvals[,c("xcoord","ycoord")]

presvals <- presvals[, !names(presvals) %in% c("xcoord","ycoord")]

# setting random seed to always create the same
# random set of points for this example
#Create a random point cloud for pseudo-absences.
coastalArea <- read_sf("DEM2mBoundaryCleaned.shp",quiet=TRUE)
backgr <- as.data.frame(st_coordinates(st_sample(x=coastalArea,size=50*nrow(obs.data),type="random",exact=TRUE)))
absvals <- extract(env.data, backgr)
absvals[is.na(absvals)] <- NA
absvals <- as.data.frame(absvals)
absvals <- cbind(backgr,absvals)
absvals <- na.omit(absvals)
backgr <- absvals[,c("X","Y")]
colnames(backgr) <- c("xcoord","ycoord")
print(paste("Number of presence points:",nrow(obs.data),"Number of background points:",nrow(backgr)))
absvals <- absvals[, !names(absvals) %in% c("X","Y")]

#Merge pseudo-absences with true presences.
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

#Drop factor-type environmental data layers for some of the downstream models.
sdmdata[,"SoCalBeachTypeAligned"] <- factor(sdmdata[,"SoCalBeachTypeAligned"],levels=1:6)
pred_nf <- dropLayer(env.data,"SoCalBeachTypeAligned")

#Construct a training and testing set for the presence data.
group <- kfold(obs.data,5)
pres_train <- obs.data[group!=1,]
pres_test <- obs.data[group==1,]

#Construct a training and testing set for the pseudo-absence data.
group <- kfold(backgr,5)
backgr_train <- backgr[group!=1,]
backgr_test <- backgr[group==1,]

#Construct presence / pseudo-absence training sets.
train <- rbind(pres_train,backgr_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backgr_train)))
envtrain <- extract(env.data,train)
envtrain <- data.frame(cbind(pa=pb_train,envtrain))
envtrain$SoCalBeachTypeAligned <- factor(envtrain$SoCalBeachTypeAligned,levels=1:6)
testpres <- data.frame(extract(env.data,pres_test))
testbackgr <- data.frame(extract(env.data,backgr_test))
testpres$SoCalBeachTypeAligned <- factor(testpres$SoCalBeachTypeAligned,levels=1:6)
testbackgr$SoCalBeachTypeAligned <- factor(testbackgr$SoCalBeachTypeAligned,levels=1:6)

#Random forest model
rf1 <- suppressWarnings(tuneRF(envtrain[,2:ncol(envtrain)],envtrain[,1],plot=FALSE,doBest=TRUE))
print(paste("Random forest variable importance",species,"data"))
importance(rf1)
print(paste("Random forest model evaluation",species,"data"))
erf <- suppressWarnings(evaluate(testpres,testbackgr,rf1))
erf
print(paste("Significance of correlation:",suppressWarnings(erf@pcor)))
print(paste("Significance of AUC:",suppressWarnings(erf@pauc)))
print(paste("Maximum Cohen's kappa score for model:",max(erf@kappa)))
# Calculate Yule's Q.
tmp <- erf@OR
tmp[!is.finite(tmp)] <- NA 
print(paste("Yule's Q score:",(max(tmp,na.rm=1)-1)/(max(tmp,na.rm=1)+1)))
#Variable response functions for random forest model.
if(species=="Grunion"){
  dev.off()
  png(paste("RandomForestResponse",species,".png",sep=""),width=2000,height=2000,res=300)
  #p1 <- partial(rf1,pred.var = "DEM2mRasterAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Elevation (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"DEM2mRasterAligned"]),max(envtrain[which(envtrain$pa==1),"DEM2mRasterAligned"])))
  p2 <- partial(rf1,pred.var = "DistanceToFreshwaterClipAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Distance to freshwater (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"DistanceToFreshwaterClipAligned"]),max(envtrain[which(envtrain$pa==1),"DistanceToFreshwaterClipAligned"])))
  #p3 <- partial(rf1,pred.var = "DistanceToSaltwaterClipAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Distance to saltwater (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"DistanceToSaltwaterClipAligned"]),max(envtrain[which(envtrain$pa==1),"DistanceToSaltwaterClipAligned"])))
  p4 <- partial(rf1,pred.var = "LogSIAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Log(SI) log(mlx)",xlim=c(min(envtrain[which(envtrain$pa==1),"LogSIAligned"]),max(envtrain[which(envtrain$pa==1),"LogSIAligned"])))
  p5 <- partial(rf1,pred.var = "Slope2mRasterAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Slope (%)",xlim=c(min(envtrain[which(envtrain$pa==1),"Slope2mRasterAligned"]),max(envtrain[which(envtrain$pa==1),"Slope2mRasterAligned"])))
  p6 <- partial(rf1,pred.var = "SoCalBeachTypeAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Beach category")
  p7 <- partial(rf1,pred.var = "SoCalBeachWidthAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Beach width (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"SoCalBeachWidthAligned"]),max(envtrain[which(envtrain$pa==1),"SoCalBeachWidthAligned"])))
  p8 <- partial(rf1,pred.var = "SVF2mAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="SVF",xlim=c(min(envtrain[which(envtrain$pa==1),"SVF2mAligned"]),max(envtrain[which(envtrain$pa==1),"SVF2mAligned"])))
  grid.arrange(p2,p4,p5,p6,p7,p8,ncol=2)
  dev.off()
}
if(species=="Plover"){
  dev.off()
  png(paste("RandomForestResponse",species,".png",sep=""),width=2000,height=2000,res=300)
  p1 <- partial(rf1,pred.var = "DEM2mRasterAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Elevation (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"DEM2mRasterAligned"]),max(envtrain[which(envtrain$pa==1),"DEM2mRasterAligned"])))
  p2 <- partial(rf1,pred.var = "DistanceToFreshwaterClipAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Distance to freshwater (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"DistanceToFreshwaterClipAligned"]),max(envtrain[which(envtrain$pa==1),"DistanceToFreshwaterClipAligned"])))
  p3 <- partial(rf1,pred.var = "DistanceToSaltwaterClipAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Distance to saltwater (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"DistanceToSaltwaterClipAligned"]),max(envtrain[which(envtrain$pa==1),"DistanceToSaltwaterClipAligned"])))
  p4 <- partial(rf1,pred.var = "LogSIAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Log(SI) log(mlx)",xlim=c(min(envtrain[which(envtrain$pa==1),"LogSIAligned"]),max(envtrain[which(envtrain$pa==1),"LogSIAligned"])))
  p5 <- partial(rf1,pred.var = "Slope2mRasterAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Slope (%)",xlim=c(min(envtrain[which(envtrain$pa==1),"Slope2mRasterAligned"]),max(envtrain[which(envtrain$pa==1),"Slope2mRasterAligned"])))
  p6 <- partial(rf1,pred.var = "SoCalBeachTypeAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Beach category")
  p7 <- partial(rf1,pred.var = "SoCalBeachWidthAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="Beach width (m)",xlim=c(min(envtrain[which(envtrain$pa==1),"SoCalBeachWidthAligned"]),max(envtrain[which(envtrain$pa==1),"SoCalBeachWidthAligned"])))
  p8 <- partial(rf1,pred.var = "SVF2mAligned",train=envtrain[,-c(1)]) %>% plotPartial(smooth=TRUE,lwd=2,ylab="predicted value",xlab="SVF",xlim=c(min(envtrain[which(envtrain$pa==1),"SVF2mAligned"]),max(envtrain[which(envtrain$pa==1),"SVF2mAligned"])))
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=2)
  dev.off()
}

#Generalized linear model
m1 <- glm(factor(pa) ~ ., data=envtrain,family = binomial(link = "logit"))
print(paste("Generalized linear model variable importance",species,"data"))
varImp(m1,scale=TRUE)
print(paste("Generalized linear model evaluation",species,"data"))
em1 <- suppressWarnings(evaluate(testpres,testbackgr,m1))
em1
print(paste("Significance of correlation:",suppressWarnings(em1@pcor)))
print(paste("Significance of AUC:",suppressWarnings(em1@pauc)))
print(paste("Maximum Cohen's kappa score for model:",max(em1@kappa)))
# Calculate Yule's Q.
tmp <- em1@OR
tmp[!is.finite(tmp)] <- NA 
print(paste("Yule's Q score:",(max(tmp,na.rm=1)-1)/(max(tmp,na.rm=1)+1)))
print("Variable coefficients generalized linear model.")
m1

#MaxEnt model
maxent()
xm <- maxent(x=env.data,p=obs.data,a=backgr,factors='SoCalBeachTypeAligned')
MaxentOutput <- var.importance(xm)
print(paste("Maximum entropy model variable importance",species,"data"))
MaxentOutput
print(paste("Maximum entropy model evaluation",species,"data"))
xmEval <- suppressWarnings(evaluate(testpres,testbackgr,xm))
xmEval
print(paste("Significance of correlation:",suppressWarnings(xmEval@pcor)))
print(paste("Significance of AUC:",suppressWarnings(xmEval@pauc)))
print(paste("Maximum Cohen's kappa score for model:",max(xmEval@kappa)))
# Calculate Yule's Q.
tmp <- xmEval@OR
tmp[!is.finite(tmp)] <- NA 
print(paste("Yule's Q score:",(max(tmp,na.rm=1)-1)/(max(tmp,na.rm=1)+1)))
#Variable response functions for Maxent model.
if(species=="Grunion"){
  dev.off()
  png(paste("MaxentResponsePlot",species,".png",sep=""),width=2000,height=2000,res=300)
  #p1 <- response(xm,var="DEM2mRasterAligned",expand=0,range="p")
  p2 <- response(xm,var="DistanceToFreshwaterClipAligned",expand=0,range="p")
  #p3 <- response(xm,var="DistanceToSaltwaterClipAligned",expand=0,range="p")
  p4 <- response(xm,var="LogSIAligned",expand=0,range="p")
  p5 <- response(xm,var="Slope2mRasterAligned",expand=0,range="p")
  p6 <- response(xm,var="SoCalBeachTypeAligned",expand=0,range="p")
  p7 <- response(xm,var="SoCalBeachWidthAligned",expand=0,range="p")
  p8 <- response(xm,var="SVF2mAligned",expand=0,range="p")
  par(mar=c(4,4,4,4))
  par(mfrow=c(3,2))
  #scatter.smooth(p1[,1],p1[,2],xlab="Elevation (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p2[,1],p2[,2],xlab="Distance to\nfreshwater (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  #scatter.smooth(p3[,1],p3[,2],xlab="Distance to\nsaltwater (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p4[,1],p4[,2],xlab="Log(SI) log(mlx)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p5[,1],p5[,2],xlab="Slope (%)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  plot(p6,col="royalblue1",lwd=2,xlab="Beach category",ylab="predicted value")
  scatter.smooth(p7[,1],p7[,2],xlab="Beach width (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p8[,1],p8[,2],xlab="SVF",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  dev.off()
}
if(species=="Plover"){
  dev.off()
  png(paste("MaxentResponsePlot",species,".png",sep=""),width=2000,height=2000,res=300)
  p1 <- response(xm,var="DEM2mRasterAligned",expand=0,range="p")
  p2 <- response(xm,var="DistanceToFreshwaterClipAligned",expand=0,range="p")
  p3 <- response(xm,var="DistanceToSaltwaterClipAligned",expand=0,range="p")
  p4 <- response(xm,var="LogSIAligned",expand=0,range="p")
  p5 <- response(xm,var="Slope2mRasterAligned",expand=0,range="p")
  p6 <- response(xm,var="SoCalBeachTypeAligned",expand=0,range="p")
  p7 <- response(xm,var="SoCalBeachWidthAligned",expand=0,range="p")
  p8 <- response(xm,var="SVF2mAligned",expand=0,range="p")
  par(mar=c(4,4,4,4))
  par(mfrow=c(4,2))
  scatter.smooth(p1[,1],p1[,2],xlab="Elevation (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p2[,1],p2[,2],xlab="Distance to\nfreshwater (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p3[,1],p3[,2],xlab="Distance to\nsaltwater (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p4[,1],p4[,2],xlab="Log(SI) log(mlx)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p5[,1],p5[,2],xlab="Slope (%)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  plot(p6,col="royalblue1",lwd=2,xlab="Beach category",ylab="predicted value")
  scatter.smooth(p7[,1],p7[,2],xlab="Beach width (m)",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  scatter.smooth(p8[,1],p8[,2],xlab="SVF",ylab="predicted value",type="l",lwd=2,col="black",lpars=c(type="l",col="royalblue1",lwd=2))
  dev.off()
}

sink()
