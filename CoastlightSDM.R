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
}
if(species=="Plover"){
  obs.data <- read.csv(file="RandomPloverPointsWGS84.csv")
}

# Drop unused column
obs.data <- obs.data[, c("xcoord", "ycoord")]
set.seed(0)
#obs.data <- sample_n(obs.data,200)

# Read in environmental map layers.
env.files <- list.files(pattern="Aligned.tif$",full.names=TRUE)
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
rf1 <- suppressWarnings(tuneRF(envtrain[,2:9],envtrain[,1],plot=FALSE,doBest=TRUE))
print(paste("Random forest variable importance",species,"data"))
importance(rf1)
print(paste("Random forest model evaluation",species,"data"))
erf <- suppressWarnings(evaluate(testpres,testbackgr,rf1))
erf

#Generalized linear model
m1 <- glm(factor(pa) ~ ., data=envtrain,family = binomial(link = "logit"))
print(paste("Generalized linear model variable importance",species,"data"))
varImp(m1,scale=TRUE)
print(paste("Generalized linear model evaluation",species,"data"))
em1 <- suppressWarnings(evaluate(testpres,testbackgr,m1))
em1

#MaxEnt model
maxent()
xm <- maxent(x=env.data,p=obs.data,a=backgr,factors='SoCalBeachTypeAligned')
MaxentOutput <- var.importance(xm)
print(paste("Maximum entropy model variable importance",species,"data"))
MaxentOutput
print(paste("Maximum entropy model evaluation",species,"data"))
xmEval <- suppressWarnings(evaluate(testpres,testbackgr,xm))
xmEval

sink()
#
write.table(MaxentOutput,"MaxentOutput.txt",quote=FALSE,sep="\t",row.names = FALSE)
#Save model graphs
pdf("MaxentResponsePlot.pdf")
response(xm)
dev.off()
pdf("MaxentResponseContributions.pdf")
plot(xm)
dev.off()

#Further investigations of the Maxent model using ENMevaluate
XMeval <- ENMevaluate(obs.data, env.data, backgr, method='randomkfold', kfolds=5, RMvalues=c(1,2), fc=c('L','LP'), parallel=FALSE, algorithm='maxent.jar',rasterPreds=FALSE,categoricals="SoCalBeachTypeAligned")
