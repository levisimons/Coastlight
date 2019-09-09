require(sp)
require(raster)
require(maptools)
require(rgdal)
require(dismo)
require(rJava)
require(arm)
require(dplyr)
require(sf)

#wd <- "~/Desktop/Coastlight/SDM"
wd <- "/home/cmb-07/sn1/alsimons/Coastlight"
setwd(wd)

# Read in grunion observation points.
obs.data <- read.csv(file="RandomGrunionPoints.csv")
# Drop unused column
obs.data <- obs.data[, c("xcoord", "ycoord")]
obs.data <- sample_n(obs.data,500)

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
set.seed(0)
#Create a random point cloud for pseudo-absences.
#backgr <- randomPoints(mask=env.data,n=100*nrow(presvals),p=presvals,ext=extent(env.data),extf=1,excludep=TRUE,tryf=300)
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

#absvals <- absvals[, !names(absvals) %in% c("xcoord","ycoord","SoCalBeachTypeAligned","SoCalBeachWidthAligned")]
#Merge pseudo-absences with true presences.
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

#Drop factor-type environmental data layers for some of the downstream models.
sdmdata[,"SoCalBeachTypeAligned"] <- as.factor(sdmdata[,"SoCalBeachTypeAligned"])
pred_nf <- dropLayer(env.data,"SoCalBeachTypeAligned")

#Construct a training and testing set for the presence data.
group <- kfold(obs.data,5)
pres_train <- obs.data[group!=1,]
pres_test <- obs.data[group==1,]

#Construct a training and testing set for the pseudo-absence data.
group <- kfold(backgr,5)
backgr_train <- backgr[group!=1,]
backgr_test <- backgr[group==1,]

#GLM
#m1 <- glm(pb ~ ., data=sdmdata,family = binomial(link = "logit"))
#summary(m1)

m2 <- bayesglm(pb ~. , data=sdmdata,family = binomial(link = "logit"), prior.scale = 2.5, maxit = 10000)
summary(m2)

#BioClim model
#bc <- bioclim(pred_nf,pres_train)
#e <- evaluate(pres_test,backgr_test,bc,pred_nf)
#print("BioClim model")
#e

#Domain model
#dm <- domain(pred_nf,pres_train)
#e <- evaluate(pres_test,backgr_test,dm,pred_nf)
#print("Domain model")
#e

#Mahalanobis distance model
#mm <- mahal(pred_nf, pres_train)
#e <- evaluate(pres_test, backgr_test, mm, pred_nf)
#print("Mahalanobis model")
#e

#MaxEnt model
maxent()
xm <- maxent(x=env.data,p=obs.data,a=backgr,factors='SoCalBeachTypeAligned')
MaxentOutput <- as.data.frame(xm@results)
MaxentOutput$Parameter <- rownames(MaxentOutput)
colnames(MaxentOutput) <- c("Value","Parameter")
write.table(MaxentOutput,"MaxentOutput.txt",quote=FALSE,sep="\t",row.names = FALSE)
