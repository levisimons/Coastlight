require(dplyr)
require(ggplot2)
require(leaflet)
require(mapview)
require(sp)
require(viridis)
##This script reads in monthly grunion spawning data and finds the spawning event for a given
##site which has the highest Walker scale value for the year.
##It then maps, by year, the beach areas of each grunion run and color codes them by
##their maximum annual Walker score.

wd <- "~/Desktop/Coastlight/Grunions"
setwd(wd)

grunionFiles <- list.files(wd,pattern="_monthly_all_months.csv")
allSitesGrunionMax <- data.frame()
for(grunionFile in grunionFiles){
  #Read in grunion spawn files by site.
  grunionSite <- read.table(grunionFile, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
  #Treat non-observations as 0.
  grunionSite[is.na(grunionSite)] <- 0
  #Find the maximum Walker scale event per year.
  grunionMax <- aggregate(Walker ~ YYYY, data = grunionSite, max)
  #Order grunion maximum spawn events by year
  grunionMax <- arrange(grunionMax,YYYY)
  colnames(grunionMax)[which(names(grunionMax) == "YYYY")] <- "Year"
  #Get site name
  siteName <- unlist(strsplit(grunionFile,split='_monthly',fixed=TRUE))[1]
  grunionMax$SiteName <- siteName
  #Combine all of the max spawn events by site and year into a single data frame.
  allSitesGrunionMax <- rbind(allSitesGrunionMax,grunionMax)
}
#Read in gps coordinates of the bounding boxes denoting each beach site.
BoundaryInput <- read.delim("Bounding_Boxes.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Merge in geographic data into grunion run data.
allSitesGrunionMaxGPS <- merge(allSitesGrunionMax,BoundaryInput,by.x=c("SiteName"),by.y=c("Bounding_Box"))

#Subset grunion data by year for mapping.
GivenYear <- 2003 #Between 2003 and 2018
allSitesGrunionMaxGPSYear <- subset(allSitesGrunionMaxGPS,Year==GivenYear)

#Map grunion runs.
MapCoordinates <- allSitesGrunionMaxGPSYear
#Set centroid coordinates of each beach area.
MapCoordinates$lat <- 0.5*(MapCoordinates$lat_min+MapCoordinates$lat_max)
MapCoordinates$lon <- 0.5*(MapCoordinates$lon_min+MapCoordinates$lon_max)
#Map data.
dev.off()
#Color code beach areas by the maximum annual Walker score.
ColorScale <- colorNumeric(palette=plasma(10),domain=allSitesGrunionMaxGPS$Walker)
CalMap = leaflet(MapCoordinates) %>% 
  addTiles() %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  setView(lng=mean(MapCoordinates$lon),lat=mean(MapCoordinates$lat),zoom=6) %>%
  addLegend("topright", pal=ColorScale,values=~Walker,title=paste("Max Walker scale<br>Year",GivenYear))

CalMap
for(site in unique(MapCoordinates$SiteName)){
  lon_min=MapCoordinates[MapCoordinates$SiteName==site,"lon_min"]
  lon_max=MapCoordinates[MapCoordinates$SiteName==site,"lon_max"]
  lat_min=MapCoordinates[MapCoordinates$SiteName==site,"lat_min"]
  lat_max=MapCoordinates[MapCoordinates$SiteName==site,"lat_max"]
  MaxWalker=MapCoordinates[MapCoordinates$SiteName==site,"Walker"]
  CalMap <- CalMap %>% addRectangles(lng1=lon_min,lat1=lat_min,lng2=lon_max,lat2=lat_max,color=~ColorScale(MaxWalker),weight=5,fillColor="transparent")
}
CalMap
