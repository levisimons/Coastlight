
#This script reads in the beach zones where grunion runs are measured and generates
#the bounding polygon coordinates for importing those areas into QGIS.

wd <- "~/Desktop/Coastlight/Grunions"
setwd(wd)

#Read in gps coordinates of the bounding boxes denoting each beach site.
BoundaryInput <- read.delim("Bounding_Boxes.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Create a text file with the bounding box coordinates for importation into QGIS and
#generating a new polygon layer.
QGISDat <- data.frame()
for(site in unique(BoundaryInput$Bounding_Box)){
  lon_min=BoundaryInput[BoundaryInput$Bounding_Box==site,"lon_min"]
  lon_max=BoundaryInput[BoundaryInput$Bounding_Box==site,"lon_max"]
  lat_min=BoundaryInput[BoundaryInput$Bounding_Box==site,"lat_min"]
  lat_max=BoundaryInput[BoundaryInput$Bounding_Box==site,"lat_max"]
  rowDat <- data.frame()
  rowDat[1,1] <- paste("POLYGON ((",lon_min," ",lat_min,", ",lon_min," ",lat_max,", ",lon_max," ",lat_max,", ",lon_min," ",lat_max,", ",lon_min," ",lat_min,"))",sep="")
  QGISDat <- rbind(QGISDat,rowDat)
}

write.table(QGISDat,"Grunion_Boxes.txt",quote=FALSE,sep="\n",row.names = FALSE,col.names = FALSE)
