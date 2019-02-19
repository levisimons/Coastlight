require(raster)
require(rgdal)
require(ggmap)
require(maps)
require(mapview)
require(mapdata)
require(munsell)
require(leaflet)
require(devtools)
require(webshot)

setwd("~/Desktop/Coastlight/CubeSat")

#Read in raster nighttime image data.
NightRaster <- raster("LA_06_0811__ground100msDarkUsed_16bit.tiff")

#Georeference the image.
xmin(NightRaster) <- -119.15
xmax(NightRaster) <- -117.35
ymin(NightRaster) <- 33.15
ymax(NightRaster) <- 34.85
crs(NightRaster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#Extract brightest raster points and map them.
NightPoints <- rasterToPoints(NightRaster,spatial=TRUE)
NightTable <- as.data.frame(NightPoints)
colnames(NightTable) <- c("Brightness","Longitude","Latitude")
MapCoordinates <- NightTable[which(NightTable$Brightness >=10*median(NightTable$Brightness)),]
CalMap = leaflet(MapCoordinates) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=rainbow(10),domain=MapCoordinates$Brightness)
CalMap %>% addCircleMarkers(color = ~ColorScale(Brightness), fill = TRUE,radius=0.1,fillOpacity = 0.1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~Brightness,title="CubeSat Brightness")

plot(NightRaster)
