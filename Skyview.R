#Script to generate a sky view factor map raster using a digital elevation map raster.
require(rgdal)
require(horizon)
require(raster)

setwd("/home/cmb-07/sn1/alsimons/Coastlight") #Working directory on the cluster

#Read in digital elevation map and convert to a raster object.
DEM <- raster("DEM2mRaster.tif")

svf(DEM,nAngles = 16,maxDist = 1000,filename = "Skyview2mRaster.tif", blockSize = 512)
