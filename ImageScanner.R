##This script counts all of the saturated pixels in the raw image files.
##These counts are used to determine which image exposure produces the fewest saturated pixels for downstream analysis.
library(hexView)
library(broman)
library(plyr)
library(ggplot2)
library(jpeg)

wd <- "~/Desktop/Coastlight"
setwd(wd)

## Code for generating the saturation pixel count per image:
# If your images are in CR2 format please first convert them to jpg format.
# Unix command for this conversion:
# for d in ./*/ ; do (cd "$d" && for i in *.CR2; do sips -s format jpeg $i --out "${i%.*}.jpg"; done); done

directories <- list.dirs('.', recursive=FALSE)

imageSaturations <- data.frame()

for(directory in directories){
  fileNames <- list.files(path=directory,pattern="*.jpg")
  for(fileName in fileNames){
    imageFile <- paste(directory,fileName,sep="/")
    imageHistogram <- data.frame(Color = names(table(as.raster(readJPEG(imageFile)))), Count = as.integer(table(as.raster(readJPEG(imageFile)))))
    whitePixelCount <- imageHistogram[which(imageHistogram$Color=="#FFFFFF"),"Count"]
    row <- data.frame()
    row[1,1] <- gsub(".jpg",".CR2",fileName) # Replace the extension in the name, but not the actual filetype, for later analysis.
    row[1,2] <- whitePixelCount
    print(row)
    imageSaturations <- rbind(imageSaturations,row)
  }
}
colnames(imageSaturations) <- c("Filename","SaturatedPixelCount")
write.table(imageSaturations,"imageSaturations.txt",quote=FALSE,sep="\t",row.names = FALSE)
## End pixel saturation count section.

## Start here if pixel saturation counts are already saved to a file.
imageSaturations <- read.table("imageSaturations.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

coastLightData <- read.table("CoastlightData.tsv",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

CoastlightMerged <-  join(coastLightData,imageSaturations,by=c("Filename"))
