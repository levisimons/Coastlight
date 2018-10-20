##This script counts all of the saturated pixels in the raw image files.
##These counts are used to determine which image exposure produces the fewest saturated pixels for downstream analysis.
library(hexView)
library(broman)

setwd("~/Desktop/Coastlight")

directories <- list.dirs('.', recursive=FALSE)

imageSaturations <- data.frame()

for(directory in directories){
  fileNames <- list.files(path=directory)
  for(fileName in fileNames){
    imageFile <- paste(directory,fileName,sep="/")
    rawImage <- readRaw(imageFile) #Read .CR2 raw images as a rawBlock.
    saturationCount <- sum(hex2dec(rawImage$fileRaw)==255) #Count the number of saturated pixels (value=255)
    row <- data.frame()
    row[1,1] <- fileName
    row[1,2] <- saturationCount
    imageSaturations <- rbind(imageSaturations,row)
    print(nrow(imageSaturations))
  }
}
colnames(imageSaturations) <- c("Filename","SaturatedPixelCount")
write.table(imageSaturations,"imageSaturations.txt",quote=FALSE,sep="\t",row.names = FALSE)
