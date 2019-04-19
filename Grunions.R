require(dplyr)
##This script reads in monthly grunion spawning data and finds the spawning event for a given
##site which has the highest Walker scale value for the year.

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
  #Get site name
  siteName <- unlist(strsplit(grunionFile,split='_',fixed=TRUE))[1]
  #Rename maximum spawn column.
  colnames(grunionMax)[which(names(grunionMax) == "Walker")] <- paste("MaxWalker",siteName,sep="")
  #Combine all of the max spawn events by site and year into a single data frame.
  if(ncol(allSitesGrunionMax) == 0){
    allSitesGrunionMax <- grunionMax
  } else{
    allSitesGrunionMax <- cbind(allSitesGrunionMax,grunionMax)
  }
}
#Remove duplicated year column from max spawn events by site and year data frame
allSitesGrunionMax <- allSitesGrunionMax[,!duplicated(colnames(allSitesGrunionMax))]
