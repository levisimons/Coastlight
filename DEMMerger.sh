#!/bin/sh
#SBATCH -t300:00:00
#SBATCH -n1
#SBATCH -c1
#SBATCH --mem=1200000m
#SBATCH -pSCCWRP
cd /home/cmb-07/sn1/alsimons/Coastlight
gdalbuildvrt mergedCoast.vrt /home/cmb-07/sn1/alsimons/Coastlight/*_compressed.tif
gdal_translate -of GTiff -co compress=lzw mergedCoast.vrt mergedCoast.tif 
