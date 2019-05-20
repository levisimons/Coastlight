#!/bin/sh
#SBATCH -t300:00:00
#SBATCH -n1
#SBATCH -c1
#SBATCH --mem=1200000m
#SBATCH -pSCCWRP
cd /home/cmb-07/sn1/alsimons/Coastlight
gdalbuildvrt mergedCoast.vrt /home/cmb-07/sn1/alsimons/Coastlight/Job470138_west_coast_2016_el_nino_dem_m6260_*.tif
gdal_translate -of GTiff -co compress=lzw mergedCoast.vrt mergedCoast.tif 
