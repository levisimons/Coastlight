# This script will run gdalwarp on each LIDAR tif file from NOAA and create a cleaned, then compressed, version with erroneous values removed.
for f in *.tif
do
        inputFile="${f}"
        outputFile="${f/%.tif/_cleaned.tif}"
        echo "${inputFile} ${outputFile}"
        gdalwarp -of ENVI -dstnodata -999999 $inputFile $outputFile
        compressedFile="${f/%.tif/_compressed.tif}"
        gdal_translate -co compress=lzw $outputFile $compressedFile
done
