#!/bin/bash

# Define resolution argument
reso=$1

# Define the array of UCPs
# ucps=('HL' 'Fb' 'Fgimp' 'H' 'Haw' 'Hsd' 'Fg')
ucps=('Fb' 'Fg' 'Fgimp' 'H' 'Hsd' 'Haw' 'Hl')

# Read lines from the file into an array
mapfile -t lines < UCPs_tif

# convert tif to NetCDF
echo "convert tif to NetCDF"

i=0
for line in "${lines[@]}"; do
    echo "Processing ${ucps[i]}..."

    gdalwarp "/tera12/yuanhua/dongwz/github/CoLMUCPs/Colm_250224_SetNoData/${line}.tif" "${line}.tif" -dstnodata -999 -ts 43200 21600 -te -180 -90 180 90 -co COMPRESS=LZW -t_srs "+proj=longlat +ellps=WGS84"
    gdal_translate -a_nodata -999 -of netCDF  -co FORMAT=NC4 -co COMPRESS=DEFLATE -co ZLEVEL=6 "${line}.tif" "${line}_.nc"

    # invert lat from -90~90 to 90~-90
    cdo -z zip_6 invertlat "${line}_.nc" "${line}.nc"
    ((i++))
done

echo "convert done"
rm -rf global_*.tif
rm -rf global_*_.nc

# resample use nearest interpolation
echo "Resample to 500m use nearest interpolation method"
# gfortran -mcmodel=large -g -fbounds-check -o mkucps CoLMUCPs_resample.F90 -I/usr/include -lnetcdf -lnetcdff
# ./mkucps

# aggregation for a model resolution
# gfortran -mcmodel=large -g -fbounds-check -o aggre_ucps CoLMUCPs_aggregation.F90 -I/usr/include -lnetcdf -lnetcdff
# ./aggre_ucps "$reso"

# split for 5x5 region
python UCPs_5x5.py
