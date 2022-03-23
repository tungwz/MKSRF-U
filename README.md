# Urban data for model
gfortran -o mkglobal -fbounds-check -g MakeGlobalUrban.F90 -I$Inc -L$Lib -lnetcdf -lnetcdff
gfortran -o mkreg -fbounds-check  -g MakeRegionUrban.F90 -I$Inc -L$Lib -lnetcdf -lnetcdff
# Aggregate 30m data to 500m
gfortran -o mkgfcc -fbounds-check -fopenmp -g Aggre_GFCC.F90 -I$Inc -L$Lib -lnetcdf -lnetcdff


NOTE:
If errors occur when run mkglobal or mkreg (i.e., 0 of dimenssion....), please first run nearest.py to fill 0 grid which
is not aggre with MODIS data and overwriting old data.
