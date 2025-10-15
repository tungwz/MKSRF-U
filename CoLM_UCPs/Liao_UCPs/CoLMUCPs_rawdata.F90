PROGRAM main

  ! gfortran -mcmodel=large -g -fbounds-check -o mkucps CoLMUCPs_resample.F90 -I/usr/include -lnetcdf -lnetcdff

  USE netcdf

  IMPLICIT NONE

  INTEGER , parameter :: r8 = selected_real_kind(12)
  INTEGER , parameter :: slon = 43200, slat = 21600
  REAL(r8), parameter :: sdelta = 1._r8/120._r8

  REAL(r8), allocatable, dimension(:)   :: lat, lon
  REAL(r8), allocatable, dimension(:,:) :: h, haw, hsd, hl, fb, fgimp, fg, fg_

  INTEGER  :: ncid, varid, lat_id, lon_id, x, y, i, j, io, jo, ivar
  INTEGER  :: ht_id, wt_id, fgimp_id, fg_id, hsd_id,  haw_id, hl_id

  CHARACTER(len=255) :: ncfile
  CHARACTER(len=255) :: infile = './UCPs_tif'

  ALLOCATE( lat(slat) )
  ALLOCATE( lon(slon) )

  ALLOCATE( h    (slon,slat) )
  ALLOCATE( haw  (slon,slat) )
  ALLOCATE( hsd  (slon,slat) )
  ALLOCATE( hl   (slon,slat) )
  ALLOCATE( fb   (slon,slat) )
  ALLOCATE( fgimp(slon,slat) )
  ALLOCATE( fg   (slon,slat) )
  ALLOCATE( fg_  (slon,slat) )

  h    (:,:) = -999.
  haw  (:,:) = -999.
  hsd  (:,:) = -999.
  hl   (:,:) = -999.
  fb   (:,:) = -999.
  fgimp(:,:) = -999.
  fg   (:,:) = -999.
  fg_  (:,:) = -999.

  OPEN(11, file=trim(infile))

  ivar = 1
  DO WHILE(.true.)

     READ(11,*,end=100) ncfile
     print*, 'Processing '//trim(ncfile)//'.nc ...'

     CALL check( nf90_open(trim(ncfile)//'.nc', nf90_nowrite, ncid) )

     CALL check( nf90_inq_varid(ncid, "Band1", varid) )

     IF (ivar == 1) THEN
        CALL check( nf90_get_var  (ncid, varid  , fb   ) )
     ELSE IF(ivar == 2) THEN
        CALL check( nf90_get_var  (ncid, varid  , fg_  ) )
     ELSE IF(ivar == 3) THEN
        CALL check( nf90_get_var  (ncid, varid  , fgimp) )
     ELSE IF(ivar == 4) THEN
        CALL check( nf90_get_var  (ncid, varid  , h    ) )
     ELSE IF(ivar == 5) THEN
        CALL check( nf90_get_var  (ncid, varid  , hsd  ) )
     ELSE IF(ivar == 6) THEN
        CALL check( nf90_get_var  (ncid, varid  , haw  ) )
     ELSE
        CALL check( nf90_get_var  (ncid, varid  , hl   ) )
     ENDIF

     CALL check( nf90_close(ncid) )

     ivar = ivar + 1

  END DO

  100 close(11)

  DO i=1,slat
     lat(i) =  90. - i*sdelta + 0.5*sdelta
  ENDDO

  DO i=1,slon
     lon(i) = -180.+ i*sdelta - 0.5*sdelta
  ENDDO

  print*, "Ceating Global 1km data ..."

  CALL check( nf90_create('CoLM_UCPs_ROOF.nc', nf90_netcdf4, ncid) )

  CALL check( nf90_def_dim(ncid, "lat", slat, y) )
  CALL check( nf90_def_dim(ncid, "lon", slon, x) )

  CALL check( nf90_def_var(ncid, "lat", NF90_DOUBLE, y, lat_id, deflate_level=6) )
  CALL check( nf90_def_var(ncid, "lon", NF90_DOUBLE, x, lon_id, deflate_level=6) )

  CALL check( nf90_put_att(ncid, lat_id , "long_name", "Latitude"        ) )
  CALL check( nf90_put_att(ncid, lat_id , "units"    , "degrees_north"   ) )
  CALL check( nf90_put_att(ncid, lon_id , "long_name", "Longitude"       ) )
  CALL check( nf90_put_att(ncid, lon_id , "units"    , "degrees_east"    ) )

  CALL check( nf90_def_var(ncid, "HT_ROOF" , NF90_FLOAT, (/x,y/), ht_id   , deflate_level=6) )
  CALL check( nf90_def_var(ncid, "PCT_ROOF", NF90_FLOAT, (/x,y/), wt_id   , deflate_level=6) )

  CALL check( nf90_put_att(ncid, ht_id   , "long_name" , "height of roof"   ) )
  CALL check( nf90_put_att(ncid, ht_id   , "units"     , "meters"           ) )
  CALL check( nf90_put_att(ncid, ht_id   , "_FillValue", -999.              ) )

  CALL check( nf90_put_att(ncid, wt_id   , "long_name" , "fraction of roof" ) )
  CALL check( nf90_put_att(ncid, wt_id   , "units"     , "unitless"         ) )
  CALL check( nf90_put_att(ncid, wt_id   , "_FillValue", -999.              ) )

  CALL check( nf90_inq_varid(ncid, "lat" , lat_id) )
  CALL check( nf90_put_var  (ncid, lat_id, lat   ) )

  CALL check( nf90_inq_varid(ncid, "lon" , lon_id) )
  CALL check( nf90_put_var  (ncid, lon_id, lon   ) )

  CALL check( nf90_inq_varid(ncid, "HT_ROOF" , ht_id   ) )
  CALL check( nf90_put_var  (ncid, ht_id     , h       ) )

  CALL check( nf90_inq_varid(ncid, "PCT_ROOF", wt_id   ) )
  CALL check( nf90_put_var  (ncid, wt_id     , fb      ) )

  CALL check( nf90_close(ncid) )

  CALL check( nf90_create('CoLM_UCPs_HL.nc', nf90_netcdf4, ncid) )

  CALL check( nf90_def_dim(ncid, "lat", slat, y) )
  CALL check( nf90_def_dim(ncid, "lon", slon, x) )

  CALL check( nf90_def_var(ncid, "lat", NF90_DOUBLE, y, lat_id, deflate_level=6) )
  CALL check( nf90_def_var(ncid, "lon", NF90_DOUBLE, x, lon_id, deflate_level=6) )

  CALL check( nf90_put_att(ncid, lat_id , "long_name", "Latitude"        ) )
  CALL check( nf90_put_att(ncid, lat_id , "units"    , "degrees_north"   ) )
  CALL check( nf90_put_att(ncid, lon_id , "long_name", "Longitude"       ) )
  CALL check( nf90_put_att(ncid, lon_id , "units"    , "degrees_east"    ) )

  CALL check( nf90_def_var(ncid, "HL_BLD", NF90_FLOAT, (/x,y/), hl_id   , deflate_level=6) )

  CALL check( nf90_put_att(ncid, hl_id   , "long_name" , "Ratio of building height to average side length") )
  CALL check( nf90_put_att(ncid, hl_id   , "units"     , "unitless"                                       ) )
  CALL check( nf90_put_att(ncid, hl_id   , "_FillValue", -999.                                            ) )

  CALL check( nf90_inq_varid(ncid, "lat" , lat_id) )
  CALL check( nf90_put_var  (ncid, lat_id, lat   ) )

  CALL check( nf90_inq_varid(ncid, "lon" , lon_id) )
  CALL check( nf90_put_var  (ncid, lon_id, lon   ) )

  CALL check( nf90_inq_varid(ncid, "HL_BLD", hl_id   ) )
  CALL check( nf90_put_var  (ncid, hl_id   , hl      ) )

  CALL check( nf90_close(ncid) )

  DO i=1,slat
     DO j=1,slon
        IF (fb(j,i)<1. .and. fb(j,i)>0. .and. (fgimp(j,i)-fb(j,i))>0) THEN
           fg(j,i) = 1. - (fgimp(j,i)-fb(j,i))/(1-fb(j,i))

           IF (fg(j,i) > 1) THEN
              print*, 'fg>1: ', fb(j,i), fg_(j,i), fg(j,i)
              fg(j,i) = 1
           ENDIF
        ELSE
           fg(j,i) = -999.
        ENDIF
     ENDDO
  ENDDO

  CALL check( nf90_create('CoLM_UCPs_Fgper.nc', nf90_netcdf4, ncid) )

  CALL check( nf90_def_dim(ncid, "lat", slat, y) )
  CALL check( nf90_def_dim(ncid, "lon", slon, x) )

  CALL check( nf90_def_var(ncid, "lat", NF90_DOUBLE, y, lat_id, deflate_level=6) )
  CALL check( nf90_def_var(ncid, "lon", NF90_DOUBLE, x, lon_id, deflate_level=6) )

  CALL check( nf90_put_att(ncid, lat_id , "long_name", "Latitude"        ) )
  CALL check( nf90_put_att(ncid, lat_id , "units"    , "degrees_north"   ) )
  CALL check( nf90_put_att(ncid, lon_id , "long_name", "Longitude"       ) )
  CALL check( nf90_put_att(ncid, lon_id , "units"    , "degrees_east"    ) )

  CALL check( nf90_def_var(ncid, "WTROAD_PERV" , NF90_FLOAT, (/x,y/), fg_id   , deflate_level=6) )

  CALL check( nf90_put_att(ncid, fg_id   , "long_name" , "fraction of pervious road") )
  CALL check( nf90_put_att(ncid, fg_id   , "units"     , "unitless"                 ) )
  CALL check( nf90_put_att(ncid, fg_id   , "_FillValue", -999.                      ) )

  CALL check( nf90_inq_varid(ncid, "lat" , lat_id) )
  CALL check( nf90_put_var  (ncid, lat_id, lat   ) )

  CALL check( nf90_inq_varid(ncid, "lon" , lon_id) )
  CALL check( nf90_put_var  (ncid, lon_id, lon   ) )

  CALL check( nf90_inq_varid(ncid, "WTROAD_PERV" , fg_id   ) )
  CALL check( nf90_put_var  (ncid, fg_id         , fg      ) )

  CALL check( nf90_close(ncid) )

  print*, "make Global 500m data SUCCESSFULLY ..."

  deallocate ( lat   )
  deallocate ( lon   )
  deallocate ( h     )
  deallocate ( haw   )
  deallocate ( hsd   )
  deallocate ( hl    )
  deallocate ( fg    )
  deallocate ( fgimp )
  deallocate ( fb    )
CONTAINS

   SUBROUTINE check(status)
      INTEGER, intent(in) :: status

      IF (status /= nf90_noerr) THEN
         PRINT *, trim( nf90_strerror(status))
         stop 2
      ENDIF
   END SUBROUTINE check

END PROGRAM main


