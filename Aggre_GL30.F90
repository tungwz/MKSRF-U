Program Aggre

 use netcdf

 IMPLICIT NONE

 INTEGER, parameter :: r8   = selected_real_kind(12)
 INTEGER, parameter :: nlat = 43200 , nlon = 86400
 REAL(r8),parameter :: lon_min = -180.0, lon_max = 180.0
 REAL(r8),parameter :: lat_min =  -90.0, lat_max =  90.0
 REAL(r8),parameter :: delta = 1._r8/240._r8
 REAL(r8),parameter :: hdelta = 0.00025
 REAL    ,parameter :: spval = -999.

 REAL(r8)       , allocatable, dimension(:)   :: lat, lon, nclat, nclon
 INTEGER(kind=2), allocatable, dimension(:,:) :: band1
 REAL(r8)       , allocatable, dimension(:,:) :: cntdata
 REAL(r8)       , allocatable, dimension(:,:) :: watdata
 REAL(r8)       , allocatable, dimension(:,:) :: outdata

 REAL(r8) :: dlat, dlon, clat, clon
 INTEGER  :: tlat, tlon
 INTEGER  :: nyo , nxo
 INTEGER  :: i   , j
 INTEGER  :: ncid, water_id, lat_id, lon_id, nid, nclat_id, nclon_id, band_id
 INTEGER  :: xx  , yy
 INTEGER  :: argn
 
 CHARACTER(len=4)   :: year = '2020'
 CHARACTER(len=255) :: ncfile
 CHARACTER(len=255) :: infile = './ncgl' 
 CHARACTER(len=255) :: watfile= '/hard/dongwz/CoLM-U/waterbody/'
 
 ALLOCATE( nclat(nlat) )
 ALLOCATE( nclon(nlon) )

 ALLOCATE( watdata(nlon,nlat) )
 ALLOCATE( cntdata(nlon,nlat) )
 ALLOCATE( outdata(nlon,nlat) )

 cntdata(:,:) = 0
 watdata(:,:) = 0
 outdata(:,:) = -999._r8

 argn = IARGC()
 IF (argn > 0) THEN
    CALL getarg(1, year)
 ENDIF

 nclat(1) = delta/2 - 90.
 DO i=2,nlat,1
    nclat(i) = (i-1)*delta + nclat(1)
 ENDDO

 nclon(1) = delta/2 - 180.
 DO i=2,nlon,1
    nclon(i) = (i-1)*delta + nclon(1)
 ENDDO
    
 CALL check( nf90_create(trim(watfile)//'Water500m_'//trim(year)//'.nc', NF90_NETCDF4, ncid) )
 
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'created by', 'Wenzong Dong Wed Mar 1 CST 2022') )

 CALL check( nf90_def_dim(ncid, 'lat', nlat, xx) )
 CALL check( nf90_def_dim(ncid, 'lon', nlon, yy) )

 CALL check( nf90_def_var(ncid, 'lat' , NF90_FLOAT     , (/xx/), nclat_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, nclat_id, 'standard_name', 'latitude'     ) )
 CALL check( nf90_put_att(ncid, nclat_id, 'long_name'    , 'latitude'     ) )
 CALL check( nf90_put_att(ncid, nclat_id, 'units'        , 'degrees_north') )

 CALL check( nf90_def_var(ncid, 'lon' , NF90_FLOAT     , (/yy/), nclon_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, nclon_id, 'standard_name', 'longitude'    ) )
 CALL check( nf90_put_att(ncid, nclon_id, 'long_name'    , 'longitude'    ) )
 CALL check( nf90_put_att(ncid, nclon_id, 'units'        , 'degrees_east' ) )

 CALL check( nf90_def_var(ncid, 'PCT_Water',  NF90_FLOAT  , (/yy,xx/), water_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, water_id, 'long_name' , 'water percent' ) )
 CALL check( nf90_put_att(ncid, water_id, 'units'     , '%'             ) )
 CALL check( nf90_put_att(ncid, water_id, '_FillValue', spval           ) )

 CALL check(nf90_enddef(ncid))
 
 OPEN(11, file=trim(infile)//trim(year))

 DO WHILE(.true.)

    READ(11,*,end=100) ncfile                         
    PRINT*, ncfile

    CALL check( nf90_open('/hard/dongwz/GFCC/waterbody'//trim(year)//'/'//trim(ncfile), NF90_NOWRITE, nid  ) )

    CALL check( nf90_inq_dimid        (nid, 'lat' , lat_id ) )
    CALL check( nf90_inquire_dimension(nid, lat_id, len=nyo) )

    CALL check( nf90_inq_dimid        (nid, 'lon' , lon_id ) )
    CALL check( nf90_inquire_dimension(nid, lon_id, len=nxo) )

    !print*, nyo,nxo
    ALLOCATE( lat(nyo) )
    ALLOCATE( lon(nxo) )

    ALLOCATE(band1(nxo,nyo) )

    CALL check( nf90_inq_varid(nid, 'lat' , lat_id) )
    CALL check( nf90_get_var  (nid, lat_id, lat   ) )

    CALL check( nf90_inq_varid(nid, 'lon' , lon_id) )
    CALL check( nf90_get_var  (nid, lon_id, lon   ) )

    CALL check( nf90_inq_varid(nid, 'Band1' , band_id) )
    CALL check( nf90_get_var  (nid, band_id, band1   ) )
 
    CALL check( nf90_close(nid) )

    DO i=1,nxo,1
       DO j=1,nyo,1
         
          clat = lat(j) + hdelta/2
          clon = lon(i) + hdelta/2
          IF (clon>=180.) THEN
             clon = lon(i)
          ENDIF

          ! 计算在哪个500m格子
          tlat = NINT((clat+ 90.)/delta+0.5)
          tlon = NINT((clon+180.)/delta+0.5)


          ! 只有当band=60时累加(水体)
          IF (band1(i,j).ne.0 .and. band1(i,j).ne.255) THEN
                IF (band1(i,j) .eq. 60) THEN
                    watdata(tlon,tlat) = watdata(tlon,tlat) + 1.
                ENDIF
                
                cntdata(tlon,tlat) = cntdata(tlon,tlat) + 1.
          ENDIF
       ENDDO
    ENDDO

    DEALLOCATE(lat  )
    DEALLOCATE(lon  )
    DEALLOCATE(band1)
 ENDDO

 100 close(11)

 DO i=1,nlon,1
    DO j=1,nlat,1  
       IF (cntdata(i,j) .ne. 0) THEN
          outdata(i,j)  = watdata(i,j)*100. / cntdata(i,j)
       ENDIF
    ENDDO
 ENDDO 

 CALL check( nf90_inq_varid(ncid, 'lat' , nclat_id) )
 CALL check( nf90_put_var  (ncid, nclat_id, nclat ) )

 CALL check( nf90_inq_varid(ncid, 'lon' , nclon_id) )
 CALL check( nf90_put_var  (ncid, nclon_id, nclon ) )

 CALL check( nf90_inq_varid(ncid, 'PCT_Water', water_id) )
 CALL check( nf90_put_var  (ncid, water_id   , outdata ) )

 CALL check( nf90_close(ncid) )

 DEALLOCATE(nclat  )
 DEALLOCATE(nclon  )
 DEALLOCATE(watdata)
 DEALLOCATE(cntdata)
 DEALLOCATE(outdata)

CONTAINS

 SUBROUTINE check(status)
    INTEGER, intent ( in) :: status

    IF(status /= nf90_noerr) THEN
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    ENDIF
 END SUBROUTINE check

END PROGRAM Aggre
