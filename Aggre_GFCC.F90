Program Aggre

 USE omp_lib
 USE netcdf

 IMPLICIT NONE

 INTEGER , parameter :: r8   = selected_real_kind(12)
 INTEGER , parameter :: nlat = 43200 , nlon = 86400
 REAL(r8), parameter :: lon_min = -180.0, lon_max = 180.0
 REAL(r8), parameter :: lat_min =  -90.0, lat_max =  90.0
 REAL(r8), parameter :: hdelta = 0.00025
 REAL(r8), parameter :: delta = 1._r8/240._r8
 REAL    , parameter :: spval = -999.

 REAL(r8), allocatable, dimension(:)   :: lat, lon, nclat, nclon
 INTEGER , allocatable, dimension(:,:) :: cntdata
 INTEGER , allocatable, dimension(:,:) :: sumdata
 REAL(r8), allocatable, dimension(:,:) :: outdata

 INTEGER(kind=2), allocatable, dimension(:,:) :: band1

 REAL(r8) :: dlat, dlon, clat, clon
 INTEGER  :: tlat, tlon
 INTEGER  :: nyo , nxo
 INTEGER  :: xid , yid
 INTEGER  :: i   , j
 INTEGER  :: ncid, tree_id, lat_id, lon_id, nid, nclat_id, nclon_id, band_id
 INTEGER  :: xx  , yy
 INTEGER  :: argn
 
 CHARACTER(len=255) :: ncfile
 CHARACTER(len=255) :: infile = './ncinfos' 
 CHARACTER(len=255) :: outfile= '/hard/dongwz/CoLM-U/GFCC500m/'
 CHARACTER(len=4  ) :: year = "2000"

 argn = IARGC()
 IF (argn > 0) THEN
    CALL getarg(1, year)
 ENDIF

 ALLOCATE( nclat(nlat) )
 ALLOCATE( nclon(nlon) )

 ALLOCATE( sumdata(nlon,nlat) )
 ALLOCATE( cntdata(nlon,nlat) )
 ALLOCATE( outdata(nlon,nlat) )

 sumdata(:,:) = 0
 cntdata(:,:) = 0
 outdata(:,:) = -999.

 nclat(1) = delta/2 - 90
 DO i=2,nlat,1
    nclat(i) = (i-1)*delta + nclat(1)
 ENDDO

 nclon(1) = delta/2 - 180
 DO i=2,nlon,1
    nclon(i) = (i-1)*delta + nclon(1)
 ENDDO

! CALL check( nf90_create(trim(outfile)//'GFCC500m_test'//trim(year)//'.nc', NF90_NETCDF4, ncid) )
 
! CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'created by', 'Wenzong Dong Tuesday Mar 1 CST 2022') )

! CALL check( nf90_def_dim(ncid, 'lat', nlat, xx) )
! CALL check( nf90_def_dim(ncid, 'lon', nlon, yy) )

! CALL check( nf90_def_var(ncid, 'lat' , NF90_FLOAT     , (/xx/), nclat_id, deflate_level=6) )
! CALL check( nf90_put_att(ncid, nclat_id, 'standard_name', 'latitude'     ) )
! CALL check( nf90_put_att(ncid, nclat_id, 'long_name'    , 'latitude'     ) )
! CALL check( nf90_put_att(ncid, nclat_id, 'units'        , 'degrees_north') )

! CALL check( nf90_def_var(ncid, 'lon' , NF90_FLOAT     , (/yy/), nclon_id, deflate_level=6) )
! CALL check( nf90_put_att(ncid, nclon_id, 'standard_name', 'longitude'   ) )
! CALL check( nf90_put_att(ncid, nclon_id, 'long_name'    , 'longitude'   ) )
! CALL check( nf90_put_att(ncid, nclon_id, 'units'        , 'degrees_east') )

! CALL check( nf90_def_var(ncid, 'PCT_Tree', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
! CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , 'forest cover'    ) )
! CALL check( nf90_put_att(ncid, tree_id   , 'units'     , '%'               ) )
! CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

! CALL check( nf90_enddef(ncid) )

 OPEN(11, file=trim(infile)//trim(year))

 DO WHILE(.true.)

    READ(11,*,end=100) ncfile                     
    PRINT*, ncfile

    CALL check( nf90_open('/hard/dongwz/GFCC/'//trim(year)//'/'//trim(ncfile), NF90_NOWRITE, nid  ) )

    CALL check( nf90_inq_dimid        (nid, 'lat' , lat_id ) )
    CALL check( nf90_inquire_dimension(nid, lat_id, len=nyo) )

    CALL check( nf90_inq_dimid        (nid, 'lon' , lon_id ) )
    CALL check( nf90_inquire_dimension(nid, lon_id, len=nxo) )

    print*, nyo,nxo
    ALLOCATE( lat(nyo) )
    ALLOCATE( lon(nxo) )

    ALLOCATE( band1(nxo,nyo) )

    CALL check( nf90_inq_varid(nid, 'lat' , lat_id) )
    CALL check( nf90_get_var  (nid, lat_id, lat   ) )

    CALL check( nf90_inq_varid(nid, 'lon' , lon_id) )
    CALL check( nf90_get_var  (nid, lon_id, lon   ) )

    CALL check( nf90_inq_varid(nid, 'Band1', band_id) )
    CALL check( nf90_get_var  (nid, band_id, band1  ) )
 
    CALL check( nf90_close(nid) )

    !print*, "strat processing"
    !write(*,*) "线程数为:",omp_get_num_procs()
!$OMP PARALLEL DO NUM_THREADS(92) &
!$OMP PRIVATE(i,j,clat,clon,tlat,tlon)
    DO i=1,nxo
       DO j=1,nyo
          
          clat = lat(j) + hdelta/2.
          clon = lon(i) + hdelta/2.

          ! 计算在哪个500m格子上
          tlat = NINT((clat+ 90.)/delta+0.5)
          tlon = NINT((clon+180.)/delta+0.5)

          ! 当band>=0且band<=100时累加，并计算有效值数量(cntdata)，用与求平均
          IF (band1(i,j)>=0 .and. band1(i,j)<=100) THEN
             sumdata(tlon,tlat) = sumdata(tlon,tlat) + band1(i,j)
             cntdata(tlon,tlat) = cntdata(tlon,tlat) + 1
          ENDIF
       ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !print*, "end"

    DEALLOCATE( lat   )
    DEALLOCATE( lon   )
    DEALLOCATE( band1 )
 ENDDO

 100 close(11)

 DO i=1,nlon,1
    DO j=1,nlat,1     
       IF (cntdata(i,j) .ne. 0) THEN
          outdata(i,j) = sumdata(i,j)*1. / cntdata(i,j)
       ENDIF
    ENDDO
 ENDDO 

 CALL check( nf90_create(trim(outfile)//'GFCC500mtest'//trim(year)//'.nc', NF90_NETCDF4, ncid) )

 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'created by', 'Wenzong Dong Tuesday Mar 1 CST 2022') )

 CALL check( nf90_def_dim(ncid, 'lat', nlat, xx) )
 CALL check( nf90_def_dim(ncid, 'lon', nlon, yy) )

 CALL check( nf90_def_var(ncid, 'lat' , NF90_FLOAT     , (/xx/), nclat_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, nclat_id, 'standard_name', 'latitude'     ) )
 CALL check( nf90_put_att(ncid, nclat_id, 'long_name'    , 'latitude'     ) )
 CALL check( nf90_put_att(ncid, nclat_id, 'units'        , 'degrees_north') )

 CALL check( nf90_def_var(ncid, 'lon' , NF90_FLOAT     , (/yy/), nclon_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, nclon_id, 'standard_name', 'longitude'   ) )
 CALL check( nf90_put_att(ncid, nclon_id, 'long_name'    , 'longitude'   ) )
 CALL check( nf90_put_att(ncid, nclon_id, 'units'        , 'degrees_east') )

 CALL check( nf90_def_var(ncid, 'PCT_Tree', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , 'forest cover'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , '%'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_enddef(ncid) )

 CALL check( nf90_inq_varid(ncid, 'lat' , nclat_id) )
 CALL check( nf90_put_var  (ncid, nclat_id, nclat ) )

 CALL check( nf90_inq_varid(ncid, 'lon' , nclon_id) )
 CALL check( nf90_put_var  (ncid, nclon_id, nclon ) )

 CALL check( nf90_inq_varid(ncid, 'PCT_Tree', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , outdata) )

 CALL check( nf90_close(ncid) )

 DEALLOCATE(nclat  )
 DEALLOCATE(nclon  )
 DEALLOCATE(sumdata)
 DEALLOCATE(cntdata)
 DEALLOCATE(outdata)

 contains

 SUBROUTINE check(status)
    INTEGER, intent (in) :: status

    IF(status /= nf90_noerr) THEN
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    END IF
 END SUBROUTINE check

END PROGRAM Aggre
