PROGRAM clmu2grid

   USE netcdf

   IMPLICIT NONE

   INTEGER , PARAMETER :: r8 = SELECTED_REAL_KIND(12)

   REAL(r8), PARAMETER :: sdelta = 1._r8/240._r8
   REAL(r8), PARAMETER :: udelta = 1._r8/120._r8
   REAL    , PARAMETER :: fv     = -999

   INTEGER , PARAMETER :: nlat = 43200, nlon = 86400
   INTEGER , PARAMETER :: nxy  = 1200
   INTEGER , PARAMETER :: rid  = 33  , den_clss = 3
   INTEGER , PARAMETER :: nxo  = 720 , nyo = 360
   INTEGER , PARAMETER :: ns   = 2   , nr  = 2
   INTEGER , PARAMETER :: ulev = 10
   INTEGER , PARAMETER :: mon  = 12
   INTEGER , PARAMETER :: npft = 16 
   
! input variables
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlats
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlatn
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlon
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlonw
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlone
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlats
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlatn
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlon
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlonw
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlone
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gfcc_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gedi_th
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: modur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gl30_wt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: harea

   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urrgid
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urden

   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: pwater, lcdata, pglacier
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: pcrop, pwetland
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: pct_urban
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: pct_pft, ppft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: htop_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: lai_pft, pftlai
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: sai_pft, pftsai

   REAL(r8), DIMENSION(den_clss,rid)      :: hwrcan, wtrf, wtrd, emrf, emwl
   REAL(r8), DIMENSION(den_clss,rid)      :: emimrd, emperd, htrf, whc, ulevimrd 
   REAL(r8), DIMENSION(den_clss,rid)      :: thrf, thwl, tbmin, tbmax
   
   REAL(r8), DIMENSION(den_clss,rid,ulev ):: cvrf, cvwl, cvimrd, &
                                                tkrf, tkwl, tkimrd
   REAL(r8), DIMENSION(den_clss,rid,nr,ns):: albrf, albwl, albimrd, albperd
   ! output variables
   INTEGER , ALLOCATABLE, DIMENSION(:,:)   :: ur_clss
   INTEGER , ALLOCATABLE, DIMENSION(:,:)   :: ur_rgid
   
   REAL(r8), ALLOCATABLE, DIMENSION(:)     :: latso, flatso
   REAL(r8), ALLOCATABLE, DIMENSION(:)     :: lonso
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: area, pct_land
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wgt_top
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: pct_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: pct_urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: hwr_can
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wt_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wt_rd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_perd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: ht_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: w_hc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: ulev_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: th_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: th_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tb_min
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tb_max

   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: hgt, avg
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: ur_dc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: pct_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ur_lai
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ur_sai
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_imrd

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_perd

   ! variable ids
   INTEGER :: reg(4)
   INTEGER :: ncid, uhd_vid, umd_vid, utbd_vid, htop_urvid
   INTEGER :: lat_vid, lon_vid, lat_dimid, lon_dimid, &
              ns_dimid, nr_dimid, ulev_dimid
   INTEGER :: ns_vid, nr_vid, pct_tcvid, pct_urvid, pct_urwtvid, ulev_vid, lev_vid
   INTEGER :: pftvid, laivid, maskid, umaskid, t_pftvid, ur_laivid, htop_pftvid
   INTEGER :: lcvid, pcropvid, pwetlandvid, pwatervid, pglaciervid, ppftvid
   INTEGER :: hlat_vid, hlon_vid, gfcc_tcvid, gedi_thvid, gl30_wtvid
   INTEGER :: urlat_vid, urlon_vid, ur_clssvid, ur_rgvid, hwr_canvid
   INTEGER :: wt_rfvid, wt_rdvid, em_rfvid, em_wlvid, em_imrdvid, em_perdvid
   INTEGER :: ht_rfvid, whcvid, cv_rfvid, cv_wlvid, cv_imrdvid, ulev_imrdvid
   INTEGER :: th_rfvid, th_wlvid, tbminvid, tbmaxvid
   INTEGER :: tk_rfvid, tk_wlvid, tk_imrdvid
   INTEGER :: alb_rfvid, alb_imrdvid, alb_perdvid, alb_wlvid
   INTEGER :: uxid, uyid, upftvid, mon_dimid, den_dimid
   INTEGER :: den_vid, mon_vid, ur_landvid, saivid, ur_saivid, ur_denvid

   CHARACTER(len=255) :: SRF_DIR='/hard/dongwz/CoLM-U/urban_5x5/'
   CHARACTER(len=255) :: RAW_DIR='/tera02/yuanhua/mksrf/srf_5x5/'
   CHARACTER(len=255) :: OUT_DIR='./'
   CHARACTER(len=255) :: REGFILE='reg_5x5'
   CHARACTER(len=255) :: filename
   CHARACTER(len=4)    , DIMENSION(4) :: creg
   CHARACTER(len=4)   :: year="2005"

   REAL(r8) :: fac(3)
   REAL(r8) :: pi, deg2rad, re, dx, dy, sumarea, sumur, sumpct
   REAL(r8) :: dll, x_delta, y_delta, wgt
   REAL(r8) :: lone1(nxo), lonw1(nxo), latn1(nyo), lats1(nyo)

   INTEGER(kind=2):: lc
   INTEGER  :: i, j, k, io, jo, m, n, ii, jj, cont, sumth, p, ip
   INTEGER  :: argn
   INTEGER  :: n_ns(2), n_nr(2), n_den(3), n_ulev(10), n_mon(12)
   INTEGER  :: XY2D(2), XY3D(3), XY4D(4), UR3D(3), UL3D(4), XY5D(5)

   LOGICAL :: fileExists, fileNotExists

   argn = IARGC()
   IF (argn > 0) THEN
      CALL getarg(1, year)
   ENDIF

   pi = 4._r8*atan(1.)
   deg2rad = pi/180._r8
   re = 6.37122e6 * 0.001

   ALLOCATE( hlat  (nxy) )
   ALLOCATE( hlats (nxy) )
   ALLOCATE( hlatn (nxy) )
   ALLOCATE( hlon  (nxy) )
   ALLOCATE( hlonw (nxy) )
   ALLOCATE( hlone (nxy) )

   ALLOCATE( harea   (nxy, nxy) )

   ALLOCATE( gfcc_tc (nxy, nxy) )
   ALLOCATE( gedi_th (nxy, nxy) )
   ALLOCATE( gl30_wt (nxy, nxy) )
   ALLOCATE( modur   (nxy, nxy) )
   ALLOCATE( urrgid  (nxy, nxy) )
   ALLOCATE( urden   (nxy, nxy) )

   ALLOCATE( latso     (nyo) )
   ALLOCATE( lonso     (nxo) )
   ALLOCATE( pct_land  (nxo, nyo) )
   ALLOCATE( ur_clss   (nxo, nyo) )
   ALLOCATE( ur_rgid   (nxo, nyo) )
   ALLOCATE( area      (nxo, nyo) )
   ALLOCATE( hgt       (nxo, nyo) )
   ALLOCATE( avg       (nxo, nyo) )
   ALLOCATE( pct_urban (nxo, nyo) )
   ALLOCATE( htop_pft  (nxo, nyo, npft     ) )
   ALLOCATE( pct_pft   (nxo, nyo, npft     ) )
   ALLOCATE( lai_pft   (nxo, nyo, npft, mon) )
   ALLOCATE( sai_pft   (nxo, nyo, npft, mon) )
   ALLOCATE( ur_dc     (nxo, nyo, den_clss ) )
   ALLOCATE( pct_ur    (nxo, nyo, den_clss ) )
   ALLOCATE( wgt_top   (nxo, nyo, den_clss ) )
   ALLOCATE( tc        (nxo, nyo, den_clss ) )
   ALLOCATE( urwt      (nxo, nyo, den_clss ) )
   ALLOCATE( htop      (nxo, nyo, den_clss ) )
   ALLOCATE( pct_tc    (nxo, nyo, den_clss ) )
   ALLOCATE( pct_urwt  (nxo, nyo, den_clss ) )
   ALLOCATE( htop_ur   (nxo, nyo, den_clss ) )
   ALLOCATE( hwr_can   (nxo, nyo, den_clss ) )
   ALLOCATE( wt_rf     (nxo, nyo, den_clss ) )
   ALLOCATE( wt_rd     (nxo, nyo, den_clss ) )
   ALLOCATE( em_rf     (nxo, nyo, den_clss ) )
   ALLOCATE( em_wl     (nxo, nyo, den_clss ) )
   ALLOCATE( em_imrd   (nxo, nyo, den_clss ) )
   ALLOCATE( em_perd   (nxo, nyo, den_clss ) )
   ALLOCATE( ht_rf     (nxo, nyo, den_clss ) )
   ALLOCATE( w_hc      (nxo, nyo, den_clss ) )
   ALLOCATE( ulev_imrd (nxo, nyo, den_clss ) )
   ALLOCATE( th_rf     (nxo, nyo, den_clss ) )
   ALLOCATE( th_wl     (nxo, nyo, den_clss ) )
   ALLOCATE( tb_min    (nxo, nyo, den_clss ) )
   ALLOCATE( tb_max    (nxo, nyo, den_clss ) )
   ALLOCATE( ur_sai    (nxo, nyo, den_clss, mon ) )
   ALLOCATE( ur_lai    (nxo, nyo, den_clss, mon ) )
   ALLOCATE( cv_rf     (nxo, nyo, den_clss, ulev) )
   ALLOCATE( cv_wl     (nxo, nyo, den_clss, ulev) )
   ALLOCATE( cv_imrd   (nxo, nyo, den_clss, ulev) )
   ALLOCATE( tk_rf     (nxo, nyo, den_clss, ulev) )
   ALLOCATE( tk_wl     (nxo, nyo, den_clss, ulev) )
   ALLOCATE( tk_imrd   (nxo, nyo, den_clss, ulev) )
   ALLOCATE( alb_rf    (nxo, nyo, den_clss, nr, ns) )
   ALLOCATE( alb_wl    (nxo, nyo, den_clss, nr, ns) )
   ALLOCATE( alb_imrd  (nxo, nyo, den_clss, nr, ns) )
   ALLOCATE( alb_perd  (nxo, nyo, den_clss, nr, ns) )

   ! initialization
   fac      (:)     = 0.
   hgt      (:,:)   = 0.
   avg      (:,:)   = 0.
   tc       (:,:,:) = 0.
   urwt     (:,:,:) = 0.
   htop     (:,:,:) = 0.
   pct_tc   (:,:,:) = 0.
   pct_urwt (:,:,:) = 0.
   htop_ur  (:,:,:) = 0.
   ur_dc    (:,:,:) = 0.

   hwr_can  (:,:,:) = 0.
   wt_rf    (:,:,:) = 0.
   wt_rd    (:,:,:) = 0.
   em_rf    (:,:,:) = 0.
   em_wl    (:,:,:) = 0.
   em_imrd  (:,:,:) = 0.
   em_perd  (:,:,:) = 0.
   ht_rf    (:,:,:) = 0.
   w_hc     (:,:,:) = 0.
   ulev_imrd(:,:,:) = 0.
   th_rf    (:,:,:) = 0.
   th_wl    (:,:,:) = 0.
   tb_min   (:,:,:) = 0.
   tb_max   (:,:,:) = 0.

   ur_lai   (:,:,:,:) = 0.
   ur_sai   (:,:,:,:) = 0.
   tk_rf    (:,:,:,:) = 0.
   tk_wl    (:,:,:,:) = 0.
   tk_imrd  (:,:,:,:) = 0.
   cv_rf    (:,:,:,:) = 0.
   cv_wl    (:,:,:,:) = 0.
   cv_imrd  (:,:,:,:) = 0.

   alb_rf   (:,:,:,:,:) = 0.
   alb_wl   (:,:,:,:,:) = 0.
   alb_imrd (:,:,:,:,:) = 0.
   alb_perd (:,:,:,:,:) = 0.

   ur_rgid (:,:) = 0.
   ur_clss (:,:) = 0.  
   hwrcan  (:,:) = 0.
   wtrf    (:,:) = 0.
   wtrd    (:,:) = 0.
   emrf    (:,:) = 0.
   emwl    (:,:) = 0.
   emimrd  (:,:) = 0.
   emperd  (:,:) = 0.
   htrf    (:,:) = 0.
   whc     (:,:) = 0.
   ulevimrd(:,:) = 0.
   thrf    (:,:) = 0.
   thwl    (:,:) = 0.
   tbmin   (:,:) = 0.
   tbmax   (:,:) = 0.

   tkrf    (:,:,:) = 0.
   tkwl    (:,:,:) = 0.
   tkimrd  (:,:,:) = 0.
   cvrf    (:,:,:) = 0.
   cvwl    (:,:,:) = 0.
   cvimrd  (:,:,:) = 0.

   albrf   (:,:,:,:) = 0.
   albwl   (:,:,:,:) = 0.
   albimrd (:,:,:,:) = 0.

   cont = 0.
   sumth= 0.

   fileNotExists = .FALSE.

   x_delta = (360._r8/nxo)*1._r8
   y_delta = (180._r8/nyo)*1._r8

   ! get urban properties data from NCAR 1km urban data
   CALL check( nf90_open('urban_properties_data.1km.210726-121520.nc', nf90_nowrite, ncid) )

   CALL check( nf90_inq_varid(ncid, "CANYON_HWR"         , hwr_canvid  ) )
   CALL check( nf90_inq_varid(ncid, "WTLUNIT_ROOF"       , wt_rfvid    ) )
   CALL check( nf90_inq_varid(ncid, "WTROAD_PERV"        , wt_rdvid    ) )
   CALL check( nf90_inq_varid(ncid, "EM_ROOF"            , em_rfvid    ) )
   CALL check( nf90_inq_varid(ncid, "EM_WALL"            , em_wlvid    ) )
   CALL check( nf90_inq_varid(ncid, "EM_IMPROAD"         , em_imrdvid  ) )
   CALL check( nf90_inq_varid(ncid, "EM_PERROAD"         , em_perdvid  ) )
   CALL check( nf90_inq_varid(ncid, "ALB_ROOF"           , alb_rfvid   ) )
   CALL check( nf90_inq_varid(ncid, "ALB_WALL"           , alb_wlvid   ) )
   CALL check( nf90_inq_varid(ncid, "ALB_IMPROAD"        , alb_imrdvid ) )
   CALL check( nf90_inq_varid(ncid, "ALB_PERROAD"        , alb_perdvid ) )
   CALL check( nf90_inq_varid(ncid, "HT_ROOF"            , ht_rfvid    ) )
   CALL check( nf90_inq_varid(ncid, "WIND_HGT_CANYON"    , whcvid      ) )
   CALL check( nf90_inq_varid(ncid, "TK_ROOF"            , tk_rfvid    ) )
   CALL check( nf90_inq_varid(ncid, "TK_WALL"            , tk_wlvid    ) )
   CALL check( nf90_inq_varid(ncid, "TK_IMPROAD"         , tk_imrdvid  ) )
   CALL check( nf90_inq_varid(ncid, "CV_ROOF"            , cv_rfvid    ) )
   CALL check( nf90_inq_varid(ncid, "CV_WALL"            , cv_wlvid    ) )
   CALL check( nf90_inq_varid(ncid, "CV_IMPROAD"         , cv_imrdvid  ) )
   CALL check( nf90_inq_varid(ncid, "NLEV_IMPROAD"       , ulev_imrdvid) )
   CALL check( nf90_inq_varid(ncid, "THICK_ROOF"         , th_rfvid    ) )
   CALL check( nf90_inq_varid(ncid, "THICK_WALL"         , th_wlvid    ) )
   CALL check( nf90_inq_varid(ncid, "T_BUILDING_MIN"     , tbminvid    ) )
   CALL check( nf90_inq_varid(ncid, "T_BUILDING_MAX"     , tbmaxvid    ) )

   CALL check( nf90_get_var(ncid, hwr_canvid  , hwrcan  ) )
   CALL check( nf90_get_var(ncid, wt_rfvid    , wtrf    ) )
   CALL check( nf90_get_var(ncid, wt_rdvid    , wtrd    ) )
   CALL check( nf90_get_var(ncid, em_rfvid    , emrf    ) )
   CALL check( nf90_get_var(ncid, em_wlvid    , emwl    ) )
   CALL check( nf90_get_var(ncid, em_imrdvid  , emimrd  ) )
   CALL check( nf90_get_var(ncid, em_perdvid  , emperd  ) )
   CALL check( nf90_get_var(ncid, alb_rfvid   , albrf   ) )
   CALL check( nf90_get_var(ncid, alb_wlvid   , albwl   ) )
   CALL check( nf90_get_var(ncid, alb_imrdvid , albimrd ) )
   CALL check( nf90_get_var(ncid, alb_perdvid , albperd ) )
   CALL check( nf90_get_var(ncid, ht_rfvid    , htrf    ) )
   CALL check( nf90_get_var(ncid, whcvid      , whc     ) )
   CALL check( nf90_get_var(ncid, tk_rfvid    , tkrf    ) )
   CALL check( nf90_get_var(ncid, tk_wlvid    , tkwl    ) )
   CALL check( nf90_get_var(ncid, tk_imrdvid  , tkimrd  ) )
   CALL check( nf90_get_var(ncid, cv_rfvid    , cvrf    ) )
   CALL check( nf90_get_var(ncid, cv_wlvid    , cvwl    ) )
   CALL check( nf90_get_var(ncid, cv_imrdvid  , cvimrd  ) )
   CALL check( nf90_get_var(ncid, ulev_imrdvid, ulevimrd) )
   CALL check( nf90_get_var(ncid, th_rfvid    , thrf    ) )
   CALL check( nf90_get_var(ncid, th_wlvid    , thwl    ) )
   CALL check( nf90_get_var(ncid, tbminvid    , tbmin   ) )
   CALL check( nf90_get_var(ncid, tbmaxvid    , tbmax   ) )

   CALL check( nf90_close(ncid) )

   PRINT *, "    reading MODIS PFTs data"

   ! read model grid PFT/lai/sai/htop data for urban lai calculating
   filename = "global_0.5x0.5.MOD"//TRIM(year)//"_v5.nc"
   inquire (file=filename, exist=fileExists)
   IF (fileExists) THEN
      ! 如果对应分辨率的全球PFTs/LAI/SAI/HTOP_PFT的文件已经生成
      ! 则从该文件中直接读取
      CALL check( nf90_open(filename, nf90_nowrite, ncid) )

      CALL check( nf90_inq_varid(ncid, "PCT_PFT"        , pftvid     ) )
      CALL check( nf90_inq_varid(ncid, "MONTHLY_PFT_LAI", laivid     ) )
      CALL check( nf90_inq_varid(ncid, "MONTHLY_PFT_SAI", saivid     ) )
      CALL check( nf90_inq_varid(ncid, "HTOP_PFT"       , htop_pftvid) )

      CALL check( nf90_get_var(ncid, pftvid     , pct_pft  ) )
      CALL check( nf90_get_var(ncid, saivid     , sai_pft  ) )
      CALL check( nf90_get_var(ncid, laivid     , lai_pft  ) )
      CALL check( nf90_get_var(ncid, htop_pftvid, htop_pft ) )

      CALL check( nf90_close(ncid) )
   ELSE
      ! 如果没有以上数据，则需要利用MakeGlobalSurface.F90的代码计算
      ! 可换成CALL Subroutine
      fileNotExists = .TRUE.

      ALLOCATE( lcdata  (nxy, nxy) )
      ALLOCATE( pcrop   (nxy, nxy) )
      ALLOCATE( pwetland(nxy, nxy) )
      ALLOCATE( pwater  (nxy, nxy) )
      ALLOCATE( pglacier(nxy, nxy) )
      ALLOCATE( ppft    (nxy, nxy, npft) )
      ALLOCATE( pftlai  (nxy, nxy, npft, mon) )
      ALLOCATE( pftsai  (nxy, nxy, npft, mon) )
   ENDIF

   ! output lat/lon
   DO i = 1, nyo
      lats1(i) = -90. + (i-1)*y_delta
      latn1(i) = -90. + i*y_delta
   ENDDO

   DO i = 1, nxo
      lonw1(i) = -180. + (i-1)*x_delta
      lone1(i) = -180. + i*x_delta
   ENDDO

   DO i = 1, nxo
      lonso(i) = -180. + i*x_delta - 0.5*x_delta
   ENDDO

   DO i = 1, nyo
      latso(i) =  90. - i*y_delta + 0.5*y_delta
   ENDDO

   ! model gird area
   DO i = 1, nyo
      dx = (lone1(1)-lonw1(1))*deg2rad
      dy = sin(latn1(i)*deg2rad) - sin(lats1(i)*deg2rad)
      area(:,i) = dx*dy*re*re
   ENDDO

   OPEN(11,FILE=REGFILE)
   OPEN(12,FILE=REGFILE)
   
   DO WHILE(.TRUE.)
      ! process global 500m raw data
      ! read 500m surface data(pct_tree/pct_water)
      ! NOTE: GEDI data not use yet
      PRINT *, "*** processing 500m raw data ***"
      PRINT *, "    reading data"
      READ(11,*,END=100) reg
      READ(12,*,END=101) creg
      filename = TRIM(SRF_DIR)//'RG_'//TRIM(creg(1))//'_'//&
                 TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.SRF'//TRIM(year)//'.nc' !TRIM(year)//'.nc'

      PRINT*, filename
      CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )

      CALL check( nf90_inq_varid(ncid, "lat"        , hlat_vid  ) )
      CALL check( nf90_inq_varid(ncid, "lon"        , hlon_vid  ) )
      CALL check( nf90_inq_varid(ncid, "PCT_Tree"   , gfcc_tcvid) )
      !CALL check( nf90_inq_varid(ncid, "Hgt_Tree"   , gedi_thvid) )
      CALL check( nf90_inq_varid(ncid, "PCT_Water"  , gl30_wtvid) )

      CALL check( nf90_get_var(ncid, hlat_vid  , hlat   ) )
      CALL check( nf90_get_var(ncid, hlon_vid  , hlon   ) )
      CALL check( nf90_get_var(ncid, gfcc_tcvid, gfcc_tc) )
      !CALL check( nf90_get_var(ncid, gedi_thvid, gedi_th) )
      CALL check( nf90_get_var(ncid, gl30_wtvid, gl30_wt) )

      CALL check( nf90_close(ncid) )

      ! read NCAR 1km urban class and region id
      PRINT *, "    reading 1km urban class data"
      filename = TRIM(SRF_DIR)//'RG_'//TRIM(creg(1))//'_'//&
                 TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.NCAR.nc'
      CALL check( nf90_open(filename, nf90_nowrite, ncid) )

      ! CALL check( nf90_inq_varid(ncid, "lat"                , urlat_vid ) )
      ! CALL check( nf90_inq_varid(ncid, "lon"                , urlon_vid ) )
      CALL check( nf90_inq_varid(ncid, "URBAN_DENSITY_CLASS", ur_clssvid) )
      CALL check( nf90_inq_varid(ncid, "REGION_ID"          , ur_rgvid  ) )

      ! CALL check( nf90_get_var(ncid, urlat_vid   , urlat ) )
      ! CALL check( nf90_get_var(ncid, urlon_vid   , urlon ) )
      CALL check( nf90_get_var(ncid, ur_clssvid  , urden ) )
      CALL check( nf90_get_var(ncid, ur_rgvid    , urrgid) )

      CALL check( nf90_close(ncid) )

      ! read 500m modis urban data
      PRINT*, "    reading MOD URBAN data"
      filename = TRIM(RAW_DIR)//'RG_'//TRIM(creg(1))//'_'//&
                 TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.MOD'//TRIM(year)//'.nc'
      ! PRINT*, filename
      IF (fileExists) THEN
         ! 当全球PFTs数据已经生成时，500m数据只需要读取以下两个变量
         CALL check( nf90_open(filename , nf90_nowrite    , ncid       ) )

         CALL check( nf90_inq_varid(ncid, "PCT_URBAN"     , upftvid    ) )
         CALL check( nf90_inq_varid(ncid, "HTOP"          , gedi_thvid ) )

         CALL check( nf90_get_var(ncid, upftvid    , modur   ) )
         CALL check( nf90_get_var(ncid, gedi_thvid , gedi_th ) )

         CALL check( nf90_close(ncid) )
      ELSE
         ! 当全球PFTs数据不存在时， 需要完整读取MOD数据进行聚合计算PFTs等
         ! 用于URBAN_LAI计算
         CALL check( nf90_open(filename , nf90_nowrite    , ncid       ) )

         CALL check( nf90_inq_varid(ncid, "PCT_URBAN"     , upftvid    ) )
         CALL check( nf90_inq_varid(ncid, "HTOP"          , gedi_thvid ) )
         CALL check( nf90_inq_varid(ncid, "PCT_PFT"       , pftvid     ) )
         CALL check( nf90_inq_varid(ncid, "MONTHLY_LAI"   , laivid     ) )
         CALL check( nf90_inq_varid(ncid, "MONTHLY_SAI"   , saivid     ) )   
         CALL check( nf90_inq_varid(ncid, "LC"            , lcvid      ) )
         CALL check( nf90_inq_varid(ncid, "PCT_CROP"      , pcropvid   ) )
         CALL check( nf90_inq_varid(ncid, "PCT_WETLAND"   , pwetlandvid) )
         CALL check( nf90_inq_varid(ncid, "PCT_WATER"     , pwatervid  ) )
         CALL check( nf90_inq_varid(ncid, "PCT_GLACIER"   , pglaciervid) )

         CALL check( nf90_get_var(ncid, lcvid      , lcdata  ) )
         CALL check( nf90_get_var(ncid, pcropvid   , pcrop   ) )
         CALL check( nf90_get_var(ncid, pwetlandvid, pwetland) )
         CALL check( nf90_get_var(ncid, pglaciervid, pglacier) )
         CALL check( nf90_get_var(ncid, pwatervid  , pwater  ) )
         CALL check( nf90_get_var(ncid, pftvid     , ppft    ) )
         CALL check( nf90_get_var(ncid, saivid     , pftsai  ) )
         CALL check( nf90_get_var(ncid, laivid     , pftlai  ) )
         CALL check( nf90_get_var(ncid, upftvid    , modur   ) )
         CALL check( nf90_get_var(ncid, gedi_thvid , gedi_th ) )

         CALL check( nf90_close(ncid) )
      ENDIF

      ! calculate the edge of small grids(500m)
      DO i = 1, nxy
         hlats(i) = reg(1) - i*sdelta
         hlatn(i) = reg(1) - (i-1)*sdelta
         hlonw(i) = reg(2) + (i-1)*sdelta
         hlone(i) = reg(2) + i*sdelta
      ENDDO

      DO i = 1, nxy
         dx = (hlone(1)-hlonw(1))*deg2rad
         dy = sin(hlatn(i)*deg2rad) - sin(hlats(i)*deg2rad)
         harea(:,i) = dx*dy*re*re
      ENDDO

      PRINT *, "    aggregating data"
      PRINT *, "    aggregating urban class data"
!!$OMP PARALLEL DO NUM_THREADS(92) &
!!$OMP PRIVATE(i,j,io,jo,uxid,m,lc,wgt,sumpct,ip)
      DO i = 1, 1200
         DO j = 1, 1200
            ! calculate io, jo
            !io = NINT((hlats(i)+sdelta/2+ 90.)/y_delta+0.5)
            !io = nyo+1-io
            io = NINT((90.-(hlats(i)+sdelta/2))/y_delta+0.5)
            jo = NINT((hlonw(j)+sdelta/2+180.)/x_delta+0.5)

            
            ! 聚合Simard树高
            ! 加权：
            ! 粗网格树高=粗网格树高+simard树高*1km格点面积
            ! 加权系数：粗网格格点面积
            IF (gedi_th(j,i) > 0) THEN
               hgt(jo,io) = hgt(jo,io) + gedi_th(j,i)*harea(j,i)
               avg(jo,io) = avg(jo,io) + harea(j,i)
            ENDIF

            IF (urden(j,i) > 0) THEN
             ! Tall-Building-Distinc urban
               IF (urden(j,i) == 1) THEN
            
               ! 加权：
               ! 粗网格城市水体(植被)覆盖度=粗网格城市水体(植被)覆盖度+500m城市格点水体(植被)覆盖度*500m城市格点面积
               ! 加权系数；粗网格城市格点面积
                  IF (gl30_wt(j,i) > 0.) THEN
                     urwt (jo,io,1) = urwt (jo,io,1) + gl30_wt(j,i)*harea(j,i)
                  ENDIF
                  IF (gfcc_tc(j,i) >0) THEN
                     tc   (jo,io,1) = tc   (jo,io,1) + gfcc_tc(j,i)*harea(j,i)
                     ! 树高加权
                     ! 粗网格城市树高=粗网格城市树高+500m城市格点植被覆盖度*城市格点树高*城市格点面积
                     ! 加权系数:城市格点植被覆盖度*城市面积
                     IF (gedi_th(j,i) > 0) THEN
                        htop(jo,io,1) = htop(jo,io,1) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)
                        wgt_top(jo,io,1) = wgt_top(jo,io,1) + gfcc_tc(j,i)*harea(j,i)
                     ENDIF
                  ENDIF

                  ur_dc    (jo,io,1) = ur_dc(jo,io,1) + harea(j,i)
                  uxid               = urrgid(j,i)
                  hwr_can  (jo,io,1) = hwr_can  (jo,io,1) + hwrcan  (1,uxid)*harea(j,i)
                  wt_rf    (jo,io,1) = wt_rf    (jo,io,1) + wtrf    (1,uxid)*harea(j,i)
                  wt_rd    (jo,io,1) = wt_rd    (jo,io,1) + wtrd    (1,uxid)*harea(j,i)
                  em_rf    (jo,io,1) = em_rf    (jo,io,1) + emrf    (1,uxid)*harea(j,i)
                  em_wl    (jo,io,1) = em_wl    (jo,io,1) + emwl    (1,uxid)*harea(j,i)
                  em_imrd  (jo,io,1) = em_imrd  (jo,io,1) + emimrd  (1,uxid)*harea(j,i)
                  em_perd  (jo,io,1) = em_perd  (jo,io,1) + emperd  (1,uxid)*harea(j,i)
                  th_rf    (jo,io,1) = th_rf    (jo,io,1) + thrf    (1,uxid)*harea(j,i)
                  th_wl    (jo,io,1) = th_wl    (jo,io,1) + thwl    (1,uxid)*harea(j,i)
                  tb_min   (jo,io,1) = tb_min   (jo,io,1) + tbmin   (1,uxid)*harea(j,i)
                  tb_max   (jo,io,1) = tb_max   (jo,io,1) + tbmax   (1,uxid)*harea(j,i)
                  ht_rf    (jo,io,1) = ht_rf    (jo,io,1) + htrf    (1,uxid)*harea(j,i)

                  alb_rf  (jo,io,1,:,:) = alb_rf  (jo,io,1,:,:) + albrf  (1,uxid,:,:)*harea(j,i)
                  alb_wl  (jo,io,1,:,:) = alb_wl  (jo,io,1,:,:) + albwl  (1,uxid,:,:)*harea(j,i)
                  alb_imrd(jo,io,1,:,:) = alb_imrd(jo,io,1,:,:) + albimrd(1,uxid,:,:)*harea(j,i)
                  alb_perd(jo,io,1,:,:) = alb_perd(jo,io,1,:,:) + albperd(1,uxid,:,:)*harea(j,i)

                  tk_rf  (jo,io,1,:) = tk_rf  (jo,io,1,:) + tkrf  (1,uxid,:)*harea(j,i)
                  tk_wl  (jo,io,1,:) = tk_wl  (jo,io,1,:) + tkwl  (1,uxid,:)*harea(j,i)
                  DO m = 1, 10
                     ! tkimrd与cvimrd有缺省值，计算需要跳过
                     IF (tkimrd(1,uxid,m) .ne. -999) THEN
                        tk_imrd(jo,io,1,m) = tk_imrd(jo,io,1,m) + tkimrd(1,uxid,m)*harea(j,i)
                     ENDIF
                     IF (cvimrd(1,uxid,m) .ne. -999.) THEN
                        cv_imrd(jo,io,1,m) = cv_imrd(jo,io,1,m) + cvimrd(1,uxid,m)*harea(j,i)
                     ENDIF
                  ENDDO
                  cv_rf  (jo,io,1,:) = cv_rf  (jo,io,1,:) + cvrf  (1,uxid,:)*harea(j,i)
                  cv_wl  (jo,io,1,:) = cv_wl  (jo,io,1,:) + cvwl  (1,uxid,:)*harea(j,i)
               ENDIF

               ! Hight-Density urban
               IF (urden(j,i) == 2 ) THEN
                  IF (gl30_wt(j,i) > 0.) THEN
                     urwt (jo,io,2) = urwt (jo,io,2) + gl30_wt(j,i)*harea(j,i)
                  ENDIF
                  IF (gfcc_tc(j,i) > 0.) THEN 
                     tc   (jo,io,2) = tc   (jo,io,2) + gfcc_tc(j,i)*harea(j,i)
                     IF (gedi_th(j,i) > 0) THEN
                        htop(jo,io,2) = htop(jo,io,2) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)
                        wgt_top(jo,io,2) = wgt_top(jo,io,2) + gfcc_tc(j,i)*harea(j,i)
                     ENDIF
                  ENDIF

                  ur_dc(jo,io,2) = ur_dc(jo,io,2) + harea(j,i)
                  uxid               = urrgid(j,i)
                  hwr_can  (jo,io,2) = hwr_can  (jo,io,2) + hwrcan  (2,uxid)*harea(j,i)
                  wt_rf    (jo,io,2) = wt_rf    (jo,io,2) + wtrf    (2,uxid)*harea(j,i)
                  wt_rd    (jo,io,2) = wt_rd    (jo,io,2) + wtrd    (2,uxid)*harea(j,i)
                  em_rf    (jo,io,2) = em_rf    (jo,io,2) + emrf    (2,uxid)*harea(j,i)
                  em_wl    (jo,io,2) = em_wl    (jo,io,2) + emwl    (2,uxid)*harea(j,i)
                  em_imrd  (jo,io,2) = em_imrd  (jo,io,2) + emimrd  (2,uxid)*harea(j,i)
                  em_perd  (jo,io,2) = em_perd  (jo,io,2) + emperd  (2,uxid)*harea(j,i)
                  th_rf    (jo,io,2) = th_rf    (jo,io,2) + thrf    (2,uxid)*harea(j,i)
                  th_wl    (jo,io,2) = th_wl    (jo,io,2) + thwl    (2,uxid)*harea(j,i)
                  tb_min   (jo,io,2) = tb_min   (jo,io,2) + tbmin   (2,uxid)*harea(j,i)
                  tb_max   (jo,io,2) = tb_max   (jo,io,2) + tbmax   (2,uxid)*harea(j,i)
                  ht_rf    (jo,io,2) = ht_rf    (jo,io,2) + htrf    (2,uxid)*harea(j,i)

                  alb_rf  (jo,io,2,:,:) = alb_rf  (jo,io,2,:,:) + albrf  (2,uxid,:,:)*harea(j,i)
                  alb_wl  (jo,io,2,:,:) = alb_wl  (jo,io,2,:,:) + albwl  (2,uxid,:,:)*harea(j,i)
                  alb_imrd(jo,io,2,:,:) = alb_imrd(jo,io,2,:,:) + albimrd(2,uxid,:,:)*harea(j,i)
                  alb_perd(jo,io,2,:,:) = alb_perd(jo,io,2,:,:) + albperd(2,uxid,:,:)*harea(j,i)

                  tk_rf  (jo,io,2,:) = tk_rf  (jo,io,2,:) + tkrf  (2,uxid,:)*harea(j,i)
                  tk_wl  (jo,io,2,:) = tk_wl  (jo,io,2,:) + tkwl  (2,uxid,:)*harea(j,i)
                  DO m = 1, 10
                     IF (tkimrd(2,uxid,m) .ne. -999) THEN
                        tk_imrd(jo,io,2,m) = tk_imrd(jo,io,2,m) + tkimrd(2,uxid,m)*harea(j,i)
                     ENDIF
                     IF (cvimrd(2,uxid,m) .ne. -999.) THEN
                        cv_imrd(jo,io,2,m) = cv_imrd(jo,io,2,m) + cvimrd(2,uxid,m)*harea(j,i)
                     ENDIF
                  ENDDO
                  cv_rf  (jo,io,2,:) = cv_rf  (jo,io,2,:) + cvrf  (2,uxid,:)*harea(j,i)
                  cv_wl  (jo,io,2,:) = cv_wl  (jo,io,2,:) + cvwl  (2,uxid,:)*harea(j,i)
               ENDIF

               ! Medium-Density urban
               IF (urden(j,i) == 3 ) THEN
                  IF (gl30_wt(j,i) > 0.) THEN
                     urwt (jo,io,3) = urwt (jo,io,3) + gl30_wt(j,i)*harea(j,i)
                  ENDIF
                  IF (gfcc_tc(j,i) > 0.) THEN
                     tc   (jo,io,3) = tc   (jo,io,3) + gfcc_tc(j,i)*harea(j,i)
                     IF (gedi_th(j,i) > 0) THEN
                        htop(jo,io,3) = htop(jo,io,3) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)
                        wgt_top(jo,io,3) = wgt_top(jo,io,3) + harea(j,i)*gfcc_tc(j,i)
                     ENDIF
                  ENDIF
                  ur_dc(jo,io,3) = ur_dc(jo,io,3) + harea(j,i)
                  uxid               = urrgid(j,i)
                  hwr_can  (jo,io,3) = hwr_can  (jo,io,3) + hwrcan  (3,uxid)*harea(j,i)
                  wt_rf    (jo,io,3) = wt_rf    (jo,io,3) + wtrf    (3,uxid)*harea(j,i)
                  wt_rd    (jo,io,3) = wt_rd    (jo,io,3) + wtrd    (3,uxid)*harea(j,i)
                  em_rf    (jo,io,3) = em_rf    (jo,io,3) + emrf    (3,uxid)*harea(j,i)
                  em_wl    (jo,io,3) = em_wl    (jo,io,3) + emwl    (3,uxid)*harea(j,i)
                  em_imrd  (jo,io,3) = em_imrd  (jo,io,3) + emimrd  (3,uxid)*harea(j,i)
                  em_perd  (jo,io,3) = em_perd  (jo,io,3) + emperd  (3,uxid)*harea(j,i)
                  th_rf    (jo,io,3) = th_rf    (jo,io,3) + thrf    (3,uxid)*harea(j,i)
                  th_wl    (jo,io,3) = th_wl    (jo,io,3) + thwl    (3,uxid)*harea(j,i)
                  tb_min   (jo,io,3) = tb_min   (jo,io,3) + tbmin   (3,uxid)*harea(j,i)
                  tb_max   (jo,io,3) = tb_max   (jo,io,3) + tbmax   (3,uxid)*harea(j,i)
                  ht_rf    (jo,io,3) = ht_rf    (jo,io,3) + htrf    (3,uxid)*harea(j,i)

                  alb_rf  (jo,io,3,:,:) = alb_rf  (jo,io,3,:,:) + albrf  (3,uxid,:,:)*harea(j,i)
                  alb_wl  (jo,io,3,:,:) = alb_wl  (jo,io,3,:,:) + albwl  (3,uxid,:,:)*harea(j,i)
                  alb_imrd(jo,io,3,:,:) = alb_imrd(jo,io,3,:,:) + albimrd(3,uxid,:,:)*harea(j,i)
                  alb_perd(jo,io,3,:,:) = alb_perd(jo,io,3,:,:) + albperd(3,uxid,:,:)*harea(j,i)

                  tk_rf  (jo,io,3,:) = tk_rf  (jo,io,3,:) + tkrf  (3,uxid,:)*harea(j,i)
                  tk_wl  (jo,io,3,:) = tk_wl  (jo,io,3,:) + tkwl  (3,uxid,:)*harea(j,i)
                  DO m = 1, 10
                     IF (tkimrd(3,uxid,m) .ne. -999.) THEN
                        tk_imrd(jo,io,3,m) = tk_imrd(jo,io,3,m) + tkimrd(3,uxid,m)*harea(j,i)
                     ENDIF
                     IF (cvimrd(3,uxid,m) .ne. -999.) THEN
                        cv_imrd(jo,io,3,m) = cv_imrd(jo,io,3,m) + cvimrd(3,uxid,m)*harea(j,i)
                     ENDIF
                  ENDDO
                  cv_rf  (jo,io,3,:) = cv_rf  (jo,io,3,:) + cvrf  (3,uxid,:)*harea(j,i)
                  cv_wl  (jo,io,3,:) = cv_wl  (jo,io,3,:) + cvwl  (3,uxid,:)*harea(j,i)
               ENDIF
            ENDIF

            ! 根据MODIS城市覆盖对城市格点补充，并将其归类为MD urban
            IF (urden(j,i)<=0 .and. modur(j,i)>0) THEN
               IF (gl30_wt(j,i) > 0.) THEN
                  urwt (jo,io,3) = urwt (jo,io,3) + gl30_wt(j,i)*harea(j,i)*modur(j,i)/100
               ENDIF
               IF (gfcc_tc(j,i) > 0.) THEN
                  tc   (jo,io,3) = tc   (jo,io,3) + gfcc_tc(j,i)*harea(j,i)*modur(j,i)/100
                  IF (gedi_th(j,i) > 0) THEN
                     htop(jo,io,3) = htop(jo,io,3) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)*modur(j,i)/100
                     wgt_top(jo,io,3) = wgt_top(jo,io,3) + harea(j,i)*gfcc_tc(j,i)*modur(j,i)/100
                  ENDIF
               ENDIF

               ur_dc(jo,io,3) = ur_dc(jo,io,3) + harea(j,i)*modur(j,i)/100
               uxid           = urrgid(j,i)
                  
               ! 部分格点MODIS与NCAR不一致(NCAR没有城市ID)，因此通过距离MODIS格点最近的NCAR城市ID赋值
               IF (reg(1)==-45 .and. reg(3)==-50 .and. reg(2)==65 .and. reg(4)==70) THEN
                  uxid = 30
               ENDIF
               ! 城市建筑属性聚合
               ! 加权：
               ! 粗网格城市属性=粗网格城市属性+细网格城市属性*细网格面积*MODIS_PCT_URBAN
               ! 加权系数：细网格面积(ur_dc)
               hwr_can  (jo,io,3) = hwr_can  (jo,io,3) + hwrcan  (3,uxid)*harea(j,i)*modur(j,i)/100
               wt_rf    (jo,io,3) = wt_rf    (jo,io,3) + wtrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               wt_rd    (jo,io,3) = wt_rd    (jo,io,3) + wtrd    (3,uxid)*harea(j,i)*modur(j,i)/100
               em_rf    (jo,io,3) = em_rf    (jo,io,3) + emrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               em_wl    (jo,io,3) = em_wl    (jo,io,3) + emwl    (3,uxid)*harea(j,i)*modur(j,i)/100
               em_imrd  (jo,io,3) = em_imrd  (jo,io,3) + emimrd  (3,uxid)*harea(j,i)*modur(j,i)/100
               em_perd  (jo,io,3) = em_perd  (jo,io,3) + emperd  (3,uxid)*harea(j,i)*modur(j,i)/100
               th_rf    (jo,io,3) = th_rf    (jo,io,3) + thrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               th_wl    (jo,io,3) = th_wl    (jo,io,3) + thwl    (3,uxid)*harea(j,i)*modur(j,i)/100
               tb_min   (jo,io,3) = tb_min   (jo,io,3) + tbmin   (3,uxid)*harea(j,i)*modur(j,i)/100
               tb_max   (jo,io,3) = tb_max   (jo,io,3) + tbmax   (3,uxid)*harea(j,i)*modur(j,i)/100
               ht_rf    (jo,io,3) = ht_rf    (jo,io,3) + htrf    (3,uxid)*harea(j,i)*modur(j,i)/100

               alb_rf  (jo,io,3,:,:) = alb_rf  (jo,io,3,:,:) + albrf  (3,uxid,:,:)*harea(j,i)*modur(j,i)/100
               alb_wl  (jo,io,3,:,:) = alb_wl  (jo,io,3,:,:) + albwl  (3,uxid,:,:)*harea(j,i)*modur(j,i)/100
               alb_imrd(jo,io,3,:,:) = alb_imrd(jo,io,3,:,:) + albimrd(3,uxid,:,:)*harea(j,i)*modur(j,i)/100
               alb_perd(jo,io,3,:,:) = alb_perd(jo,io,3,:,:) + albperd(3,uxid,:,:)*harea(j,i)*modur(j,i)/100

               tk_rf  (jo,io,3,:) = tk_rf  (jo,io,3,:) + tkrf  (3,uxid,:)*harea(j,i)*modur(j,i)/100
               tk_wl  (jo,io,3,:) = tk_wl  (jo,io,3,:) + tkwl  (3,uxid,:)*harea(j,i)*modur(j,i)/100
               DO m = 1, 10
                  IF (tkimrd(3,uxid,m) .ne. -999.) THEN
                     tk_imrd(jo,io,3,m) = tk_imrd(jo,io,3,m) + tkimrd(3,uxid,m)*harea(j,i)*modur(j,i)/100
                  ENDIF
                  IF (cvimrd(3,uxid,m) .ne. -999.) THEN
                     cv_imrd(jo,io,3,m) = cv_imrd(jo,io,3,m) + cvimrd(3,uxid,m)*harea(j,i)*modur(j,i)/100
                  ENDIF
               ENDDO
               cv_rf  (jo,io,3,:) = cv_rf  (jo,io,3,:) + cvrf  (3,uxid,:)*harea(j,i)*modur(j,i)/100
               cv_wl  (jo,io,3,:) = cv_wl  (jo,io,3,:) + cvwl  (3,uxid,:)*harea(j,i)*modur(j,i)/100
            ENDIF

            ! MakeGlobalSurface.F90, 可改成CALL Subroutine
            ! https://github.com/sunan126/MKSRF.git MakeGlobalSurface.F90
            IF (fileNotExists) THEN
               lc = lcdata(j,i)
               IF (lc /= 0) THEN

                  ! aggregate land
                  pct_land(jo,io) = pct_land(jo,io) + harea(j,i)
                  ! aggregate on PFT level
                  ! 隐含的假设: 水体、冰川、城市、湿地和作物可由专业化的数据或高分辨率数据获得
                  ! 下面的代码具有兼容性，但目前并没有进行如上更新
                  ! exclude water, ice, urban, wetland, and crop
                  wgt = max(0., 1.-pwater(j,i)-pglacier(j,i)-modur(j,i)-pwetland(j,i)-pcrop(j,i))

                  ! adjust nature PFT's total percentage to 100%
                  sumpct = sum(ppft(j,i,1:(npft-1)))
                  IF (sumpct > 0.) THEN
                     ppft(j,i,1:(npft-1)) = ppft(j,i,1:(npft-1))/sumpct
                  ENDIF

                  DO ip = 1, npft-1
                     pct_pft(jo,io,ip) = pct_pft(jo,io,ip) + &
                        harea(j,i)*wgt*ppft(j,i,ip)
                     htop_pft(jo,io,ip) = htop_pft(jo,io,ip) + &
                        harea(j,i)*wgt*ppft(j,i,ip)*gedi_th(j,i)
                     lai_pft(jo,io,ip,:) = lai_pft(jo,io,ip,:) + &
                        harea(j,i)*wgt*ppft(j,i,ip)*pftlai(j,i,ip,:)
                     sai_pft(jo,io,ip,:) = sai_pft(jo,io,ip,:) + &
                        harea(j,i)*wgt*ppft(j,i,ip)*pftsai(j,i,ip,:)
                  ENDDO

                  ! for crop (npft=16, is crop)
                  pct_pft(jo,io,npft) = pct_pft(jo,io,npft) + &
                     harea(j,i)*ppft(j,i,npft)
                  htop_pft(jo,io,npft) = htop_pft(jo,io,npft) + &
                     harea(j,i)*ppft(j,i,npft)*gedi_th(j,i)
                  lai_pft(jo,io,npft,:) = lai_pft(jo,io,npft,:) + &
                     harea(j,i)*ppft(j,i,npft)*pftlai(j,i,npft,:)
                  sai_pft(jo,io,npft,:) = sai_pft(jo,io,npft,:) + &
                     harea(j,i)*ppft(j,i,npft)*pftsai(j,i,npft,:)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!!$OMP END PARALLEL DO
   ENDDO

   100 close(11)
   101 close(12)

   DO i = 1, nyo
      DO j = 1, nxo
         DO k = 1, 3
               IF (ur_dc(j,i,k) > 0) THEN
                  ! calculate urban tree cover
                  pct_tc  (j,i,k) = tc  (j,i,k) / ur_dc(j,i,k) !* 100
                  ! calculate urban water cover
                  pct_urwt(j,i,k) = urwt(j,i,k) / ur_dc(j,i,k) !* 100
                  IF (wgt_top(j,i,k) > 0.) THEN
                  ! calculate urban tree height
                     htop_ur (j,i,k) = htop(j,i,k) / wgt_top(j,i,k)!tc   (j,i,k-1)
                  ENDIF
               
                  hwr_can  (j,i,k) = hwr_can  (j,i,k) / ur_dc(j,i,k)
                  wt_rf    (j,i,k) = wt_rf    (j,i,k) / ur_dc(j,i,k)
                  wt_rd    (j,i,k) = wt_rd    (j,i,k) / ur_dc(j,i,k)
                  em_rf    (j,i,k) = em_rf    (j,i,k) / ur_dc(j,i,k)
                  em_wl    (j,i,k) = em_wl    (j,i,k) / ur_dc(j,i,k)
                  em_imrd  (j,i,k) = em_imrd  (j,i,k) / ur_dc(j,i,k)
                  em_perd  (j,i,k) = em_perd  (j,i,k) / ur_dc(j,i,k)
                  th_rf    (j,i,k) = th_rf    (j,i,k) / ur_dc(j,i,k)
                  th_wl    (j,i,k) = th_wl    (j,i,k) / ur_dc(j,i,k)
                  tb_min   (j,i,k) = tb_min   (j,i,k) / ur_dc(j,i,k)
                  tb_max   (j,i,k) = tb_max   (j,i,k) / ur_dc(j,i,k)
                  ht_rf    (j,i,k) = ht_rf    (j,i,k) / ur_dc(j,i,k)
    
                  alb_rf  (j,i,k,:,:) = alb_rf  (j,i,k,:,:) / ur_dc(j,i,k)
                  alb_wl  (j,i,k,:,:) = alb_wl  (j,i,k,:,:) / ur_dc(j,i,k)
                  alb_imrd(j,i,k,:,:) = alb_imrd(j,i,k,:,:) / ur_dc(j,i,k)
                  alb_perd(j,i,k,:,:) = alb_perd(j,i,k,:,:) / ur_dc(j,i,k)
    
                  tk_rf  (j,i,k,:) = tk_rf  (j,i,k,:) / ur_dc(j,i,k)
                  tk_wl  (j,i,k,:) = tk_wl  (j,i,k,:) / ur_dc(j,i,k)
                  DO m = 1, 10
                     IF (tk_imrd(j,i,k,m) >= 0.) THEN
                        tk_imrd(j,i,k,m) = tk_imrd(j,i,k,m) / ur_dc(j,i,k)
                     ENDIF
                     IF (cv_imrd(j,i,k,m) >= 0.) THEN
                        cv_imrd(j,i,k,m) = cv_imrd(j,i,k,m) / ur_dc(j,i,k)
                     ENDIF
                  ENDDO
                  cv_rf  (j,i,k,:) = cv_rf  (j,i,k,:) / ur_dc(j,i,k)
                  cv_wl  (j,i,k,:) = cv_wl  (j,i,k,:) / ur_dc(j,i,k)
               ENDIF

               IF (fileNotExists) THEN
                  DO ip = 1, npft
                     IF (pct_pft(j,i,ip) > 0) THEN
                        lai_pft(j,i,ip,:) =  lai_pft(j,i,ip,:)/pct_pft(j,i,ip)
                        sai_pft(j,i,ip,:) =  sai_pft(j,i,ip,:)/pct_pft(j,i,ip)
                        htop_pft(j,i,ip)   = htop_pft(j,i,ip)  /pct_pft(j,i,ip)
                     ENDIF
                  ENDDO

                  IF (pct_land(j,i) > 0) THEN
                  ! PFT level
                     pct_pft    (j,i,:) = pct_pft(j,i,:) / pct_land(j,i) * 100.
                  ENDIF
               ENDIF
         ENDDO

         IF (avg(j,i) > 0) THEN
            hgt(j,i) = hgt(j,i) / avg(j,i)
         ENDIF 

         sumur = ur_dc(j,i,1) + ur_dc(j,i,2) + ur_dc(j,i,3)
         IF (sumur > 0.) THEN
            ! calculate 3 types urban cover, sum(pct_ur(j,i,:))=100
            pct_ur(j,i,1) = ur_dc(j,i,1) / sumur * 100.
            pct_ur(j,i,2) = ur_dc(j,i,2) / sumur * 100.
            pct_ur(j,i,3) = ur_dc(j,i,3) / sumur * 100.
         ENDIF

         ! 检查3类城市是否超过100%
         IF (sum(pct_ur(j,i,1:3)) > 1e-6 .and. abs(sum(pct_ur(j,i,1:3))-100) > 1e-3) THEN
            PRINT *, 'urban_pct > 100'
            PRINT *, pct_ur(j,i,1:3)
         ENDIF
      ENDDO
   ENDDO   
   
   PRINT *, "********************************"   
   DO i = 1, nyo 
      DO j = 1, nxo
         DO k =1, 3
            ! check for htop_ur of urban grid
            ! 如果该城市格点有植被覆盖却没有树高，则将聚合过程中生成的该格点的htop数据
            ! 赋值为城市树高
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               htop_ur(j,i,k) = hgt(j,i)
            ENDIF

            ! 如果经过上一步仍没有树高数据
            ! 则将该点同纬度的所有树高求平均赋值城市树高
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               DO p = 1, nxo
                  IF (hgt(p,i) > 0) THEN
                     sumth = sumth + hgt(p,i)
                     cont  = cont  + 1
                  ENDIF
               ENDDO

               IF (cont > 0) THEN
                  htop_ur(j,i,k) = sumth/cont
               ENDIF

               cont  = 0
               sumth = 0
            ENDIF
         ENDDO

         ! calculate urban lai
         !-----------------------------------------------------
         !             ___npft
         ! urban_lai = \   pct_pft*pft_lai*(htop_gedi/htop_pft)
         !             /__ 1
         !-----------------------------------------------------
         ! 只需要考虑树(1-8类PFTs)
         DO m = 2, 9
            IF (htop_pft(j,i,m) > 0) THEN
               fac(1) = htop_ur(j,i,1) / htop_pft(j,i,m)
               fac(2) = htop_ur(j,i,2) / htop_pft(j,i,m)
               fac(3) = htop_ur(j,i,3) / htop_pft(j,i,m)

               wgt = pct_pft(j,i,m)/sum(pct_pft(j,i,2:9))               
               DO k = 1, 12       
                  ur_lai (j,i,1,k) = ur_lai(j,i,1,k)+wgt*lai_pft(j,i,m,k)*min(1.5,fac(1))
                  ur_lai (j,i,2,k) = ur_lai(j,i,2,k)+wgt*lai_pft(j,i,m,k)*min(1.5,fac(2))
                  ur_lai (j,i,3,k) = ur_lai(j,i,3,k)+wgt*lai_pft(j,i,m,k)*min(1.5,fac(3))
                  ur_sai (j,i,1,k) = ur_sai(j,i,1,k)+wgt*sai_pft(j,i,m,k)*min(1.5,fac(1))
                  ur_sai (j,i,2,k) = ur_sai(j,i,2,k)+wgt*sai_pft(j,i,m,k)*min(1.5,fac(2))
                  ur_sai (j,i,3,k) = ur_sai(j,i,3,k)+wgt*sai_pft(j,i,m,k)*min(1.5,fac(3))
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   PRINT *, "*****************************************"

   DO i = 1, 2
      n_nr(i) = i
   ENDDO

   DO i = 1, 2
      n_ns(i) = i
   ENDDO

   DO i = 1, 3
      n_den(i) = i
   ENDDO

   DO i = 1, 10
      n_ulev(i) = i
   ENDDO

   DO i = 1, 12
      n_mon(i) = i
   ENDDO

   
   CALL check( nf90_create(TRIM(OUT_DIR)//"colm_urban_data_modis_v5.6_"//trim(year)//".nc", NF90_NETCDF4, ncid) )

   CALL check( nf90_def_dim(ncid, "lat"     , nyo     , lat_dimid ) )
   CALL check( nf90_def_dim(ncid, "lon"     , nxo     , lon_dimid ) )
   CALL check( nf90_def_dim(ncid, "numsolar", ns      , ns_dimid  ) )
   CALL check( nf90_def_dim(ncid, "numrad"  , nr      , nr_dimid  ) )
   CALL check( nf90_def_dim(ncid, "ulev"    , ulev    , ulev_dimid) )
   CALL check( nf90_def_dim(ncid, "density" , den_clss, den_dimid ) )
   CALL check( nf90_def_dim(ncid, "month"   , mon     , mon_dimid ) )

   CALL check( nf90_def_var(ncid, "lat"     , NF90_FLOAT, lat_dimid , lat_vid ) )
   CALL check( nf90_def_var(ncid, "lon"     , NF90_FLOAT, lon_dimid , lon_vid ) )
   CALL check( nf90_def_var(ncid, "numsolar", NF90_INT  , ns_dimid  , ns_vid  ) )
   CALL check( nf90_def_var(ncid, "numrad"  , NF90_INT  , nr_dimid  , nr_vid  ) )
   CALL check( nf90_def_var(ncid, "ulev"    , NF90_INT  , ulev_dimid, lev_vid ) )
   CALL check( nf90_def_var(ncid, "density" , NF90_INT  , den_dimid , den_vid ) )
   CALL check( nf90_def_var(ncid, "month"   , NF90_INT  , mon_dimid , mon_vid ) )

   CALL check( nf90_put_att(ncid, lat_vid , "long_name", "Latitude"        ) )
   CALL check( nf90_put_att(ncid, lat_vid , "units"    , "degrees_north"   ) )
   CALL check( nf90_put_att(ncid, lon_vid , "long_name", "Longitude"       ) )
   CALL check( nf90_put_att(ncid, lon_vid , "units"    , "degrees_east"    ) )
   CALL check( nf90_put_att(ncid, ns_vid  , "long_name", "solar band"      ) )
   CALL check( nf90_put_att(ncid, ns_vid  , "units"    , "1-dir,2-diff"    ) )
   CALL check( nf90_put_att(ncid, nr_vid  , "long_name", "radiation band"  ) )
   CALL check( nf90_put_att(ncid, nr_vid  , "units"    , "1-vis,2-nir"     ) )
   CALL check( nf90_put_att(ncid, lev_vid , "long_name", "urban layer"     ) )
   CALL check( nf90_put_att(ncid, lev_vid , "units"    , "urban layer"     ) )
   CALL check( nf90_put_att(ncid, den_vid , "long_name", "urban density"   ) )
   CALL check( nf90_put_att(ncid, den_vid , "units"    , "1-tbd,2-hd,3-md" ) )
   CALL check( nf90_put_att(ncid, mon_vid , "long_name", "month"           ) )
   CALL check( nf90_put_att(ncid, mon_vid , "units"    , "month"           ) )

   XY2D = (/ lon_dimid, lat_dimid /)
   XY3D = (/ lon_dimid, lat_dimid, den_dimid /)
   CALL check( nf90_def_var(ncid, "URBAN_TREE_PCT" , NF90_FLOAT, XY3D, pct_tcvid  ) )
   CALL check( nf90_def_var(ncid, "URBAN_WATER_PCT", NF90_FLOAT, XY3D, pct_urwtvid) )
   CALL check( nf90_def_var(ncid, "URBAN_TREE_TOP" , NF90_FLOAT, XY3D, htop_urvid ) )
   CALL check( nf90_def_var(ncid, "CANYON_HWR"     , NF90_FLOAT, XY3D, hwr_canvid ) )
   CALL check( nf90_def_var(ncid, "WTLUNIT_ROOF"   , NF90_FLOAT, XY3D, wt_rfvid   ) )
   CALL check( nf90_def_var(ncid, "WTROAD_PERV"    , NF90_FLOAT, XY3D, wt_rdvid   ) )
   CALL check( nf90_def_var(ncid, "EM_ROOF"        , NF90_FLOAT, XY3D, em_rfvid   ) )
   CALL check( nf90_def_var(ncid, "EM_WALL"        , NF90_FLOAT, XY3D, em_wlvid   ) )
   CALL check( nf90_def_var(ncid, "EM_IMPROAD"     , NF90_FLOAT, XY3D, em_imrdvid ) )
   CALL check( nf90_def_var(ncid, "EM_PERROAD"     , NF90_FLOAT, XY3D, em_perdvid ) )
   CALL check( nf90_def_var(ncid, "HT_ROOF"        , NF90_FLOAT, XY3D, ht_rfvid   ) )
   CALL check( nf90_def_var(ncid, "THICK_ROOF"     , NF90_FLOAT, XY3D, th_rfvid   ) )
   CALL check( nf90_def_var(ncid, "THICK_WALL"     , NF90_FLOAT, XY3D, th_wlvid   ) )
   CALL check( nf90_def_var(ncid, "T_BUILDING_MIN" , NF90_FLOAT, XY3D, tbminvid   ) )
   CALL check( nf90_def_var(ncid, "T_BUILDING_MAX" , NF90_FLOAT, XY3D, tbmaxvid   ) )

   XY4D = (/ lon_dimid, lat_dimid, den_dimid, ulev_dimid /)
   UR3D = (/ lon_dimid, lat_dimid, den_dimid  /)
   CALL check( nf90_def_var(ncid, "URBAN_PCT" , NF90_FLOAT, UR3D, pct_urvid ) )
   CALL check( nf90_def_var(ncid, "TK_ROOF"   , NF90_FLOAT, XY4D, tk_rfvid  ) )
   CALL check( nf90_def_var(ncid, "TK_WALL"   , NF90_FLOAT, XY4D, tk_wlvid  ) )
   CALL check( nf90_def_var(ncid, "TK_IMPROAD", NF90_FLOAT, XY4D, tk_imrdvid) )
   CALL check( nf90_def_var(ncid, "CV_ROOF"   , NF90_FLOAT, XY4D, cv_rfvid  ) )
   CALL check( nf90_def_var(ncid, "CV_WALL"   , NF90_FLOAT, XY4D, cv_wlvid  ) )
   CALL check( nf90_def_var(ncid, "CV_IMPROAD", NF90_FLOAT, XY4D, cv_imrdvid) )

   UL3D = (/ lon_dimid, lat_dimid, den_dimid, mon_dimid /)
   CALL check( nf90_def_var(ncid, "URBAN_TREE_LAI", NF90_FLOAT, UL3D, ur_laivid ) )
   CALL check( nf90_def_var(ncid, "URBAN_TREE_SAI", NF90_FLOAT, UL3D, ur_saivid ) )

   XY5D = (/ lon_dimid, lat_dimid, den_dimid, nr_dimid, ns_dimid /)
   CALL check( nf90_def_var(ncid, "ALB_ROOF"   , NF90_FLOAT, XY5D, alb_rfvid   ) )
   CALL check( nf90_def_var(ncid, "ALB_WALL"   , NF90_FLOAT, XY5D, alb_wlvid   ) )
   CALL check( nf90_def_var(ncid, "ALB_IMPROAD", NF90_FLOAT, XY5D, alb_imrdvid ) )
   CALL check( nf90_def_var(ncid, "ALB_PERROAD", NF90_FLOAT, XY5D, alb_perdvid ) )

   CALL check( nf90_put_att(ncid, pct_urvid , "units"     , "%"                                          ) )
   CALL check( nf90_put_att(ncid, pct_urvid , "long_name" , "Percentage of each urban type (density)"    ) )
   CALL check( nf90_put_att(ncid, pct_urvid , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, pct_tcvid, "units"     , "%"                                           ) )
   CALL check( nf90_put_att(ncid, pct_tcvid, "long_name" , "Urban percent tree cover"                    ) )
   CALL check( nf90_put_att(ncid, pct_tcvid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, htop_urvid, "units"     , "m"                                          ) )
   CALL check( nf90_put_att(ncid, htop_urvid, "long_name" , "Urban tree top height"                      ) )
   CALL check( nf90_put_att(ncid, htop_urvid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, pct_urwtvid, "units"     , "%"                                         ) )
   CALL check( nf90_put_att(ncid, pct_urwtvid, "long_name" , "Urban percent water cover"                 ) )
   CALL check( nf90_put_att(ncid, pct_urwtvid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, ur_laivid, "units"     , "m^2/m^2"                                     ) )
   CALL check( nf90_put_att(ncid, ur_laivid, "long_name" , "Urban tree monthly lai"                      ) )
   CALL check( nf90_put_att(ncid, ur_laivid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, ur_saivid, "units"     , "m^2/m^2"                                     ) )
   CALL check( nf90_put_att(ncid, ur_saivid, "long_name" , "Urban tree monthly sai"                      ) )
   CALL check( nf90_put_att(ncid, ur_saivid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, hwr_canvid, "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, hwr_canvid, "long_name" , "canyon height to width ratio"               ) )
   CALL check( nf90_put_att(ncid, hwr_canvid, "_FillValue", fv   ) )

   CALL check( nf90_put_att(ncid, wt_rfvid  , "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, wt_rfvid  , "long_name" , "fraction of roof"                           ) )
   CALL check( nf90_put_att(ncid, wt_rfvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, wt_rdvid  , "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, wt_rdvid  , "long_name" , "fraction of pervious road"                  ) )
   CALL check( nf90_put_att(ncid, wt_rdvid  , "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, em_rfvid  , "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, em_rfvid  , "long_name" , "emissivity of roof"                         ) )
   CALL check( nf90_put_att(ncid, em_rfvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, em_wlvid  , "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, em_wlvid  , "long_name" , "emissivity of wall"                         ) )
   CALL check( nf90_put_att(ncid, em_wlvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, em_imrdvid, "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, em_imrdvid, "long_name" , "emissivity of impervious road"              ) )
   CALL check( nf90_put_att(ncid, em_imrdvid, "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, em_perdvid, "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, em_perdvid, "long_name" , "emissivity of pervious road"                ) )
   CALL check( nf90_put_att(ncid, em_perdvid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, ht_rfvid  , "units"     , "meters"                                     ) )
   CALL check( nf90_put_att(ncid, ht_rfvid  , "long_name" , "height of roof"                             ) )
   CALL check( nf90_put_att(ncid, ht_rfvid  , "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, th_rfvid  , "units"     , "m"                                          ) )
   CALL check( nf90_put_att(ncid, th_rfvid  , "long_name" , "thickness of roof"                          ) )
   CALL check( nf90_put_att(ncid, th_rfvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, th_wlvid  , "units"     , "m"                                          ) )
   CALL check( nf90_put_att(ncid, th_wlvid  , "long_name" , "thickness of wall"                          ) )
   CALL check( nf90_put_att(ncid, th_wlvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, tbminvid  , "units"     , "K"                                          ) )
   CALL check( nf90_put_att(ncid, tbminvid  , "long_name" , "minimum interior building temperature"      ) )
   CALL check( nf90_put_att(ncid, tbminvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, tbmaxvid  , "units"     , "K"                                          ) )
   CALL check( nf90_put_att(ncid, tbmaxvid  , "long_name" , "maximum interior building temperature"      ) )
   CALL check( nf90_put_att(ncid, tbmaxvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, tk_rfvid  , "units"     , "W/m*K"                                      ) )
   CALL check( nf90_put_att(ncid, tk_rfvid  , "long_name" , "thermal conductivity of roof"               ) )
   CALL check( nf90_put_att(ncid, tk_rfvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, tk_wlvid  , "units", "W/m*K"                                           ) )
   CALL check( nf90_put_att(ncid, tk_wlvid  , "long_name", "thermal conductivity of wall"                ) )
   CALL check( nf90_put_att(ncid, tk_wlvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, tk_imrdvid, "units"     , "W/m*K"                                      ) )
   CALL check( nf90_put_att(ncid, tk_imrdvid, "long_name" , "thermal conductivity of impervious road"    ) )
   CALL check( nf90_put_att(ncid, tk_imrdvid, "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, cv_rfvid  , "units"     , "J/m^3*K"                                    ) )
   CALL check( nf90_put_att(ncid, cv_rfvid  , "long_name" , "volumetric heat capacity of roof"           ) )
   CALL check( nf90_put_att(ncid, cv_rfvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, cv_wlvid  , "units"     , "J/m^3*K"                                    ) )
   CALL check( nf90_put_att(ncid, cv_wlvid  , "long_name" , "volumetric heat capacity of wall"           ) )
   CALL check( nf90_put_att(ncid, cv_wlvid  , "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, cv_imrdvid, "units"     , "J/m^3*K"                                    ) )
   CALL check( nf90_put_att(ncid, cv_imrdvid, "long_name" , "volumetric heat capacity of impervious road") )
   CALL check( nf90_put_att(ncid, cv_imrdvid, "_FillValue", fv) )

   CALL check( nf90_put_att(ncid, alb_rfvid , "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, alb_rfvid , "long_name" , "albedo of roof"                             ) )
   CALL check( nf90_put_att(ncid, alb_rfvid , "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, alb_wlvid , "units"     , "unitless"                                   ) )
   CALL check( nf90_put_att(ncid, alb_wlvid , "long_name" , "albedo of wall"                             ) )
   CALL check( nf90_put_att(ncid, alb_wlvid , "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, alb_imrdvid, "units"     , "unitless"                                  ) )
   CALL check( nf90_put_att(ncid, alb_imrdvid, "long_name" , "albedo of impervious road"                 ) )
   CALL check( nf90_put_att(ncid, alb_imrdvid, "_FillValue", fv) )
 
   CALL check( nf90_put_att(ncid, alb_perdvid, "units"     , "unitless"                                  ) )
   CALL check( nf90_put_att(ncid, alb_perdvid, "long_name" , "albedo of pervious road"                   ) )
   CALL check( nf90_put_att(ncid, alb_perdvid, "_FillValue", fv) )

   CALL check( nf90_enddef(ncid) )
   
   CALL check( nf90_inq_varid(ncid, "lat"            , urlat_vid  ) )
   CALL check( nf90_put_var  (ncid, urlat_vid        , latso     ) )

   CALL check( nf90_inq_varid(ncid, "lon"            , urlon_vid  ) )
   CALL check( nf90_put_var  (ncid, urlon_vid        , lonso      ) )

   call check( nf90_inq_varid(ncid, "numsolar"       , ns_vid     ) )
   CALL check( nf90_put_var  (ncid, ns_vid           , n_ns       ) )

   CALL check( nf90_inq_varid(ncid, "numrad"         , nr_vid     ) )
   CALL check( nf90_put_var  (ncid, nr_vid           , n_nr       ) )

   CALL check( nf90_inq_varid(ncid, "month"          , mon_vid    ) )
   CALL check( nf90_put_var  (ncid, mon_vid          , n_mon      ) )

   CALL check( nf90_inq_varid(ncid, "density"        , den_vid    ) )
   CALL check( nf90_put_var  (ncid, den_vid          , n_den      ) )

   CALL check( nf90_inq_varid(ncid, "ulev"           , lev_vid    ) )
   CALL check( nf90_put_var  (ncid, lev_vid          , n_ulev     ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_PCT" , pct_tcvid  ) )
   CALL check( nf90_put_var  (ncid, pct_tcvid        , pct_tc    ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_WATER_PCT", pct_urwtvid) )
   CALL check( nf90_put_var  (ncid, pct_urwtvid      , pct_urwt  ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_TOP" , htop_urvid ) )
   CALL check( nf90_put_var  (ncid, htop_urvid       , htop_ur   ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_LAI" , ur_laivid  ) )
   CALL check( nf90_put_var  (ncid, ur_laivid        , ur_lai    ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_SAI" , ur_saivid  ) )
   CALL check( nf90_put_var  (ncid, ur_saivid        , ur_sai    ) )

   CALL check( nf90_inq_varid(ncid, "CANYON_HWR"     , hwr_canvid ) )
   CALL check( nf90_put_var  (ncid,hwr_canvid        , hwr_can   ) )

   CALL check( nf90_inq_varid(ncid, "WTLUNIT_ROOF"   , wt_rfvid   ) )
   CALL check( nf90_put_var  (ncid, wt_rfvid         , wt_rf     ) )

   CALL check( nf90_inq_varid(ncid, "WTROAD_PERV"    , wt_rdvid   ) )
   CALL check( nf90_put_var  (ncid, wt_rdvid         , wt_rd     ) )

   CALL check( nf90_inq_varid(ncid, "EM_ROOF"        , em_rfvid   ) )
   CALL check( nf90_put_var  (ncid, em_rfvid         , em_rf     ) )

   CALL check( nf90_inq_varid(ncid, "EM_WALL"        , em_wlvid   ) )
   CALL check( nf90_put_var  (ncid, em_wlvid         , em_wl     ) )

   CALL check( nf90_inq_varid(ncid, "EM_IMPROAD"     , em_imrdvid ) )
   CALL check( nf90_put_var  (ncid, em_imrdvid       , em_imrd   ) )

   CALL check( nf90_inq_varid(ncid, "EM_PERROAD"     , em_perdvid ) )
   CALL check( nf90_put_var  (ncid, em_perdvid       , em_perd   ) )

   CALL check( nf90_inq_varid(ncid, "HT_ROOF"        , ht_rfvid   ) )
   CALL check( nf90_put_var  (ncid, ht_rfvid         , ht_rf     ) )

   !CALL check( nf90_inq_varid(ncid, "WIND_HGT_CANYON",whcvid   ) )
   !CALL check( nf90_put_var(ncid,whcvid,w_hc                   ) )

   CALL check( nf90_inq_varid(ncid, "THICK_ROOF"     , th_rfvid   ) )
   CALL check( nf90_put_var  (ncid, th_rfvid         , th_rf     ) )

   CALL check( nf90_inq_varid(ncid, "THICK_WALL"     , th_wlvid   ) )
   CALL check( nf90_put_var  (ncid, th_wlvid         , th_wl     ) )

   CALL check( nf90_inq_varid(ncid, "T_BUILDING_MIN" , tbminvid   ) )
   CALL check( nf90_put_var  (ncid, tbminvid         , tb_min    ) )

   CALL check( nf90_inq_varid(ncid, "T_BUILDING_MAX" , tbmaxvid   ) )
   CALL check( nf90_put_var  (ncid, tbmaxvid         , tb_max    ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_PCT"      , pct_urvid  ) )
   CALL check( nf90_put_var  (ncid, pct_urvid        , pct_ur    ) )

   CALL check( nf90_inq_varid(ncid, "TK_ROOF"        , tk_rfvid   ) )
   CALL check( nf90_put_var  (ncid, tk_rfvid         , tk_rf     ) )

   CALL check( nf90_inq_varid(ncid, "TK_WALL"        , tk_wlvid   ) )
   CALL check( nf90_put_var  (ncid, tk_wlvid         , tk_wl     ) )

   CALL check( nf90_inq_varid(ncid, "TK_IMPROAD"     , tk_imrdvid ) )
   CALL check( nf90_put_var  (ncid, tk_imrdvid       , tk_imrd   ) )

   CALL check( nf90_inq_varid(ncid, "CV_ROOF"        , cv_rfvid   ) )
   CALL check( nf90_put_var  (ncid, cv_rfvid         , cv_rf     ) )

   CALL check( nf90_inq_varid(ncid, "CV_WALL"        , cv_wlvid   ) )
   CALL check( nf90_put_var  (ncid, cv_wlvid         , cv_wl     ) )

   CALL check( nf90_inq_varid(ncid, "CV_IMPROAD"     , cv_imrdvid ) )
   CALL check( nf90_put_var  (ncid, cv_imrdvid       , cv_imrd   ) )

   CALL check( nf90_inq_varid(ncid, "ALB_ROOF"       , alb_rfvid  ) )
   CALL check( nf90_put_var  (ncid, alb_rfvid        , alb_rf    ) )

   CALL check( nf90_inq_varid(ncid, "ALB_WALL"       , alb_wlvid  ) )
   CALL check( nf90_put_var  (ncid, alb_wlvid        , alb_wl    ) )

   CALL check( nf90_inq_varid(ncid, "ALB_IMPROAD"    , alb_imrdvid) )
   CALL check( nf90_put_var  (ncid, alb_imrdvid      , alb_imrd  ) )

   CALL check( nf90_inq_varid(ncid, "ALB_PERROAD"    , alb_perdvid) )
   CALL check( nf90_put_var  (ncid, alb_perdvid      , alb_perd  ) )

   CALL check( nf90_close(ncid) )

   PRINT*, "*** SUCCESS write surface file ***"
   

   DEALLOCATE( hlat    )
   DEALLOCATE( hlats   )
   DEALLOCATE( hlatn   )
   DEALLOCATE( hlon    )
   DEALLOCATE( hlonw   )
   DEALLOCATE( hlone   )
   DEALLOCATE( harea   )
   DEALLOCATE( gfcc_tc )
   DEALLOCATE( gedi_th )
   DEALLOCATE( gl30_wt )

   DEALLOCATE( urrgid  )

   DEALLOCATE( latso     )
   DEALLOCATE( lonso     )
   DEALLOCATE( area      )
   DEALLOCATE( pct_urban )
   DEALLOCATE( htop_pft  )
   DEALLOCATE( pct_pft   )
   DEALLOCATE( lai_pft   )
   DEALLOCATE( ur_dc     )
   DEALLOCATE( pct_ur    )

   DEALLOCATE( pct_tc    )
   DEALLOCATE( pct_urwt  )
   DEALLOCATE( htop_ur   )
   DEALLOCATE( ur_clss   )
   DEALLOCATE( ur_rgid   )
   DEALLOCATE( ur_lai    )
   DEALLOCATE( hwr_can   )
   DEALLOCATE( wt_rf     )
   DEALLOCATE( wt_rd     )
   DEALLOCATE( em_rf     )
   DEALLOCATE( em_wl     )
   DEALLOCATE( em_imrd   )
   DEALLOCATE( em_perd   )
   DEALLOCATE( ht_rf     )
   DEALLOCATE( w_hc      )
   DEALLOCATE( ulev_imrd )
   DEALLOCATE( th_rf     )
   DEALLOCATE( th_wl     )
   DEALLOCATE( tb_min    )
   DEALLOCATE( tb_max    )
   DEALLOCATE( cv_rf     )
   DEALLOCATE( cv_wl     )
   DEALLOCATE( cv_imrd   )
   DEALLOCATE( tk_rf     )
   DEALLOCATE( tk_wl     )
   DEALLOCATE( tk_imrd   )
   DEALLOCATE( alb_rf    )
   DEALLOCATE( alb_wl    )
   DEALLOCATE( alb_imrd  )
   DEALLOCATE( alb_perd  )

CONTAINS 
   
   SUBROUTINE check(status)
      INTEGER, intent(in) :: status

      IF (status /= nf90_noerr) THEN
         PRINT *, trim( nf90_strerror(status))
         stop 2
      ENDIF 
   END SUBROUTINE check

END PROGRAM clmu2grid
