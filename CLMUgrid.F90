PROGRAM clmu2grid

   use netcdf

   IMPLICIT NONE

   INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12)

   REAL(r8), PARAMETER :: sdelta = 1._r8/240._r8
   REAL(r8), PARAMETER :: udelta = 1._r8/120._r8
   REAL    , PARAMETER :: fv     = -999

   INTEGER, PARAMETER :: plat = 600 , plon = 600
   INTEGER, PARAMETER :: ulat = 600 , ulon = 600
   INTEGER, PARAMETER :: nlat = 1200, nlon = 1200
   INTEGER, PARAMETER :: rid = 33   , den_clss = 3
   INTEGER, PARAMETER :: nxo = 720  , nyo = 360
   INTEGER, PARAMETER :: ns  = 2    , nr  = 2
   INTEGER, PARAMETER :: ulev= 10
   INTEGER, PARAMETER :: mon = 12
   INTEGER, PARAMETER :: npft= 16 
   
! input variables
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: tlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: tlon
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlats
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlatn
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlon
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlonw
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlone
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: nrlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlats
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlatn
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlon
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlonw
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: urlone
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gfcc_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gedi_th, ngedi_th
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: modur, nmodur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gl30_wt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: mask
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: harea
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: uarea
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: larea

   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urclss
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urrgid
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urden

   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: umd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: uhd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: utbd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: umask

   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: npct_urban
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: pct_urban
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: ugt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: uvg
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: pct_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: npct_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: htop_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: nhtop_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: lai_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: nlai_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: sai_pft
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: nsai_pft

   REAL(r8), DIMENSION(3,rid)      :: hwrcan, wtrf, wtrd, emrf, emwl
   REAL(r8), DIMENSION(3,rid)      :: emimrd, emperd, htrf, whc, ulevimrd 
   REAL(r8), DIMENSION(3,rid)      :: thrf, thwl, tbmin, tbmax
   
   REAL(r8), DIMENSION(3,rid,10)   :: cvrf, cvwl, cvimrd, &
                                             tkrf, tkwl, tkimrd
   REAL(r8), DIMENSION(3,rid,nr,ns):: albrf, albwl, albimrd, albperd
   ! output variables
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: ur_clss
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: ur_rgid
   
   REAL(r8), ALLOCATABLE, DIMENSION(:)     :: latso, flatso
   REAL(r8), ALLOCATABLE, DIMENSION(:)     :: lonso
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wgt_top
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: cnt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: area
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: lurrgid
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: nmd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: ntbd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: pct_tc, fpct_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: pct_urwt, fpct_urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop_ur, fhtop_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: hwr_can, fhwr_can
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wt_rf, fwt_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wt_rd, fwt_rd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_rf, fem_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_wl, fem_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_imrd, fem_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_perd, fem_perd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: ht_rf, fht_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: w_hc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: ulev_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: th_rf, fth_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: th_wl, fth_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tb_min, ftb_min
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tb_max, ftb_max

   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: pct_tland
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: pct_wland
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: pct_hland
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: ur_land
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: hgt, avg
   
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: ur_dc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: pct_ur, fpct_ur
  ! REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: lai_wgt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ur_lai, fur_lai
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ur_sai, fur_sai
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_rf, fcv_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_wl, fcv_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_imrd, fcv_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_rf, ftk_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_wl, ftk_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_imrd, ftk_imrd

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_rf, falb_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_wl, falb_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_imrd, falb_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_perd, falb_perd

   ! variable ids
   INTEGER :: ncid, uhd_vid, umd_vid, utbd_vid, htop_urvid, reg
   INTEGER :: lat_vid, lon_vid, lat_dimid, lon_dimid, ns_dimid, nr_dimid, ulev_dimid
   INTEGER :: ns_vid, nr_vid, pct_tcvid, pct_urvid, pct_urwtvid, ulev_vid, lev_vid
   INTEGER :: pftvid, laivid, maskid, umaskid, t_pftvid, ur_laivid, htop_pftvid
   INTEGER :: hlat_vid, hlon_vid, gfcc_tcvid, gedi_thvid, gl30_wtvid
   INTEGER :: urlat_vid, urlon_vid, ur_clssvid, ur_rgvid, hwr_canvid
   INTEGER :: wt_rfvid, wt_rdvid, em_rfvid, em_wlvid, em_imrdvid, em_perdvid
   INTEGER :: ht_rfvid, whcvid, cv_rfvid, cv_wlvid, cv_imrdvid, ulev_imrdvid
   INTEGER :: th_rfvid, th_wlvid, tbminvid, tbmaxvid
   INTEGER :: tk_rfvid, tk_wlvid, tk_imrdvid
   INTEGER :: alb_rfvid, alb_imrdvid, alb_perdvid, alb_wlvid
   INTEGER :: uxid, uyid, upftvid, mon_dimid, den_dimid
   INTEGER :: den_vid, mon_vid, ur_landvid, saivid, ur_saivid, ur_denvid

   CHARACTER(len=255) :: SRF_DIR='/hard/dongwz/github/CLMUgrid/v4/srf_5x5/'
   CHARACTER(len=255) :: RAW_DIR='/hard/yuanhua/mksrf/srf_5x5/'
   CHARACTER(len=255) :: REGFILE='5x5_F90'
   CHARACTER(len=255) :: filename
   CHARACTER(len=4)    , DIMENSION(4) :: creg

   REAL(r8) :: fac(3)
   REAL(r8) :: pi, deg2rad, re, dx, dy, sumarea, sumur, sumpct
   REAL(r8) :: dll, x_delta, y_delta, wgt
   !REAL(r8) :: lone1(nxy1), lonw1(nxy1), latn1(nxy1), lats1(nxy1)
   REAL(r8) :: lone2(nxo), lonw2(nxo), latn2(nyo), lats2(nyo)
   REAL(r8) :: lone3(nxo), lonw3(nxo), latn3(nyo), lats3(nyo)

   INTEGER  :: i, j, k, io, jo, m, n, ii, jj, ir, jr, cont, sumth, p

   !REAL(r8) :: wgt(3)
   INTEGER  :: n_ns(2), n_nr(2), n_den(3), n_ulev(10), n_mon(12)
   INTEGER  :: XY2D(2), XY3D(3), XY4D(4), UR3D(3), UL3D(4), XY5D(5)

   pi = 4._r8*atan(1.)
   deg2rad = pi/180._r8
   re = 6.37122e6 * 0.001

   allocate( tlat  (nlat) )
   allocate( tlon  (nlon) )
   allocate( hlat  (nlat) )
   allocate( hlats (nlat) )
   allocate( hlatn (nlat) )
   allocate( hlon  (nlon) )
   allocate( hlonw (nlon) )
   allocate( hlone (nlon) )
   allocate( urlat (plat) )
   allocate( nrlat (plat) )
   allocate( urlats(plat) )
   allocate( urlatn(plat) )
   allocate( urlon (plon) )
   allocate( urlonw(plon) )
   allocate( urlone(plon) )

   allocate( harea   (nlon, nlat) )
   allocate( mask    (nlon, nlat) )
   allocate( gfcc_tc (nlon, nlat) )
   allocate( gedi_th (nlon, nlat) )
   allocate( ngedi_th(nlon, nlat) )
   allocate( gl30_wt (nlon, nlat) )
   allocate( modur   (nlon, nlat) )
   allocate( nmodur  (nlon, nlat) )

   allocate( uarea  (plon, plat) )
   allocate( umask  (plon, plat) )
   allocate( umd    (ulon, ulat) )
   allocate( nmd    (nlon, nlat) )
   allocate( uhd    (ulon, ulat) )
   allocate( lurrgid(nxo , nyo ) )
   allocate( utbd   (ulon, ulat) )
   allocate( ntbd   (ulon, ulat) )
   allocate( urclss (plon, plat) )
   allocate( urrgid (plon, plat) )
   allocate( urden  (plon, plat) )

   allocate( latso     (nyo) )
   allocate( lonso     (nxo) )

   allocate( area      (nxo, nyo) )
   allocate( hgt       (nxo, nyo) )
   allocate( avg       (nxo, nyo) )
   allocate( pct_urban (nxo, nyo) )
   allocate( npct_urban(nxo, nyo) )

   allocate( htop_pft  (nxo, nyo, npft     ) )
   allocate( nhtop_pft (nxo, nyo, npft     ) )
   allocate( pct_pft   (nxo, nyo, npft     ) )
   allocate( npct_pft  (nxo, nyo, npft     ) )
   allocate( lai_pft   (nxo, nyo, npft, mon) )
   allocate( nlai_pft  (nxo, nyo, npft, mon) )
   allocate( sai_pft   (nxo, nyo, npft, mon) )
   allocate( nsai_pft  (nxo, nyo, npft, mon) )
   allocate( ur_dc     (nxo, nyo, den_clss ) )
   allocate( pct_ur    (nxo, nyo, den_clss ) )

   allocate( pct_tland(nxo, nyo) )
   allocate( pct_wland(nxo, nyo) )
   allocate( pct_hland(nxo, nyo) )
   allocate( ur_land  (nxo, nyo) )
   allocate( ur_clss  (nxo, nyo) )
   allocate( ur_rgid  (nxo, nyo) )

   allocate( ugt      (nxo, nyo, 3) )
   allocate( uvg      (nxo, nyo, 3) )
   allocate( wgt_top  (nxo, nyo, 3) )
   allocate( cnt      (nxo, nyo, 3) )

   allocate( tc       (nxo, nyo, den_clss) )
   allocate( urwt     (nxo, nyo, den_clss) )
   allocate( htop     (nxo, nyo, den_clss) )
   allocate( pct_tc   (nxo, nyo, den_clss) )
   allocate( pct_urwt (nxo, nyo, den_clss) )
   allocate( htop_ur  (nxo, nyo, den_clss) )
   allocate( hwr_can  (nxo, nyo, den_clss) )
   allocate( wt_rf    (nxo, nyo, den_clss) )
   allocate( wt_rd    (nxo, nyo, den_clss) )
   allocate( em_rf    (nxo, nyo, den_clss) )
   allocate( em_wl    (nxo, nyo, den_clss) )
   allocate( em_imrd  (nxo, nyo, den_clss) )
   allocate( em_perd  (nxo, nyo, den_clss) )
   allocate( ht_rf    (nxo, nyo, den_clss) )
   allocate( w_hc     (nxo, nyo, den_clss) )
   allocate( ulev_imrd(nxo, nyo, den_clss) )
   allocate( th_rf    (nxo, nyo, den_clss) )
   allocate( th_wl    (nxo, nyo, den_clss) )
   allocate( tb_min   (nxo, nyo, den_clss) )
   allocate( tb_max   (nxo, nyo, den_clss) )
   allocate( ur_sai   (nxo, nyo, den_clss, mon ) )
   allocate( ur_lai   (nxo, nyo, den_clss, mon ) )
   allocate( cv_rf    (nxo, nyo, den_clss, ulev) )
   allocate( cv_wl    (nxo, nyo, den_clss, ulev) )
   allocate( cv_imrd  (nxo, nyo, den_clss, ulev) )
   allocate( tk_rf    (nxo, nyo, den_clss, ulev) )
   allocate( tk_wl    (nxo, nyo, den_clss, ulev) )
   allocate( tk_imrd  (nxo, nyo, den_clss, ulev) )
   allocate( alb_rf   (nxo, nyo, den_clss, nr, ns) )
   allocate( alb_wl   (nxo, nyo, den_clss, nr, ns) )
   allocate( alb_imrd (nxo, nyo, den_clss, nr, ns) )
   allocate( alb_perd (nxo, nyo, den_clss, nr, ns) )

   allocate( flatso    (nyo) )
   allocate( fpct_tc   (nxo, nyo, den_clss) )
   allocate( fpct_urwt (nxo, nyo, den_clss) )
   allocate( fhtop_ur  (nxo, nyo, den_clss) )
   allocate( fpct_ur   (nxo, nyo, den_clss) )
   allocate( fhwr_can  (nxo, nyo, den_clss) )
   allocate( fwt_rf    (nxo, nyo, den_clss) )
   allocate( fwt_rd    (nxo, nyo, den_clss) )
   allocate( fem_rf    (nxo, nyo, den_clss) )
   allocate( fem_wl    (nxo, nyo, den_clss) )
   allocate( fem_imrd  (nxo, nyo, den_clss) )
   allocate( fem_perd  (nxo, nyo, den_clss) )
   allocate( fht_rf    (nxo, nyo, den_clss) )
   allocate( fth_rf    (nxo, nyo, den_clss) )
   allocate( fth_wl    (nxo, nyo, den_clss) )
   allocate( ftb_min   (nxo, nyo, den_clss) )
   allocate( ftb_max   (nxo, nyo, den_clss) )
   allocate( fur_sai   (nxo, nyo, den_clss, mon ) )
   allocate( fur_lai   (nxo, nyo, den_clss, mon ) )
   allocate( fcv_rf    (nxo, nyo, den_clss, ulev) )
   allocate( fcv_wl    (nxo, nyo, den_clss, ulev) )
   allocate( fcv_imrd  (nxo, nyo, den_clss, ulev) )
   allocate( ftk_rf    (nxo, nyo, den_clss, ulev) )
   allocate( ftk_wl    (nxo, nyo, den_clss, ulev) )
   allocate( ftk_imrd  (nxo, nyo, den_clss, ulev) )
   allocate( falb_rf   (nxo, nyo, den_clss, nr, ns) )
   allocate( falb_wl   (nxo, nyo, den_clss, nr, ns) )
   allocate( falb_imrd (nxo, nyo, den_clss, nr, ns) )
   allocate( falb_perd (nxo, nyo, den_clss, nr, ns) )

   ! initialization
   fac      (:)     = 0.
   cnt      (:,:,:) = 0.
   ur_land  (:,:)   = 0.
   pct_tland(:,:)   = 0.
   pct_wland(:,:)   = 0.
   pct_hland(:,:)   = 0.
   utbd     (:,:)   = 0.
   uhd      (:,:)   = 0.
   umd      (:,:)   = 0.
   hgt      (:,:)   = 0.
   avg      (:,:)   = 0.
   tc       (:,:,:) = 0.
   urwt     (:,:,:) = 0.
   htop     (:,:,:) = 0.
   pct_tc   (:,:,:) = 0.
   pct_urwt (:,:,:) = 0.
   htop_ur  (:,:,:) = 0.
   ur_dc    (:,:,:) = 0.

   ugt      (:,:,:) = 0.
   uvg      (:,:,:) = 0.
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

   fpct_tc   (:,:,:) = 0.
   fpct_ur   (:,:,:) = 0.
   fhtop_ur  (:,:,:) = 0.
   fpct_urwt (:,:,:) = 0.
   fhwr_can  (:,:,:) = 0.
   fwt_rf    (:,:,:) = 0.
   fwt_rd    (:,:,:) = 0.
   fem_rf    (:,:,:) = 0.
   fem_wl    (:,:,:) = 0.
   fem_imrd  (:,:,:) = 0.
   fem_perd  (:,:,:) = 0.
   fht_rf    (:,:,:) = 0.
   fth_rf    (:,:,:) = 0.
   fth_wl    (:,:,:) = 0.
   ftb_min   (:,:,:) = 0.
   ftb_max   (:,:,:) = 0.

   fur_lai   (:,:,:,:) = 0.
   fur_sai   (:,:,:,:) = 0.
   ftk_rf    (:,:,:,:) = 0.
   ftk_wl    (:,:,:,:) = 0.
   ftk_imrd  (:,:,:,:) = 0.
   fcv_rf    (:,:,:,:) = 0.
   fcv_wl    (:,:,:,:) = 0.
   fcv_imrd  (:,:,:,:) = 0.

   falb_rf   (:,:,:,:,:) = 0.
   falb_wl   (:,:,:,:,:) = 0.
   falb_imrd (:,:,:,:,:) = 0.
   falb_perd (:,:,:,:,:) = 0.

   fpct_tc   (:,:,:) = 0.
   fpct_urwt (:,:,:) = 0.
   fhtop_ur  (:,:,:) = 0.

   ur_rgid(:,:) = 0.
   ur_clss(:,:) = 0.
   
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

   x_delta = (360._r8/nxo)*1._r8
   y_delta = (180._r8/nyo)*1._r8

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

   CALL check( nf90_open("urban_properties_data.5deg.211121-164510.nc", nf90_nowrite, ncid) )

   CALL check( nf90_inq_varid(ncid, "REGION_ID", ur_rgvid) )
   CALL check( nf90_get_var(ncid  , ur_rgvid   , lurrgid ) )

   CALL check( nf90_close(ncid) )

   print *, "    reading MODIS PFTs data"

   CALL check( nf90_open("global_0.5x0.5.MOD2005_V4.5.nc", nf90_nowrite, ncid) )

   CALL check( nf90_inq_varid(ncid, "PCT_PFT"        , pftvid     ) )
   CALL check( nf90_inq_varid(ncid, "MONTHLY_PFT_LAI", laivid     ) )
   CALL check( nf90_inq_varid(ncid, "MONTHLY_PFT_SAI", saivid     ) )
   ! CALL check( nf90_inq_varid(ncid, "PCT_URBAN"      , upftvid    ) )
   CALL check( nf90_inq_varid(ncid, "HTOP_PFT"       , htop_pftvid) )

   CALL check( nf90_get_var(ncid, pftvid     , npct_pft  ) )
   ! CALL check( nf90_get_var(ncid, upftvid    , npct_urban) )
   CALL check( nf90_get_var(ncid, saivid     , nsai_pft  ) )
   CALL check( nf90_get_var(ncid, laivid     , nlai_pft  ) )
   CALL check( nf90_get_var(ncid, htop_pftvid, nhtop_pft ) )

   CALL check( nf90_close(ncid) )

   ! data flip
   ! 调整PFTs数据与1km数据一致 pft(lon,n,npft) = pft(lon,lat-n+1,npft) lat(0)=90 -> lat(0)=-90
   DO i = 1, nyo
      DO j = 1, nxo
         k = nyo - i + 1
         ! pct_urban(j,i)     = npct_urban(j,k)
         htop_pft (j,i,:)   = nhtop_pft (j,k,:)
         pct_pft  (j,i,:)   = npct_pft  (j,k,:)
         lai_pft  (j,i,:,:) = nlai_pft  (j,k,:,:)
         sai_pft  (j,i,:,:) = nsai_pft  (j,k,:,:)
      ENDDO
   ENDDO

   DO i = 1, nyo
      lats2(i) = -90. + (i-1)*y_delta
      latn2(i) = -90. + i*y_delta
   ENDDO

   DO i = 1, nxo
      lonw2(i) = -180. + (i-1)*x_delta
      lone2(i) = -180. + i*x_delta
   ENDDO

   DO i = 1, nxo
      lonso(i) = -180. + i*x_delta - 0.5*x_delta
   ENDDO

   DO i = 1, nyo
      latso(i) =  -90. + i*y_delta - 0.5*y_delta
   ENDDO

   DO i = 1, nyo
      dx = (lone2(1)-lonw2(1))*deg2rad
      dy = sin(latn2(i)*deg2rad) - sin(lats2(i)*deg2rad)
      area(:,i) = dx*dy*re*re
   ENDDO

   OPEN(11,FILE=REGFILE)
   OPEN(12,FILE=REGFILE)
   
   DO WHILE(.TRUE.)
     ! process global 500m raw data
     print *, "*** processing 500m raw data ***"
     print *, "    reading data"
     READ(11,*,END=100) reg
     READ(12,*,END=101) creg
     filename = TRIM(SRF_DIR)//'RG_'//TRIM(creg(1))//'_'//&
                TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.SRF.nc'

     print*, filename
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

     ! filename = TRIM(ROOT_DIR)//TRIM(RAW_DIR)//'RG_'//TRIM(creg(1))//'_'//&
     !            TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.RAW2005.nc'
     ! CALL check( nf90_open(filename, nf90_nowrite, ncid) )

     ! CALL check( nf90_inq_varid(ncid, "tree_height", gedi_thvid) )
     ! CALL check( nf90_get_var  (ncid, gedi_thvid   , gedi_th   ) )

     ! CALL check( nf90_close(ncid) )

     print *, "    reading 1km urban class data"

     !print *, "    reading 1km urban properties data"
     filename = TRIM(SRF_DIR)//'RG_'//TRIM(creg(1))//'_'//&
                TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.NCAR.nc'
     CALL check( nf90_open(filename, nf90_nowrite, ncid) )

     CALL check( nf90_inq_varid(ncid, "lat"                , urlat_vid ) )
     CALL check( nf90_inq_varid(ncid, "lon"                , urlon_vid ) )
    !CALL check( nf90_inq_varid(ncid, "LANDMASK"           , umaskid     ) )
     CALL check( nf90_inq_varid(ncid, "URBAN_DENSITY_CLASS", ur_clssvid) )
     CALL check( nf90_inq_varid(ncid, "REGION_ID"          , ur_rgvid  ) )
     !CALL check( nf90_inq_varid(ncid, "URBAN_DENSITY"      , ur_denvid   ) 

     CALL check( nf90_get_var(ncid, urlat_vid   , urlat ) )
     CALL check( nf90_get_var(ncid, urlon_vid   , urlon ) )
     !CALL check( nf90_get_var(ncid, umaskid     , umask   ) )
     CALL check( nf90_get_var(ncid, ur_clssvid  , urden ) )
     CALL check( nf90_get_var(ncid, ur_rgvid    , urrgid) )

     CALL check( nf90_close(ncid) )

     print*, "    reading MOD URBAN data"
     filename = TRIM(RAW_DIR)//'RG_'//TRIM(creg(1))//'_'//&
                TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.MOD2005.nc'
     print*, filename

     CALL check( nf90_open(filename, nf90_nowrite, ncid      ) )
     CALL check( nf90_inq_varid(ncid, "PCT_URBAN", upftvid   ) )
     CALL check( nf90_inq_varid(ncid, "HTOP"     , gedi_thvid) )

     CALL check( nf90_get_var  (ncid, upftvid    , nmodur    ) )
     CALL check( nf90_get_var  (ncid, gedi_thvid , ngedi_th  ) )
     CALL check( nf90_close(ncid) )

     DO i = 1, 1200
        k = 1200 - i + 1
        ! pct_urban(j,i)     = npct_urban(j,k)
        modur  (:,i) = nmodur  (:,k)
        gedi_th(:,i) = ngedi_th(:,k)
     ENDDO

     ! calculate output lat/lon
     

     ! calculate the edge of small grids
     DO i = 1, nlat
        hlats(i) = hlat(1) - sdelta/2. + (i-1)*sdelta
        hlatn(i) = hlat(1) - sdelta/2. + i*sdelta
     ENDDO
   
     DO i = 1, nlon
        hlonw(i) = hlon(1) - sdelta/2. + (i-1)*sdelta
        hlone(i) = hlon(1) - sdelta/2. + i*sdelta
     ENDDO
   
     DO i = 1, nlat
        dx = (hlone(1)-hlonw(1))*deg2rad
        dy = sin(hlatn(i)*deg2rad) - sin(hlats(i)*deg2rad)
        harea(:,i) = dx*dy*re*re
     ENDDO

      DO i = 1, plat
         urlatn(i) = urlat(1) - udelta/2. + i*udelta
         urlats(i) = urlat(1) - udelta/2. + (i-1)*udelta
         ! urlat(i)  = -90. + i*udelta - 0.5*udelta
      ENDDO

      DO i = 1, plon
         urlonw(i) = urlon(1) - udelta/2. + (i-1)*udelta
         urlone(i) = urlon(1) - udelta/2. + i*udelta
         ! urlon(i)  = -180. + i*udelta - 0.5*udelta
      ENDDO

      DO i = 1, plat
         dx = (urlone(1)-urlonw(1))*deg2rad
         dy = sin(urlatn(i)*deg2rad) - sin(urlats(i)*deg2rad)
         uarea(:,i) = dx*dy*re*re
      ENDDO

      print *, "    aggregating data"
      print *, "    aggregating urban class data"
      DO i = 1, 1200
         DO j = 1, 1200
            ! calculate io, jo
            ! ir = NINT((hlats(i)+sdelta/2+ 90.)/udelta+0.5)
            ! jr = NINT((hlonw(j)+sdelta/2+180.)/udelta+0.5)
            ir = CEILING(i*1._r8/2)
            jr = CEILING(j*1._r8/2)
            io = NINT((hlats(i)+sdelta/2+ 90.)/y_delta+0.5)
            jo = NINT((hlonw(j)+sdelta/2+180.)/x_delta+0.5)

            !聚合Simard树高
            !加权：
            !粗网格树高=粗网格树高+simard树高*1km格点面积
            !加权系数：粗网格格点面积
            IF (gedi_th(j,i) > 0) THEN
               hgt(jo,io) = hgt(jo,io) + gedi_th(j,i)*harea(j,i)
               avg(jo,io) = avg(jo,io) + harea(j,i)
            ENDIF

            IF (urden(jr,ir) > 0) THEN
             !Tall-Building-Distinc urban
               IF (urden(jr,ir) == 1) THEN
            
               !加权：
               !粗网格城市水体(植被)覆盖度=粗网格城市水体(植被)覆盖度+500m城市格点水体(植被)覆盖度*500m城市格点面积
               !加权系数；粗网格城市格点面积
                  IF (gl30_wt(j,i) > 0.) THEN
                     urwt (jo,io,1) = urwt (jo,io,1) + gl30_wt(j,i)*harea(j,i)
                  ENDIF
                  IF (gfcc_tc(j,i) >0) THEN
                     tc   (jo,io,1) = tc   (jo,io,1) + gfcc_tc(j,i)*harea(j,i)
                     !树高加权
                     !粗网格城市树高=粗网格城市树高+500m城市格点植被覆盖度*城市格点树高*城市格点面积
                     !加权系数:城市格点植被覆盖度*城市面积
                     IF (gedi_th(j,i) > 0) THEN
                        htop(jo,io,1) = htop(jo,io,1) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)
                        wgt_top(jo,io,1) = wgt_top(jo,io,1) + gfcc_tc(j,i)*harea(j,i)
                     ENDIF
                  ENDIF
               ENDIF

               IF (urden(jr,ir) == 2 ) THEN
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
               ENDIF

               IF (urden(jr,ir) == 3 ) THEN
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
               ENDIF
            ENDIF

            IF (urden(jr,ir)<=0 .and. modur(j,i)>0) THEN
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
            ENDIF
         ENDDO
      ENDDO

      print *, "*** processing CLM 1km urban data ***"
      print *, "    assigning building properites for model grid"
      DO i = 1, 600
         DO j = 1, 600
            ! calculate io, jo
            io = NINT((urlats(i)+udelta/2._r8+ 90.)/y_delta+0.5)
            jo = NINT((urlonw(j)+udelta/2._r8+180.)/x_delta+0.5)

            ! IF (urden(j,i) > 0) THEN
            !    ur_land(jo,io) = ur_land(jo,io) + uarea(j,i)
            ! END IF

            ! 城市建筑属性聚合
            ! 加权：
            ! 粗网格城市属性=粗网格城市属性+细网格城市属性*细网格面积
            ! 加权系数：细网格面积 
            IF (urden(j,i) == 1) THEN
               ur_dc    (jo,io,1) = ur_dc(jo,io,1) + uarea(j,i)
               uxid               = urrgid(j,i)
               hwr_can  (jo,io,1) = hwr_can  (jo,io,1) + hwrcan  (1,uxid)*uarea(j,i)
               wt_rf    (jo,io,1) = wt_rf    (jo,io,1) + wtrf    (1,uxid)*uarea(j,i)
               wt_rd    (jo,io,1) = wt_rd    (jo,io,1) + wtrd    (1,uxid)*uarea(j,i)
               em_rf    (jo,io,1) = em_rf    (jo,io,1) + emrf    (1,uxid)*uarea(j,i)
               em_wl    (jo,io,1) = em_wl    (jo,io,1) + emwl    (1,uxid)*uarea(j,i)
               em_imrd  (jo,io,1) = em_imrd  (jo,io,1) + emimrd  (1,uxid)*uarea(j,i)
               em_perd  (jo,io,1) = em_perd  (jo,io,1) + emperd  (1,uxid)*uarea(j,i)
               !ulev_imrd(jo,io,1) = ulevimrd(1,uxid)
               th_rf    (jo,io,1) = th_rf    (jo,io,1) + thrf    (1,uxid)*uarea(j,i)
               th_wl    (jo,io,1) = th_wl    (jo,io,1) + thwl    (1,uxid)*uarea(j,i)
               tb_min   (jo,io,1) = tb_min   (jo,io,1) + tbmin   (1,uxid)*uarea(j,i)
               tb_max   (jo,io,1) = tb_max   (jo,io,1) + tbmax   (1,uxid)*uarea(j,i)
               ht_rf    (jo,io,1) = ht_rf    (jo,io,1) + htrf    (1,uxid)*uarea(j,i)
               !w_hc     (jo,io,1) = whc     (1,uxid)

               alb_rf  (jo,io,1,:,:) = alb_rf  (jo,io,1,:,:) + albrf  (1,uxid,:,:)*uarea(j,i)
               alb_wl  (jo,io,1,:,:) = alb_wl  (jo,io,1,:,:) + albwl  (1,uxid,:,:)*uarea(j,i)
               alb_imrd(jo,io,1,:,:) = alb_imrd(jo,io,1,:,:) + albimrd(1,uxid,:,:)*uarea(j,i)
               alb_perd(jo,io,1,:,:) = alb_perd(jo,io,1,:,:) + albperd(1,uxid,:,:)*uarea(j,i)

               tk_rf  (jo,io,1,:) = tk_rf  (jo,io,1,:) + tkrf  (1,uxid,:)*uarea(j,i)
               tk_wl  (jo,io,1,:) = tk_wl  (jo,io,1,:) + tkwl  (1,uxid,:)*uarea(j,i)
               DO m = 1, 10
                  ! tkimrd与cvimrd有缺省值，计算需要跳过
                  IF (tkimrd(1,uxid,m) .ne. -999) THEN
                     tk_imrd(jo,io,1,m) = tk_imrd(jo,io,1,m) + tkimrd(1,uxid,m)*uarea(j,i)
                  ENDIF
                  IF (cvimrd(1,uxid,m) .ne. -999.) THEN
                     cv_imrd(jo,io,1,m) = cv_imrd(jo,io,1,m) + cvimrd(1,uxid,m)*uarea(j,i)
                  ENDIF
               ENDDO
               cv_rf  (jo,io,1,:) = cv_rf  (jo,io,1,:) + cvrf  (1,uxid,:)*uarea(j,i)
               cv_wl  (jo,io,1,:) = cv_wl  (jo,io,1,:) + cvwl  (1,uxid,:)*uarea(j,i)
            ENDIF

            IF (urden(j,i) == 2) THEN
               ur_dc(jo,io,2) = ur_dc(jo,io,2) + uarea(j,i)
               uxid               = urrgid(j,i)
               hwr_can  (jo,io,2) = hwr_can  (jo,io,2) + hwrcan  (2,uxid)*uarea(j,i)
               wt_rf    (jo,io,2) = wt_rf    (jo,io,2) + wtrf    (2,uxid)*uarea(j,i)
               wt_rd    (jo,io,2) = wt_rd    (jo,io,2) + wtrd    (2,uxid)*uarea(j,i)
               em_rf    (jo,io,2) = em_rf    (jo,io,2) + emrf    (2,uxid)*uarea(j,i)
               em_wl    (jo,io,2) = em_wl    (jo,io,2) + emwl    (2,uxid)*uarea(j,i)
               em_imrd  (jo,io,2) = em_imrd  (jo,io,2) + emimrd  (2,uxid)*uarea(j,i)
               em_perd  (jo,io,2) = em_perd  (jo,io,2) + emperd  (2,uxid)*uarea(j,i)
               !ulev_imrd(jo,io,1) = ulevimrd(1,uxid)
               th_rf    (jo,io,2) = th_rf    (jo,io,2) + thrf    (2,uxid)*uarea(j,i)
               th_wl    (jo,io,2) = th_wl    (jo,io,2) + thwl    (2,uxid)*uarea(j,i)
               tb_min   (jo,io,2) = tb_min   (jo,io,2) + tbmin   (2,uxid)*uarea(j,i)
               tb_max   (jo,io,2) = tb_max   (jo,io,2) + tbmax   (2,uxid)*uarea(j,i)
               ht_rf    (jo,io,2) = ht_rf    (jo,io,2) + htrf    (2,uxid)*uarea(j,i)
               !w_hc     (jo,io,1) = whc     (1,uxid)

               alb_rf  (jo,io,2,:,:) = alb_rf  (jo,io,2,:,:) + albrf  (2,uxid,:,:)*uarea(j,i)
               alb_wl  (jo,io,2,:,:) = alb_wl  (jo,io,2,:,:) + albwl  (2,uxid,:,:)*uarea(j,i)
               alb_imrd(jo,io,2,:,:) = alb_imrd(jo,io,2,:,:) + albimrd(2,uxid,:,:)*uarea(j,i)
               alb_perd(jo,io,2,:,:) = alb_perd(jo,io,2,:,:) + albperd(2,uxid,:,:)*uarea(j,i)

               tk_rf  (jo,io,2,:) = tk_rf  (jo,io,2,:) + tkrf  (2,uxid,:)*uarea(j,i)
               tk_wl  (jo,io,2,:) = tk_wl  (jo,io,2,:) + tkwl  (2,uxid,:)*uarea(j,i)
               DO m = 1, 10
                  IF (tkimrd(2,uxid,m) .ne. -999) THEN
                     tk_imrd(jo,io,2,m) = tk_imrd(jo,io,2,m) + tkimrd(2,uxid,m)*uarea(j,i)
                  ENDIF
                  IF (cvimrd(2,uxid,m) .ne. -999.) THEN
                     cv_imrd(jo,io,2,m) = cv_imrd(jo,io,2,m) + cvimrd(2,uxid,m)*uarea(j,i)
                  ENDIF
               ENDDO
               cv_rf  (jo,io,2,:) = cv_rf  (jo,io,2,:) + cvrf  (2,uxid,:)*uarea(j,i)
               cv_wl  (jo,io,2,:) = cv_wl  (jo,io,2,:) + cvwl  (2,uxid,:)*uarea(j,i)
            ENDIF

            IF (urden(j,i) == 3) THEN
               ur_dc(jo,io,3) = ur_dc(jo,io,3) + uarea(j,i)
               uxid               = urrgid(j,i)
               hwr_can  (jo,io,3) = hwr_can  (jo,io,3) + hwrcan  (3,uxid)*uarea(j,i)
               wt_rf    (jo,io,3) = wt_rf    (jo,io,3) + wtrf    (3,uxid)*uarea(j,i)
               wt_rd    (jo,io,3) = wt_rd    (jo,io,3) + wtrd    (3,uxid)*uarea(j,i)
               em_rf    (jo,io,3) = em_rf    (jo,io,3) + emrf    (3,uxid)*uarea(j,i)
               em_wl    (jo,io,3) = em_wl    (jo,io,3) + emwl    (3,uxid)*uarea(j,i)
               em_imrd  (jo,io,3) = em_imrd  (jo,io,3) + emimrd  (3,uxid)*uarea(j,i)
               em_perd  (jo,io,3) = em_perd  (jo,io,3) + emperd  (3,uxid)*uarea(j,i)
               !ulev_imrd(jo,io,1) = ulevimrd(1,uxid)
               th_rf    (jo,io,3) = th_rf    (jo,io,3) + thrf    (3,uxid)*uarea(j,i)
               th_wl    (jo,io,3) = th_wl    (jo,io,3) + thwl    (3,uxid)*uarea(j,i)
               tb_min   (jo,io,3) = tb_min   (jo,io,3) + tbmin   (3,uxid)*uarea(j,i)
               tb_max   (jo,io,3) = tb_max   (jo,io,3) + tbmax   (3,uxid)*uarea(j,i)
               ht_rf    (jo,io,3) = ht_rf    (jo,io,3) + htrf    (3,uxid)*uarea(j,i)
               !w_hc     (jo,io,1) = whc     (1,uxid)

               alb_rf  (jo,io,3,:,:) = alb_rf  (jo,io,3,:,:) + albrf  (3,uxid,:,:)*uarea(j,i)
               alb_wl  (jo,io,3,:,:) = alb_wl  (jo,io,3,:,:) + albwl  (3,uxid,:,:)*uarea(j,i)
               alb_imrd(jo,io,3,:,:) = alb_imrd(jo,io,3,:,:) + albimrd(3,uxid,:,:)*uarea(j,i)
               alb_perd(jo,io,3,:,:) = alb_perd(jo,io,3,:,:) + albperd(3,uxid,:,:)*uarea(j,i)

               tk_rf  (jo,io,3,:) = tk_rf  (jo,io,3,:) + tkrf  (3,uxid,:)*uarea(j,i)
               tk_wl  (jo,io,3,:) = tk_wl  (jo,io,3,:) + tkwl  (3,uxid,:)*uarea(j,i)
               DO m = 1, 10
                  IF (tkimrd(3,uxid,m) .ne. -999.) THEN
                     tk_imrd(jo,io,3,m) = tk_imrd(jo,io,3,m) + tkimrd(3,uxid,m)*uarea(j,i)
                  ENDIF
                  IF (cvimrd(3,uxid,m) .ne. -999.) THEN
                     cv_imrd(jo,io,3,m) = cv_imrd(jo,io,3,m) + cvimrd(3,uxid,m)*uarea(j,i)
                  ENDIF
               ENDDO
               cv_rf  (jo,io,3,:) = cv_rf  (jo,io,3,:) + cvrf  (3,uxid,:)*uarea(j,i)
               cv_wl  (jo,io,3,:) = cv_wl  (jo,io,3,:) + cvwl  (3,uxid,:)*uarea(j,i)
            ENDIF
         ENDDO
      ENDDO

      ! MODIS城市补全，归类为MD类型城市
      DO i=1,1200
         DO j=1,1200
            ir = CEILING(i*1._r8/2)
            jr = CEILING(j*1._r8/2)
            io = NINT((hlat(i)+ 90.)/y_delta+0.5)
            jo = NINT((hlon(j)+180.)/x_delta+0.5)

            ! IF (io>360 .or. jo>720)THEN
            !    print*, i,j
            !    print*, io,jo
            !    print*, hlonw(j)
            ! ENDIF

            IF (urden(jr,ir)<=0 .and. modur(j,i)>0) THEN
               ur_dc(jo,io,3) = ur_dc(jo,io,3) + harea(j,i)*modur(j,i)/100
               uxid           = urrgid(jr,ir)
               IF (uxid == 0) THEN
                  uxid = lurrgid(jo,io)
                  IF (uxid == 0) THEN
                     print*, io, jo, modur(j,i)  
                     print*,i,j     
                     IF (io == 81 .and. jo == 498) THEN
                        uxid = 30
                     ENDIF
                     IF (io == 176 .and. jo == 198) THEN
                        uxid = 31
                     ENDIF
                     IF (io == 187 .and. jo == 601) THEN
                        uxid = 27
                     ENDIF
                     IF (io == 206 .and. jo == 220) THEN
                        uxid = 31
                     ENDIF
                     IF (io == 206 .and. jo == 221) THEN
                        uxid = 31
                     ENDIF
                     IF (io == 274 .and. jo == 248) THEN
                        uxid = 6
                     ENDIF
                  ENDIF
               ENDIF
               hwr_can  (jo,io,3) = hwr_can  (jo,io,3) + hwrcan  (3,uxid)*harea(j,i)*modur(j,i)/100
               wt_rf    (jo,io,3) = wt_rf    (jo,io,3) + wtrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               wt_rd    (jo,io,3) = wt_rd    (jo,io,3) + wtrd    (3,uxid)*harea(j,i)*modur(j,i)/100
               em_rf    (jo,io,3) = em_rf    (jo,io,3) + emrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               em_wl    (jo,io,3) = em_wl    (jo,io,3) + emwl    (3,uxid)*harea(j,i)*modur(j,i)/100
               em_imrd  (jo,io,3) = em_imrd  (jo,io,3) + emimrd  (3,uxid)*harea(j,i)*modur(j,i)/100
               em_perd  (jo,io,3) = em_perd  (jo,io,3) + emperd  (3,uxid)*harea(j,i)*modur(j,i)/100
               !ulev_imrd(jo,io,1) = ulevimrd(1,uxid)
               th_rf    (jo,io,3) = th_rf    (jo,io,3) + thrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               th_wl    (jo,io,3) = th_wl    (jo,io,3) + thwl    (3,uxid)*harea(j,i)*modur(j,i)/100
               tb_min   (jo,io,3) = tb_min   (jo,io,3) + tbmin   (3,uxid)*harea(j,i)*modur(j,i)/100
               tb_max   (jo,io,3) = tb_max   (jo,io,3) + tbmax   (3,uxid)*harea(j,i)*modur(j,i)/100
               ht_rf    (jo,io,3) = ht_rf    (jo,io,3) + htrf    (3,uxid)*harea(j,i)*modur(j,i)/100
               !w_hc     (jo,io,1) = whc     (1,uxid)

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
         ENDDO
      ENDDO
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
         ENDDO

         IF (avg(j,i) > 0) THEN
            hgt(j,i) = hgt(j,i) / avg(j,i)
         ENDIF 

         sumur = ur_dc(j,i,1) + ur_dc(j,i,2) + ur_dc(j,i,3)
         IF (sumur > 0.) THEN
            ! calculate 3 type urban cover, sum(pct_ur(j,i,:))=100
            pct_ur(j,i,1) = ur_dc(j,i,1) / sumur * 100.
            pct_ur(j,i,2) = ur_dc(j,i,2) / sumur * 100.
            pct_ur(j,i,3) = ur_dc(j,i,3) / sumur * 100.
         ENDIF

         ! 检查3类城市是否超过100%
         IF (sum(pct_ur(j,i,1:3)) > 1e-6 .and. abs(sum(pct_ur(j,i,1:3))-100) > 1e-3) THEN
            print *, 'urban_pct > 100'
            print *, pct_ur(j,i,1:3)
         ENDIF

      ENDDO
   ENDDO
   
   
   print *, "********************************"
   
   DO i = 1, nyo 
      DO j = 1, nxo
         ! calculate urban lai
         !-----------------------------------------------------
         !             ___npft
         ! urban_lai = \   pct_pft*pft_lai*(htop_gedi/htop_pft)
         !             /__ 1
         !-----------------------------------------------------      
         ! 只需要考虑树(1-8类PFTs)
         DO k =1, 3
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               htop_ur(j,i,k) = hgt(j,i)
            ENDIF

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

         DO m = 2, 9
            IF (htop_pft(j,i,m) > 0) THEN
               fac(1) = htop_ur(j,i,1) / htop_pft(j,i,m)!htop_ur(j,i,1)!htop_pft(j,i,m)
               fac(2) = htop_ur(j,i,2) / htop_pft(j,i,m)!htop_ur(j,i,1)!htop_pft(j,i,m)
               fac(3) = htop_ur(j,i,3) / htop_pft(j,i,m)!htop_ur(j,i,1)!htop_pft(j,i,m)

               wgt = pct_pft(j,i,m)/sum(pct_pft(j,i,2:9))
               
               DO k = 1, 12       
                  ur_lai (j,i,1,k) = ur_lai(j,i,1,k)+wgt*lai_pft(j,i,m,k)*min(1.5,fac(1))!/100!*min(fac,1.)/100.
                  ur_lai (j,i,2,k) = ur_lai(j,i,2,k)+wgt*lai_pft(j,i,m,k)*min(1.5,fac(2))!/100
                  ur_lai (j,i,3,k) = ur_lai(j,i,3,k)+wgt*lai_pft(j,i,m,k)*min(1.5,fac(3))!/100
                  ur_sai (j,i,1,k) = ur_sai(j,i,1,k)+wgt*sai_pft(j,i,m,k)*min(1.5,fac(1))!/100!*min(fac,1.)/100.
                  ur_sai (j,i,2,k) = ur_sai(j,i,2,k)+wgt*sai_pft(j,i,m,k)*min(1.5,fac(2))!/100
                  ur_sai (j,i,3,k) = ur_sai(j,i,3,k)+wgt*sai_pft(j,i,m,k)*min(1.5,fac(3))!/100 
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   ! data flip lat(0)=-90 -> lat(0)=90
   DO i = 1, nyo
      DO j = 1, nxo
         k = nyo - i + 1
         
         flatso (i)     = latso(k)
         fpct_ur(j,i,:) = pct_ur(j,k,:)
         fpct_tc(j,i,:) = pct_tc(j,k,:)
         fpct_urwt(j,i,:)   = pct_urwt(j,k,:)
         fhtop_ur (j,i,:)   = htop_ur(j,k,:)
         fur_lai  (j,i,:,:) = ur_lai(j,k,:,:)
         fur_sai  (j,i,:,:) = ur_sai(j,k,:,:)
         fhwr_can (j,i,:) = hwr_can(j,k,:)
         fwt_rf   (j,i,:) = wt_rf(j,k,:)
         fwt_rd   (j,i,:) = wt_rd(j,k,:)
         fem_rf   (j,i,:) = em_rf(j,k,:)
         fem_wl   (j,i,:) = em_wl(j,k,:)
         fem_imrd (j,i,:) = em_imrd(j,k,:)
         fem_perd (j,i,:) = em_perd(j,k,:)
         fht_rf   (j,i,:) = ht_rf(j,k,:)
         fth_rf   (j,i,:) = th_rf(j,k,:)
         fth_wl   (j,i,:) = th_wl(j,k,:)
         ftb_min  (j,i,:) = tb_min(j,k,:)
         ftb_max  (j,i,:) = tb_max(j,k,:)
         ftk_rf   (j,i,:,:) = tk_rf(j,k,:,:) 
         ftk_wl   (j,i,:,:) = tk_wl(j,k,:,:)
         ftk_imrd (j,i,:,:) = tk_imrd(j,k,:,:)
         fcv_rf   (j,i,:,:) = cv_rf(j,k,:,:)
         fcv_wl   (j,i,:,:) = cv_wl(j,k,:,:)
         fcv_imrd (j,i,:,:) = cv_imrd(j,k,:,:)
         falb_rf  (j,i,:,:,:) = alb_rf(j,k,:,:,:)
         falb_wl  (j,i,:,:,:) = alb_wl(j,k,:,:,:)
         falb_imrd(j,i,:,:,:) = alb_imrd(j,k,:,:,:)
         falb_perd(j,i,:,:,:) = alb_perd(j,k,:,:,:)
      ENDDO
   ENDDO

   print *, "*****************************************"

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

   CALL check( nf90_create("scolm_urban_data_modis.nc", NF90_NETCDF4, ncid) )

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
   !CALL check( nf90_def_var(ncid, "REGION_ID"      , NF90_INT  , XY2D, ur_rgvid   ) )
   CALL check( nf90_def_var(ncid, "CANYON_HWR"     , NF90_FLOAT, XY3D, hwr_canvid ) )
   CALL check( nf90_def_var(ncid, "WTLUNIT_ROOF"   , NF90_FLOAT, XY3D, wt_rfvid   ) )
   CALL check( nf90_def_var(ncid, "WTROAD_PERV"    , NF90_FLOAT, XY3D, wt_rdvid   ) )
   CALL check( nf90_def_var(ncid, "EM_ROOF"        , NF90_FLOAT, XY3D, em_rfvid   ) )
   CALL check( nf90_def_var(ncid, "EM_WALL"        , NF90_FLOAT, XY3D, em_wlvid   ) )
   CALL check( nf90_def_var(ncid, "EM_IMPROAD"     , NF90_FLOAT, XY3D, em_imrdvid ) )
   CALL check( nf90_def_var(ncid, "EM_PERROAD"     , NF90_FLOAT, XY3D, em_perdvid ) )
   CALL check( nf90_def_var(ncid, "HT_ROOF"        , NF90_FLOAT, XY3D, ht_rfvid   ) )
   !CALL check( nf90_def_var(ncid, "WIND_HGT_CANYON", NF90_FLOAT, XY3D, whcvid     ) )
   !CALL check( nf90_def_var(ncid, "NLEV_IMPROAD"   , NF90_FLOAT, XY3D, ulev_vid   ) )
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

   CALL check( nf90_put_att(ncid, pct_urvid , "units"     , "%"                                           ) )
   CALL check( nf90_put_att(ncid, pct_urvid , "long_name" , "Percentage of each urban type (density)"     ) )
   CALL check( nf90_put_att(ncid, pct_urvid , "_FillValue", fv                   ) )

   !CALL check( nf90_put_att(ncid, ur_rgvid  , "units"     , "unitless"           ) )
   !CALL check( nf90_put_att(ncid, ur_rgvid  , "long_name" , "Region ID"          ) )

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
 
   !CALL check( nf90_put_att(ncid, whcvid    , "units"     , "meters"                                     ) )
   !CALL check( nf90_put_att(ncid, whcvid    , "long_name" , "height of wind in canyon"                   ) )
   !CALL check( nf90_put_att(ncid, whcvid    , "_FillValue", fv) )

   !CALL check( nf90_put_att(ncid, ulev_vid  , "units"     , "unitless"                                   ) )
   !CALL check( nf90_put_att(ncid, ulev_vid  , "long_name" , "number of impervious road layers"           ) )
   !CALL check( nf90_put_att(ncid, ulev_vid  , "_FillValue", fv) )

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
   
   CALL check( nf90_inq_varid(ncid, "lat"     , urlat_vid       ) )
   CALL check( nf90_put_var  (ncid, urlat_vid , flatso          ) )

   CALL check( nf90_inq_varid(ncid, "lon"     , urlon_vid       ) )
   CALL check( nf90_put_var  (ncid, urlon_vid , lonso           ) )

   call check( nf90_inq_varid(ncid, "numsolar", ns_vid         ) )
   CALL check( nf90_put_var  (ncid, ns_vid    , n_ns           ) )

   CALL check( nf90_inq_varid(ncid, "numrad"  , nr_vid           ) )
   CALL check( nf90_put_var  (ncid, nr_vid    , n_nr             ) )

   CALL check( nf90_inq_varid(ncid, "month"   , mon_vid         ) )
   CALL check( nf90_put_var  (ncid, mon_vid   , n_mon           ) )

   CALL check( nf90_inq_varid(ncid, "density" , den_vid         ) )
   CALL check( nf90_put_var  (ncid, den_vid   , n_den           ) )

   CALL check( nf90_inq_varid(ncid, "ulev"    , lev_vid) )
   CALL check( nf90_put_var  (ncid, lev_vid   , n_ulev) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_PCT", pct_tcvid      ) )
   CALL check( nf90_put_var  (ncid, pct_tcvid , fpct_tc         ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_WATER_PCT" , pct_urwtvid   ) )
   CALL check( nf90_put_var  (ncid, pct_urwtvid, fpct_urwt      ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_TOP" , htop_urvid     ) )
   CALL check( nf90_put_var  (ncid, htop_urvid, fhtop_ur        ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_LAI", ur_laivid     ) )
   CALL check( nf90_put_var  (ncid, ur_laivid  , fur_lai        ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_TREE_SAI", ur_saivid     ) )
   CALL check( nf90_put_var  (ncid, ur_saivid  , fur_sai        ) )

   !CALL check( nf90_inq_varid(ncid, "REGION_ID", ur_rgvid      ) )
   !CALL check( nf90_put_var  (ncid,ur_rgvid, ur_rgid           ) )

   CALL check( nf90_inq_varid(ncid, "CANYON_HWR", hwr_canvid   ) )
   CALL check( nf90_put_var  (ncid,hwr_canvid,fhwr_can          ) )

   CALL check( nf90_inq_varid(ncid, "WTLUNIT_ROOF",wt_rfvid    ) )
   CALL check( nf90_put_var(ncid,wt_rfvid,fwt_rf                ) )

   CALL check( nf90_inq_varid(ncid, "WTROAD_PERV",wt_rdvid     ) )
   CALL check( nf90_put_var(ncid,wt_rdvid, fwt_rd               ) )

   CALL check( nf90_inq_varid(ncid, "EM_ROOF",em_rfvid         ) )
   CALL check( nf90_put_var(ncid,em_rfvid,fem_rf                ) )

   CALL check( nf90_inq_varid(ncid, "EM_WALL",em_wlvid         ) )
   CALL check( nf90_put_var(ncid,em_wlvid,fem_wl                ) )

   CALL check( nf90_inq_varid(ncid, "EM_IMPROAD",em_imrdvid    ) )
   CALL check( nf90_put_var(ncid,em_imrdvid,fem_imrd            ) )

   CALL check( nf90_inq_varid(ncid, "EM_PERROAD",em_perdvid    ) )
   CALL check( nf90_put_var(ncid,em_perdvid,fem_perd            ) )

   CALL check( nf90_inq_varid(ncid, "HT_ROOF",ht_rfvid         ) )
   CALL check( nf90_put_var(ncid,ht_rfvid,fht_rf                ) )

   !CALL check( nf90_inq_varid(ncid, "WIND_HGT_CANYON",whcvid   ) )
   !CALL check( nf90_put_var(ncid,whcvid,w_hc                   ) )

   CALL check( nf90_inq_varid(ncid, "THICK_ROOF",th_rfvid      ) )
   CALL check( nf90_put_var(ncid,th_rfvid,fth_rf                ) )

   CALL check( nf90_inq_varid(ncid, "THICK_WALL",th_wlvid      ) )
   CALL check( nf90_put_var(ncid,th_wlvid,fth_wl                ) )

   CALL check( nf90_inq_varid(ncid, "T_BUILDING_MIN",tbminvid  ) )
   CALL check( nf90_put_var(ncid,tbminvid,ftb_min               ) )

   CALL check( nf90_inq_varid(ncid, "T_BUILDING_MAX",tbmaxvid  ) )
   CALL check( nf90_put_var(ncid,tbmaxvid,ftb_max               ) )

   !CALL check( nf90_inq_varid(ncid, "NLEV_IMPROAD",ulev_imrdvid) )
   !CALL check( nf90_put_var(ncid,ulev_imrdvid,ulev_imrd        ) )

   CALL check( nf90_inq_varid(ncid, "URBAN_PCT", pct_urvid     ) )
   CALL check( nf90_put_var(ncid,pct_urvid,fpct_ur              ) )

   CALL check( nf90_inq_varid(ncid, "TK_ROOF",tk_rfvid         ) )
   CALL check( nf90_put_var(ncid,tk_rfvid,ftk_rf                ) )

   CALL check( nf90_inq_varid(ncid, "TK_WALL",tk_wlvid         ) )
   CALL check( nf90_put_var(ncid,tk_wlvid,ftk_wl                ) )

   CALL check( nf90_inq_varid(ncid, "TK_IMPROAD",tk_imrdvid    ) )
   CALL check( nf90_put_var(ncid,tk_imrdvid,ftk_imrd            ) )

   CALL check( nf90_inq_varid(ncid, "CV_ROOF",cv_rfvid         ) )
   CALL check( nf90_put_var(ncid,cv_rfvid,fcv_rf                ) )

   CALL check( nf90_inq_varid(ncid, "CV_WALL",cv_wlvid         ) )
   CALL check( nf90_put_var(ncid,cv_wlvid,fcv_wl                ) )

   CALL check( nf90_inq_varid(ncid, "CV_IMPROAD",cv_imrdvid    ) )
   CALL check( nf90_put_var(ncid,cv_imrdvid,fcv_imrd            ) )

   CALL check( nf90_inq_varid(ncid, "ALB_ROOF",alb_rfvid       ) )
   CALL check( nf90_put_var(ncid,alb_rfvid,falb_rf              ) )

   CALL check( nf90_inq_varid(ncid, "ALB_WALL",alb_wlvid       ) )
   CALL check( nf90_put_var(ncid,alb_wlvid,falb_wl              ) )

   CALL check( nf90_inq_varid(ncid, "ALB_IMPROAD",alb_imrdvid  ) )
   CALL check( nf90_put_var(ncid,alb_imrdvid,falb_imrd          ) )

   CALL check( nf90_inq_varid(ncid, "ALB_PERROAD",alb_perdvid  ) )
   CALL check( nf90_put_var(ncid,alb_perdvid,falb_perd          ) )

   CALL check( nf90_close(ncid) )

   !CALL check( nf90_create("hgt.nc", NF90_NETCDF4, ncid) )

   !CALL check( nf90_def_dim(ncid, "lat"     , nyo     , lat_dimid ) )
   !CALL check( nf90_def_dim(ncid, "lon"     , nxo     , lon_dimid ) )

   !CALL check( nf90_def_var(ncid, "lat", nf90_double, (/lat_dimid/), lat_vid) )
   !CALL check( nf90_def_var(ncid, "lon", nf90_double, (/lon_dimid/), lon_vid) )
   !CALL check( nf90_def_var(ncid, "hgt", nf90_double, (/lon_dimid,lat_dimid/), p) )
   !CALL check( nf90_inq_varid(ncid, "lat", lat_vid) )
   !CALL check( nf90_put_var(ncid, lat_vid, latso) )
   !CALL check( nf90_inq_varid(ncid, "lon", lon_vid) )
   !CALL check( nf90_put_var(ncid, lon_vid, lonso) )
   !CALL check( nf90_inq_varid(ncid, "hgt", p) )
   !CALL checK( nf90_put_var(ncid, p, hgt) )
   !CALL check( nf90_close(ncid) )
   print*, "*** SUCCESS write surface file ***"
   
   deallocate( tlat   )
   deallocate( tlon   )
   deallocate( hlat   )
   deallocate( hlats  )
   deallocate( hlatn  )
   deallocate( hlon   )
   deallocate( hlonw  )
   deallocate( hlone  )
   deallocate( urlat  )
   deallocate( nrlat  )
   deallocate( urlats )
   deallocate( urlatn )
   deallocate( urlon  )
   deallocate( urlonw )
   deallocate( urlone )
   deallocate( harea  )
   deallocate( mask   )
   deallocate( gfcc_tc)
   deallocate( gedi_th)
   deallocate( gl30_wt)

   deallocate( uarea  )
   deallocate( umask  )
   deallocate( umd    )
   deallocate( nmd    )
   deallocate( uhd    )
   deallocate( lurrgid    )
   deallocate( utbd   )
   deallocate( ntbd   )
   deallocate( urclss )
   deallocate( urrgid )

   deallocate( latso     )
   deallocate( lonso     )
   deallocate( area      )
   deallocate( pct_urban )
   deallocate( npct_urban)
   deallocate( htop_pft  )
   deallocate( nhtop_pft )
   deallocate( pct_pft   )
   deallocate( npct_pft  )
   deallocate( lai_pft   )
   deallocate( nlai_pft  )
   deallocate( ur_dc     )
   deallocate( pct_ur    )

   deallocate( cnt       )
   deallocate( pct_tc    )
   deallocate( pct_tland )
   deallocate( pct_wland )
   deallocate( pct_hland )
   deallocate( pct_urwt  )
   deallocate( htop_ur   )
   deallocate( ur_land   )
   deallocate( ur_clss   )
   deallocate( ur_rgid   )
   deallocate( ur_lai    )
   deallocate( hwr_can   )
   deallocate( wt_rf     )
   deallocate( wt_rd     )
   deallocate( em_rf     )
   deallocate( em_wl     )
   deallocate( em_imrd   )
   deallocate( em_perd   )
   deallocate( ht_rf     )
   deallocate( w_hc      )
   deallocate( ulev_imrd )
   deallocate( th_rf     )
   deallocate( th_wl     )
   deallocate( tb_min    )
   deallocate( tb_max    )
   deallocate( cv_rf     )
   deallocate( cv_wl     )
   deallocate( cv_imrd   )
   deallocate( tk_rf     )
   deallocate( tk_wl     )
   deallocate( tk_imrd   )
   deallocate( alb_rf    )
   deallocate( alb_wl    )
   deallocate( alb_imrd  )
   deallocate( alb_perd  )

CONTAINS 
   
   SUBROUTINE check(status)
      INTEGER, intent(in) :: status

      IF (status /= nf90_noerr) THEN
         print *, trim( nf90_strerror(status))
         stop 2
      ENDIF 
   END SUBROUTINE check

END PROGRAM clmu2grid
