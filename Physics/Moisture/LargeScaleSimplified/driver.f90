!TOP LEVEL DRIVER, LOADS IN VARIABLES AND SENDS TO SUBROUTINE
!AS REQUIRED BY TAF, ALL CALCULATIONS ARE DONE IN A SUBROUTINE

program driver

!Modules to Include
use netcdf

IMPLICIT NONE

LOGICAL, PARAMETER :: write_data = .true.

! THE PARAMETERS
INTEGER, PARAMETER :: kmax=72

INTEGER, PARAMETER :: LON_RES = 144
INTEGER, PARAMETER :: imax=144, jmax=91
REAl, PARAMETER :: DT = 1800.0

!ARES TO STUDY
INTEGER, PARAMETER :: i_colmin = 1, j_colmin = 1
INTEGER, PARAMETER :: i_colmax = imax, j_colmax = jmax

!INTEGER, PARAMETER :: LON_RES = 576
!INTEGER, PARAMETER :: imax=576, jmax=361
!REAl, PARAMETER :: DT = 1200.0

INTEGER :: i, j, k, kx, l

!MODEL VARIABLES FROM DATA SET
REAL :: ak(kmax+1), bk(kmax+1), PREF(kmax+1)

!MODEL VARIABLES PROVIDED FROM RUN
integer :: rcode, ncID, dimID, varID
integer :: nlat, nlon, nlev!, ntim
integer :: x_dimid, y_dimid, z_dimid
integer :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8

character(len=64) :: input_file, output_file
character(len=64) :: charbuf

integer, allocatable, dimension(:) :: lev
real, allocatable, dimension(:) :: lat, lon

real, allocatable, dimension(:,:,:) :: Q, TH, U, V
real, allocatable, dimension(:,:) :: Ps, Ts, FRLAND
integer, allocatable, dimension(:,:) :: KCBL

REAL, DIMENSION(imax,jmax,kmax) :: Q_OUT, TH_OUT, U_OUT, V_OUT, P_OUT

REAL, DIMENSION(imax,jmax) :: ls_precip

!==============================================================================
!              BEGIN READ IN MODEL VARIABLES AND CONSTANTS       
!==============================================================================

! READ IN REFERENCE PRESSURE AND COMP SIGE AND ICMIN
open (11,file='pref.txt')
do k=1,73
   read (11,'(f9.4)') PREF(k)
enddo
 close (11)
PREF = PREF/0.01 !Convert from hPa to Pa

! READ IN COORDINATE DETAILS
open (12,file='akbk72.txt')
do k=1,73
   read (12,'(i3,2f13.6)') kx, ak(k), bk(k)
enddo
 close (12)

! READ IN MODEL VARIABLES FROM NETCDFs
CALL ChDir("/discover/nobackup/drholdaw/inputs_forRAS/high_vs_low")

! Read a nc4 file
if (LON_RES == 144) THEN
   input_file = 'high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4'
elseif (LON_RES == 576) THEN
   input_file = 'high_vs_low_576361.geosgcm_rasinputs3d.20100702_1200z.nc4'
endif
rcode = nf90_open(trim(adjustl(input_file)), nf90_nowrite, ncID)

rcode = nf90_inq_dimid(ncID, "lev", dimID)
rcode = nf90_inquire_dimension(ncID, dimID, charbuf, nlev)
rcode = nf90_inq_dimid(ncID, "lon", dimID)
rcode = nf90_inquire_dimension(ncID, dimID, charbuf, nlon)
rcode = nf90_inq_dimid(ncID, "lat", dimID)
rcode = nf90_inquire_dimension(ncID, dimID, charbuf, nlat)

allocate(lat(nlat))
rcode = nf90_inq_varid(ncID, 'lat', varID)
rcode = nf90_get_var(ncID, varID, lat) 

allocate(lon(nlon))
rcode = nf90_inq_varid(ncID, 'lon', varID)
rcode = nf90_get_var(ncID, varID, lon) 

allocate(lev(nlev))
rcode = nf90_inq_varid(ncID, 'lev', varID)
rcode = nf90_get_var(ncID, varID, lev) 

allocate(Q(nlon,nlat,nlev))
rcode = nf90_inq_varid(ncID, 'Q_preras', varID)
rcode = nf90_get_var(ncID, varID, Q)

allocate(TH(nlon,nlat,nlev))
rcode = nf90_inq_varid(ncID, 'TH_preras', varID)
rcode = nf90_get_var(ncID, varID, TH)

allocate(U(nlon,nlat,nlev))
rcode = nf90_inq_varid(ncID, 'U_preras', varID)
rcode = nf90_get_var(ncID, varID, U)

allocate(V(nlon,nlat,nlev))
rcode = nf90_inq_varid(ncID, 'V_preras', varID)
rcode = nf90_get_var(ncID, varID, V)

rcode = nf90_close(ncid)

if (LON_RES == 144) THEN
   input_file = 'high_vs_low_14491.geosgcm_rasinputs2d.20100702_1200z.nc4'
elseif (LON_RES == 576) THEN
   input_file = 'high_vs_low_576361.geosgcm_rasinputs2d.20100702_1200z.nc4'
endif
rcode = nf90_open(trim(adjustl(input_file)), nf90_nowrite, ncID)

allocate(Ts(nlon,nlat))
rcode = nf90_inq_varid(ncID, 'TS_preras', varID)
rcode = nf90_get_var(ncID, varID, Ts)

allocate(KCBL(nlon,nlat))
rcode = nf90_inq_varid(ncID, 'KCBL', varID)
rcode = nf90_get_var(ncID, varID, KCBL)

rcode = nf90_close(ncid)

if (LON_RES == 144) THEN
   input_file = 'high_vs_low_14491.geosgcm_surf.20100702_1330z.nc4'
elseif (LON_RES == 576) THEN
   input_file = 'high_vs_low_576361.geosgcm_surf.20100702_1330z.nc4'
endif
rcode = nf90_open(trim(adjustl(input_file)), nf90_nowrite, ncID)

allocate(FRLAND(nlon,nlat))
rcode = nf90_inq_varid(ncID, 'FRLAND', varID)
rcode = nf90_get_var(ncID, varID, FRLAND)

rcode = nf90_close(ncid)
if (LON_RES == 144) THEN
   input_file = 'high_vs_low_14491.geosgcm_prog.20100702_1200z.nc4'
elseif (LON_RES == 576) THEN
   input_file = 'high_vs_low_576361.geosgcm_prog.20100702_1200z.nc4'
endif
rcode = nf90_open(trim(adjustl(input_file)), nf90_nowrite, ncID)

allocate(Ps(nlon,nlat))
rcode = nf90_inq_varid(ncID, 'PS', varID)
rcode = nf90_get_var(ncID, varID, Ps)
Ps = 0.01*Ps

rcode = nf90_close(ncid)

CALL ChDir('/home/drholdaw/moist_adjoint/largecloud')

!==============================================================================
!                  END READ IN MODEL VARIABLES AND CONSTANTS       
!==============================================================================

call global_driver(imax,jmax,kmax,i_colmin,j_colmin,i_colmax,j_colmax,DT,U,V,TH,Q,Ps,KCBL,Ts,FRLAND,ak,bk,PREF,&
                   U_OUT,V_OUT,TH_OUT,Q_OUT,ls_precip,P_OUT)

if (write_data == .true.) then

   CALL ChDir("/discover/nobackup/drholdaw/ExperimentData/moistoutput")

   output_file = '3DFields_largecloud.nc4'
   rcode = nf90_create(trim(adjustl(output_file)), NF90_CLOBBER, ncid)
   rcode = nf90_def_dim(ncid, "x", jmax, x_dimid)
   rcode = nf90_def_dim(ncid, "y", imax, y_dimid)
   rcode = nf90_def_dim(ncid, "z", kmax, z_dimid)

   rcode = nf90_def_var(ncid, "Q",  NF90_DOUBLE,(/ y_dimid, x_dimid , z_dimid /), varid1)
   rcode = nf90_def_var(ncid, "TH", NF90_DOUBLE,(/ y_dimid, x_dimid , z_dimid /), varid3)
   rcode = nf90_def_var(ncid, "U",  NF90_DOUBLE,(/ y_dimid, x_dimid , z_dimid /), varid4)
   rcode = nf90_def_var(ncid, "V",  NF90_DOUBLE,(/ y_dimid, x_dimid , z_dimid /), varid5)
   rcode = nf90_def_var(ncid, "P",  NF90_DOUBLE,(/ y_dimid, x_dimid , z_dimid /), varid6)
   rcode = nf90_def_var(ncid, "R",  NF90_DOUBLE,(/ y_dimid, x_dimid /), varid7)

   rcode = nf90_enddef(ncid)
              
   rcode = nf90_put_var(ncid, varid1, Q_out(:,:,:))
   rcode = nf90_put_var(ncid, varid3, TH_out(:,:,:))
   rcode = nf90_put_var(ncid, varid4, U_out(:,:,:))
   rcode = nf90_put_var(ncid, varid5, V_out(:,:,:))
   rcode = nf90_put_var(ncid, varid6, P_out(:,:,:))
   rcode = nf90_put_var(ncid, varid7, ls_precip(:,:))


   rcode = nf90_close(ncid)
   
   CALL ChDir("/home/moist_adjoint/largecloud")

endif



!DEALLOCATE TO FREE MEMORY
deallocate(lat)
deallocate(lon)
deallocate(lev)
deallocate(Q)
deallocate(TH)
deallocate(U)
deallocate(V)
deallocate(Ts)
deallocate(KCBL)
deallocate(FRLAND)
deallocate(Ps)


end program





