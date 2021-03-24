program Jacobian_con


!Modules to Include
use MAPL_ConstantsMod
use netcdf

IMPLICIT NONE

LOGICAL, PARAMETER :: write_data = .true., COMP_EVALS = .false., COMP_RMSE = .true., SAVEHQ = .true.

! THE PARAMETERS
INTEGER, PARAMETER :: kmax=72

INTEGER, PARAMETER :: LON_RES = 144
INTEGER, PARAMETER :: imax=144, jmax=91
REAl, PARAMETER :: DT = 1800.0

INTEGER :: i, j, k, kx, l, pert_point

INTEGER :: ii, jj, kk, ll

INTEGER :: icol, jcol

!POINT TO STUDY
INTEGER, PARAMETER :: global_points = imax*jmax
INTEGER :: num_active_points
INTEGER :: active_points(imax*jmax,6)

!MODEL VARIABLES FROM DATA SET
REAL :: ak(kmax+1), bk(kmax+1), PREF(kmax+1)

!MODEL VARIABLES PROVIDED FROM RUN
integer :: rcode, ncID, dimID, varID
integer :: nlat, nlon, nlev!, ntim
integer :: x_dimid, y_dimid, z_dimid
integer :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8
integer :: varid9, varid10, varid11, varid12, varid13, varid14, varid15, varid16

character(len=64) :: input_file, output_file
character(len=64) :: charbuf

integer, allocatable, dimension(:) :: lev
real, allocatable, dimension(:) :: lat, lon

real, allocatable, dimension(:,:,:) :: U, V, TH, Q
real,              dimension(imax,jmax,0:kmax) :: P
real, allocatable, dimension(:,:) :: Ps, Ts, FRLAND
integer, allocatable, dimension(:,:) :: KCBL

REAL, DIMENSION(imax,jmax,kmax) :: PRC3_RAS
REAL, DIMENSION(imax,jmax,kmax) :: PRC3_LSC

INTEGER, PARAMETER :: n_ras_params = 25
REAL :: RASPARAMS(n_ras_params)

!STANDARD DEVIATION DATA
REAL, DIMENSION(kmax) :: Ps_sd, Ph_lev, TH_sd, T_sd, Q_sd, U_sd, V_sd
REAL, DIMENSION(0:kmax) :: P_sd
INTEGER, PARAMETER :: pert_vec_size = 14
REAL :: pert_vec(pert_vec_size), pert

!SAVE INPUT VARIABLES
real, dimension(imax,jmax,kmax) :: U_in, V_in, TH_in, Q_in
real, dimension(imax,jmax,0:kmax) :: P_in
real, dimension(imax,jmax) :: Ps_in

!SAVE INPUT VARIABLES
real, dimension(kmax) :: U_pert_in, V_pert_in, TH_pert_in, Q_pert_in
real, dimension(0:kmax) :: P_pert_in

!TENDANCIES THROUGH FIRST CALL
real, dimension(imax,jmax,kmax) :: DUDT, DVDT, DTHDT, DQDT

!PERTURBED INPUTS
real, dimension(kmax) :: U_pert, V_pert, TH_pert, Q_pert
real :: pert_amount_var
real, dimension(0:kmax) :: P_pert

REAL, DIMENSION(kmax) :: PRC3_RAS_pert, PRC3_LSC_pert

!TENDANCIES THROUGH PERT CALL
real, dimension(kmax) :: DUDT_pert, DVDT_pert, DTHDT_pert, DQDT_pert

real, dimension(pert_vec_size,5*kmax,5*kmax) :: JACOBIAN

!EIGENVALUE COMPUTATION
real, dimension(4*kmax,4*kmax) :: EYE, A
REAL :: REVALS_MAX(global_points,pert_vec_size+1), IEVALS_MAX(global_points,pert_vec_size+1)

INTEGER, PARAMETER :: N = 4*kmax
INTEGER, PARAMETER :: LDA = N, LDVL = N, LDVR = N
INTEGER, PARAMETER :: LWMAX = 10000

INTEGER :: INFO, LWORK

REAL :: VL(LDVL,N), VR(LDVR,N), WR(N), WI(N), WORK(LWMAX)

EXTERNAL :: DGEEV

!RMSE COMPUTATION
REAL :: sum2, sumn, RMSE(global_points,pert_vec_size)


!==============================================================================
!              BEGIN READ IN MODEL VARIABLES AND CONSTANTS       
!==============================================================================

CALL Chdir("/home/drholdaw/Lin_Moist_Physics/Inputs")

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

!LOAD IN STANDARD DEVIATION DATA
open(13,file='errorvariance72')

DO i = 1,kmax
   read (13,'(f9.3,f8.3,f7.3,f12.8,2f8.2)') Ph_lev(i), Ps_sd(i), T_sd(i), Q_sd(i), U_sd(i), V_sd(i)
endDO
 close(13)
!Approximate P_sd through layer
P_sd(1:kmax) = (Ps_sd(kmax)/Ph_lev(kmax))*Ph_lev

Ps_sd = Ps_sd*0.01
P_sd = P_sd*0.01
Ph_lev = Ph_lev*0.01

!Convert to give TH s.d.
TH_sd = (1000.0/Ph_lev)**MAPL_KAPPA * T_sd

!Amounts to perturb by (pert_amoutn * sd)
pert_vec =  (/ 0.0001, -2.0, -1.0, -0.5, -0.1, -0.01, -0.001, -0.0001, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0 /)

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

CALL ChDir('/home/drholdaw/Lin_Moist_Physics/Jacobian_Pert')

!DERIVE P, WHICH IS WILL BE IN THE TRAJECTORY OF THE MODEL
DO k = 0,kmax
   P(:,:,k) = ak(k+1)*.01 + bk(k+1)*Ps(:,:)
endDO

!SAVE PRE CALL PRONOSTIC VARIABLES
U_in = U
V_in = V
TH_in = TH
Q_in = Q
P_in = P
Ps_in = Ps

 !RASPARAMS
 RASPARAMS( 1) = 1.000
 RASPARAMS( 2) = 0.05
 !if (FRLAND<0.1) then
    RASPARAMS( 3) = 2.5e-3   ! ocean value
 !else
 !   RASPARAMS( 3) = 2.5e-3  ! land value
 !endif
 RASPARAMS( 4) = 8.0e-4
 RASPARAMS( 5) = 1800.
 RASPARAMS( 6) = 43200.0
 RASPARAMS( 7) = -300.
 RASPARAMS( 8) = 4.0 
 RASPARAMS( 9) = 0.0 
 RASPARAMS(10) = 200.
 RASPARAMS(11) = 7.5e-4
 RASPARAMS(12) = 1.0 
 RASPARAMS(13) =-1.0 
 RASPARAMS(14) = 1.3 
 RASPARAMS(15) = 1.3 
 RASPARAMS(16) = 263.
 RASPARAMS(17) = 0.5 
 RASPARAMS(18) = 1.0 
 RASPARAMS(19) = 0.0 
 RASPARAMS(20) = 0.1 
 RASPARAMS(21) = 0.8 
 RASPARAMS(22) = 1.0
 !if (imax==144) then
    RASPARAMS(23) = 4000.0   !Low Res
 !elseif (imax==576) then
 !   RASPARAMS(23) = 1000.0   !High Res
 !endif
 RASPARAMS(24) = 0.5
 RASPARAMS(25) = 0.65

!CALL THE NONLINEAR SCHEME, COLUMN-BY-COLUMN, FOR THE WHOLE GLOBE
!DO icol = 1,imax
!   DO jcol = 1,jmax

DO icol = 116,116
   DO jcol = 61,61


      call convection_driver(kmax,DT,KCBL(icol,jcol),Ts(icol,jcol),FRLAND(icol,jcol),PREF,&
                             U(icol,jcol,:),V(icol,jcol,:),TH(icol,jcol,:),Q(icol,jcol,:),P(icol,jcol,:),&
                             PRC3_RAS(icol,jcol,:),PRC3_LSC(icol,jcol,:),RASPARAMS,TH_in(icol,jcol,:),Q_in(icol,jcol,:))

   endDO
endDO

!COMPUTE GLOBAL RATES OF CHANGE
DUDt  = (U - U_in)/DT
DVDt  = (V - V_in)/DT
DTHDt = (TH - TH_in)/DT
DQDt  = (Q - Q_in)/DT

print*, DUDt(116,61,:)

IF (SAVEHQ == .true.) THEN
   CALL ChDir('/home/drholdaw/Lin_Moist_Physics/RAS_Jacobian_Paper/mfiles')   
   output_file = 'HQ.nc4'
   rcode = nf90_create(trim(adjustl(output_file)), NF90_CLOBBER, ncid)
   rcode = nf90_def_dim(ncid, "x", jmax, x_dimid) 
   rcode = nf90_def_dim(ncid, "z", kmax, z_dimid) 
   rcode = nf90_def_dim(ncid, "y", imax, y_dimid)
   rcode = nf90_def_var(ncid, "DUDt", NF90_DOUBLE,  (/ y_dimid, x_dimid, z_dimid /), varid1)
   rcode = nf90_def_var(ncid, "DVDt", NF90_DOUBLE,  (/ y_dimid, x_dimid, z_dimid /), varid2)
   rcode = nf90_def_var(ncid, "DTHDt", NF90_DOUBLE, (/ y_dimid, x_dimid, z_dimid /), varid3)
   rcode = nf90_def_var(ncid, "DQDt" , NF90_DOUBLE, (/ y_dimid, x_dimid, z_dimid /), varid4)   
   rcode = nf90_enddef(ncid)
   rcode = nf90_put_var(ncid, varid1, DUDt)
   rcode = nf90_put_var(ncid, varid2, DVDt)
   rcode = nf90_put_var(ncid, varid3, DTHDt)
   rcode = nf90_put_var(ncid, varid4, DQDt)
   rcode = nf90_close(ncid)
endIF


!FIND POINTS WHERE SCHEME WAS ACTIVE
num_active_points = 0
DO i = 1,imax
   DO j = 1,jmax
      DO k = 1,kmax
         
         IF (DTHDt(i,j,k) .ge. 10e-5) THEN
         !IF (DTHDt(i,j,k) .ne. 0) THEN
            
            num_active_points = num_active_points + 1
            active_points(num_active_points,1) = i
            active_points(num_active_points,2) = j
            active_points(num_active_points,3) = k
            active_points(num_active_points,4) = KCBL(i,j)
            active_points(num_active_points,5) = KCBL(i,j) - k

            exit

         endIF
      endDO
   endDO
endDO

DO i = 1,num_active_points
   !print*, i, active_points(i,1), active_points(i,2), active_points(i,3), active_points(i,4), active_points(i,5)
endDO
!print*, 'Number of Active Points = ', num_active_points


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   MAKE PERTURBATIONS FOR ACTIVE POINTS AND COMPUTE JACOBIAN   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

JACOBIAN = 0.0
REVALS_MAX = 0.0
IEVALS_MAX = 0.0
RMSE = 0.0

DO pert_point = 1,1!num_active_points
!DO pert_point = 3153,3153 !num_active_points
!DO pert_point = 2011,2011!num_active_points

    !print*, active_points(pert_point,1), active_points(pert_point,2), active_points(pert_point,3), active_points(pert_point,4), active_points(pert_point,5)
    !print*, 'Lon = ', lon(icol), 'Lat = ', lat(jcol)

    !print*, pert_point

    icol = active_points(pert_point,1)
    jcol = active_points(pert_point,2)
   
   DO ii = 1,14 !Perturb different amounts (1 - 7)
   !DO ii = 1,1!10,9!14 !Perturb different amounts (1 - 7)

      pert = pert_vec(ii)
 
      !print*, 'Pert Size = ', pert
         
      DO kk = 1,4 !PERTURB EACH VARIABLE

         DO ll = 1,kmax !Perturb level by level
         !DO ll = 45,45 !Perturb level by level
 
            pert_amount_var = 0.0
               
            IF (kk == 1) THEN !Perturb U

               U_pert = U_in(icol,jcol,:)
               U_pert(ll) = U_pert(ll) + pert*U_sd(ll)
               pert_amount_var = U_pert(ll) - U_in(icol,jcol,ll)
                     
               V_pert  = V_in(icol,jcol,:)
               TH_pert = TH_in(icol,jcol,:)
               Q_pert  = Q_in(icol,jcol,:)
               P_pert  = P_in(icol,jcol,:)
                       
            elseIF (kk == 2) THEN !Perturb V

               V_pert = V_in(icol,jcol,:)
               V_pert(ll) = V_pert(ll) + pert*V_sd(ll)
               pert_amount_var = V_pert(ll) - V_in(icol,jcol,ll)
                
               U_pert  = U_in(icol,jcol,:)
               TH_pert = TH_in(icol,jcol,:)
               Q_pert  = Q_in(icol,jcol,:)
               P_pert  = P_in(icol,jcol,:)
                       
            elseIF (kk == 3) THEN !Perturb TH

               TH_pert = TH_in(icol,jcol,:)
               TH_pert(ll) = TH_pert(ll) + pert*TH_sd(ll)
               pert_amount_var = TH_pert(ll) - TH_in(icol,jcol,ll)
                 
               U_pert = U_in(icol,jcol,:)
               V_pert = V_in(icol,jcol,:)
               Q_pert = Q_in(icol,jcol,:)
               P_pert = P_in(icol,jcol,:)

            elseIF (kk == 4) THEN !Perturb Q

               Q_pert = Q_in(icol,jcol,:)
               Q_pert(ll) = Q_pert(ll) + pert*Q_sd(ll)
               pert_amount_var = Q_pert(ll) - Q_in(icol,jcol,ll)

               U_pert  = U_in(icol,jcol,:)
               V_pert  = V_in(icol,jcol,:)
               TH_pert = TH_in(icol,jcol,:)
               P_pert  = P_in(icol,jcol,:)

            elseIF (kk == 5) THEN !Perturb P

               P_pert = P_in(icol,jcol,:)
               P_pert(k) = P_pert(ll) + pert*P_sd(ll)
               pert_amount_var = P_pert(ll) - P_in(icol,jcol,ll)

               U_pert  = U_in(icol,jcol,:)
               V_pert  = V_in(icol,jcol,:)
               TH_pert = TH_in(icol,jcol,:)
               Q_pert  = Q_in(icol,jcol,:)

            endIF

            !SAVE PRECALL PROGNOSTIC VARIABLES
            U_pert_in = U_pert
            V_pert_in = V_pert
            TH_pert_in = TH_pert
            Q_pert_in = Q_pert
            P_pert_in = P_pert


            call convection_driver(kmax,DT,KCBL(icol,jcol),Ts(icol,jcol),FRLAND(icol,jcol),PREF,&
                                   U_pert,V_pert,TH_pert,Q_pert,P_pert,&
                                   PRC3_RAS_pert,PRC3_LSC_pert,RASPARAMS,TH_in(icol,jcol,:),Q_in(icol,jcol,:))

               
            !COMPUTE VELOCITY CHANGE, HEATING AND MOISTENING RATES
            DUDt_pert  = (U_pert - U_pert_in)/DT
            DVDt_pert  = (V_pert - V_pert_in)/DT
            DTHDt_pert = (TH_pert - TH_pert_in)/DT
            DQDt_pert  = (Q_pert - Q_pert_in)/DT
!print*,     DTHDt_pert -  DTHdt(icol,jcol,:)
            JACOBIAN(ii,0*kmax+1:1*kmax,(kk-1)*kmax+ll) = (DUDt_pert - DUdt(icol,jcol,:))/pert_amount_var
            JACOBIAN(ii,1*kmax+1:2*kmax,(kk-1)*kmax+ll) = (DVDt_pert - DVdt(icol,jcol,:))/pert_amount_var
            JACOBIAN(ii,2*kmax+1:3*kmax,(kk-1)*kmax+ll) = (DTHDt_pert - DTHdt(icol,jcol,:))/pert_amount_var
            JACOBIAN(ii,3*kmax+1:4*kmax,(kk-1)*kmax+ll) = (DQDt_pert - DQdt(icol,jcol,:))/pert_amount_var

         endDO !ll
      endDo !kk         
         
   endDO !ii, Iteration over different perturbation amounts

   !Compute the eigenvalues of the Jacobian
   IF (COMP_EVALS == .true.) THEN

      DO k = 1,pert_vec_size

         EYE = 0.0
         A = 0.0
         DO j = 1,5*kmax
            EYE(j,j) = 1.0
         endDO
         
         !Convert to forward operator.     
         A = EYE + DT*JACOBIAN(k,1:4*kmax,1:4*kmax)

         LWORK = -1
         work = 0.0
         info = 0

         CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
   
         LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

         !Solve eigenproblem.
         CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

         IF ( INFO.GT.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         END IF

         REVALS_MAX(pert_point,k) = maxval(WR)
         IEVALS_MAX(pert_point,k) = maxval(WI)

      endDO
 
      REVALS_MAX(pert_point,pert_vec_size+1) = maxval(REVALS_MAX(pert_point,1:pert_vec_size))
      IEVALS_MAX(pert_point,pert_vec_size+1) = maxval(IEVALS_MAX(pert_point,1:pert_vec_size))
 
      print*, 'Max Eigenvalue = ', REVALS_MAX(pert_point,:)

   endIF

   !Compute RMS error between Jacobians
   IF (COMP_RMSE == .true.) THEN

      DO k = 2,pert_vec_size

         sum2 = 0.
         sumn = 0.

         DO i = 2*kmax+1,4*kmax
            DO j = 2*kmax+1,4*kmax

               IF ( (JACOBIAN(1,i,j) .ne. 0) .or. (JACOBIAN(k,i,j) .ne. 0) ) Then

                  sum2 = sum2 + (JACOBIAN(1,i,j) - JACOBIAN(k,i,j))**2
                  sumn = sumn+1.

               endIF

            endDO
         endDO

         RMSE(pert_point,k-1) = sqrt(sum2/sumn)

      enddo

      RMSE(pert_point,pert_vec_size) = maxval(RMSE(pert_point,1:pert_vec_size-1))

      !print*, RMSE(pert_point,:)

   ENDIF


endDO ! Over active points.

CALL ChDir('/home/drholdaw/Lin_Moist_Physics/RAS_Jacobian_Paper/mfiles')

! Write Data to nc4 file
IF (WRITE_DATA == .true.) THEN
   output_file = 'JACOBIAN_DEEP.nc4'
   rcode = nf90_create(trim(adjustl(output_file)), NF90_CLOBBER, ncid)
   rcode = nf90_def_dim(ncid, "x", kmax*5, x_dimid) 
   rcode = nf90_def_dim(ncid, "y", kmax*5, y_dimid)
   rcode = nf90_def_dim(ncid, "z", pert_vec_size, z_dimid)
   rcode = nf90_def_var(ncid, "J", NF90_DOUBLE, (/ z_dimid, y_dimid, x_dimid /), varid)
   rcode = nf90_def_var(ncid, "KCBL", NF90_DOUBLE, varid1)
   rcode = nf90_def_var(ncid, "TOP", NF90_DOUBLE, varid2)
   rcode = nf90_def_var(ncid, "DEPTH", NF90_DOUBLE, varid3)
   rcode = nf90_enddef(ncid)
   rcode = nf90_put_var(ncid, varid, JACOBIAN)
   rcode = nf90_put_var(ncid, varid1, active_points(1,4))
   rcode = nf90_put_var(ncid, varid2, active_points(1,3))
   rcode = nf90_put_var(ncid, varid3, active_points(1,5))
   rcode = nf90_close(ncid)
endIF


IF (COMP_EVALS == .true.) THEN
   output_file = 'EVALS.nc4'
   rcode = nf90_create(trim(adjustl(output_file)), NF90_CLOBBER, ncid)
   rcode = nf90_def_dim(ncid, "x", pert_vec_size+1, x_dimid) 
   rcode = nf90_def_dim(ncid, "z", 6, z_dimid) 
   rcode = nf90_def_dim(ncid, "y", num_active_points, y_dimid)
   rcode = nf90_def_var(ncid, "mu", NF90_DOUBLE, (/ y_dimid, x_dimid /), varid)
   rcode = nf90_def_var(ncid, "Points", NF90_DOUBLE, (/ y_dimid, z_dimid /), varid1)   
   rcode = nf90_enddef(ncid)
   rcode = nf90_put_var(ncid, varid, REVALS_MAX(1:num_active_points,:))
   rcode = nf90_put_var(ncid, varid1, active_points(1:num_active_points,:))
   rcode = nf90_close(ncid)
endIF

IF (COMP_RMSE == .true.) THEN
   output_file = 'RMSE.nc4'
   rcode = nf90_create(trim(adjustl(output_file)), NF90_CLOBBER, ncid)
   rcode = nf90_def_dim(ncid, "x", pert_vec_size, x_dimid) 
   rcode = nf90_def_dim(ncid, "z", 6, z_dimid) 
   rcode = nf90_def_dim(ncid, "y", num_active_points, y_dimid)
   rcode = nf90_def_var(ncid, "RMSE", NF90_DOUBLE, (/ y_dimid, x_dimid /), varid)
   rcode = nf90_def_var(ncid, "Points", NF90_DOUBLE, (/ y_dimid, z_dimid /), varid1)   
   rcode = nf90_enddef(ncid)
   rcode = nf90_put_var(ncid, varid, RMSE(1:num_active_points,:))
   rcode = nf90_put_var(ncid, varid1, active_points(1:num_active_points,:))
   rcode = nf90_close(ncid)
endIF


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





