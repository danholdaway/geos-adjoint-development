subroutine global_driver(imax,jmax,kmax,icol_min,jcol_min,icol_max,jcol_max,DT,U,V,TH,Q,Ps,KCBL,Ts,FRLAND,ak,bk,PREF,&
                   U_OUT,V_OUT,TH_OUT,Q_OUT,ls_precip,Ph)

IMPLICIT NONE

!MAPL CONSTANTS
 real, parameter :: MAPL_PI_R8  = 3.14159265358979323846
 real, parameter :: MAPL_PI     = MAPL_PI_R8
 real, parameter :: MAPL_GRAV   = 9.80                   ! m^2/s
 real, parameter :: MAPL_RADIUS = 6376.0E3               ! m
 real, parameter :: MAPL_OMEGA  = 2.0*MAPL_PI/86164.0    ! 1/s
 real, parameter :: MAPL_ALHL   = 2.4665E6               ! J/kg @15C
 real, parameter :: MAPL_ALHF   = 3.3370E5               ! J/kg
 real, parameter :: MAPL_ALHS   = MAPL_ALHL+MAPL_ALHF    ! J/kg
 real, parameter :: MAPL_STFBOL = 5.6734E-8              ! W/(m^2 K^4)
 real, parameter :: MAPL_AIRMW  = 28.97                  ! kg/Kmole
 real, parameter :: MAPL_H2OMW  = 18.01                  ! kg/Kmole
 real, parameter :: MAPL_O3MW   = 47.9982                ! kg/Kmole
 real, parameter :: MAPL_RUNIV  = 8314.3                 ! J/(Kmole K)
 real, parameter :: MAPL_KAPPA  = 2.0/7.0                ! --
 real, parameter :: MAPL_RVAP   = MAPL_RUNIV/MAPL_H2OMW  ! J/(kg K)
 real, parameter :: MAPL_RGAS   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
 real, parameter :: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K)
 real, parameter :: MAPL_P00    = 100000.0               ! Pa
 real, parameter :: MAPL_CAPICE = 2000.                  ! J/(K kg)
 real, parameter :: MAPL_CAPWTR = 4218.                  ! J/(K kg)
 real, parameter :: MAPL_RHOWTR = 1000.                  ! kg/m^3
 real, parameter :: MAPL_NUAIR  = 1.533E-5               ! m^2/S (@ 18C)
 real, parameter :: MAPL_TICE   = 273.16                 ! K
 real, parameter :: MAPL_SRFPRS = 98470                  ! Pa
 real, parameter :: MAPL_KARMAN = 0.40                   ! --
 real, parameter :: MAPL_USMIN  = 1.00                   ! m/s
 real, parameter :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0 
 real, parameter :: MAPL_AVOGAD = 6.023E26
 real, parameter :: PMIN_DET    = 3000.0
 real, parameter :: RKAP        = MAPL_RGAS/MAPL_CP
 real, parameter :: PMIN_CBL    = 50000.0                !Constants for BL bit

!!INPUTS!!
INTEGER :: imax, jmax, kmax
REAL :: DT

REAL, DIMENSION(imax,jmax,kmax) :: U, V, TH, Q
REAL, DIMENSION(imax,jmax) :: Ps, Ts, FRLAND
INTEGER, DIMENSION(imax,jmax) :: KCBL

REAL, DIMENSION(kmax+1) :: ak, bk, PREF

INTEGER :: icol_min, jcol_min
INTEGER :: icol_max, jcol_max

!!OUTPUTS!!
REAL, DIMENSION(imax,jmax,kmax) :: U_RASOUT, V_RASOUT, TH_RASOUT, Q_RASOUT
REAL, DIMENSION(imax,jmax,kmax) :: U_LSCOUT, V_LSCOUT, TH_LSCOUT, Q_LSCOUT
REAL, DIMENSION(imax,jmax,kmax) :: U_OUT, V_OUT, TH_OUT, Q_OUT

!!LOCALS!!
INTEGER :: i, j, k, kx, iras, jras, ilsc, jlsc, l

INTEGER, PARAMETER :: IDIM = 1, IRUN = 1

INTEGER :: N_DTL, ICMIN, SEEDRAS(2)
REAL :: SIGE(kmax+1)

INTEGER :: KCBL_ras(IDIM)

!DERIVED MODEL VARIABLES
REAL :: P(imax,jmax,kmax+1), Ph(imax,jmax,kmax), Pi(imax,jmax,kmax+1), Pih(imax,jmax,kmax), TEMP(imax,jmax,kmax)
REAL :: LONS(imax), LATS(jmax) 
INTEGER :: IRAS1(imax), JRAS1(jmax), IRAS_ras(IDIM), JRAS_ras(IDIM)


!!!!!!!!!!!
!!! RAS !!!
!!!!!!!!!!!

!MODEL VARIABLES READY FOR RAS - PICKS OUT COLUMN AND ALLOWS FOR TRUNCATION ABOVE CONVECTION THROUGH kmax
REAL :: T_ras(IDIM,kmax), Q_ras(IDIM,kmax), TH_ras(IDIM,kmax), U_ras(IDIM,kmax), V_ras(IDIM,kmax)
REAL :: Ps_ras(IDIM), P_ras(IDIM,kmax+1), Pi_ras(IDIM,kmax+1), Ph_ras(IDIM,kmax), Pih_ras(IDIM,kmax)
REAL :: ZLE_ras(IDIM,0:kmax), ZLO_ras(IDIM,kmax), ZCBLx_ras(IDIM)
REAL :: QSS_ras(IDIM,kmax), DQS_ras(IDIM,kmax)

!BOUNDARY LAYER TERMS
REAL :: WGT0(IDIM,kmax), WGT1(IDIM,kmax), MXDIAMx(IDIM), TPERT_ras(IDIM), QPERT_ras(IDIM)
REAL, PARAMETER :: CBL_TPERT = 1.0, CBL_QPERT = 0.0, CBL_TPERT_MXOCN = 2.0, CBL_TPERT_MXLND = 4.0
REAL :: QSSFC(imax,jmax), QSSFC_ras(IDIM)

!RAS-2 VARIABLES, NOT CURRENTLY REQUIRED
REAL :: GZLE(IDIM,kmax+1), GZLO(IDIM,kmax), QLCN(IDIM,kmax), QICN(IDIM,kmax)

!COLUMN OUTPUTS - AS RETURNED FROM RAS
REAL :: CNV_DQLDT(IDIM,kmax), CNV_MF0(IDIM,kmax), CNV_MFD(IDIM,kmax), CNV_MFC(IDIM,kmax)
REAL :: CNV_PRC3(IDIM,kmax), CNV_UPDF(IDIM,kmax), CNV_CVW(IDIM,kmax), CNV_QC(IDIM,kmax)
REAL :: CLCN(IDIM,kmax), HHO(IDIM,kmax), HSO(IDIM,kmax), RASPRCP(IDIM)

!OTHER PARAMETERS
INTEGER, PARAMETER :: n_ras_params = 25
REAL :: RASPARAMS(n_ras_params)
INTEGER, PARAMETER :: ITRCR = 1


!!!!!!!!!!!
!!! LSC !!!
!!!!!!!!!!!

!MODEL VARIABLES READY FOR LSC 
REAL :: TH_lsc(kmax), Q_lsc(kmax)
REAL :: Ps_lsc, P_lsc(kmax+1), Ph_lsc(kmax), Pih_lsc(kmax), lsprecip_lsc

REAL, DIMENSION(imax,jmax) :: ls_precip

!!! BEGIN CALCS !!!

SIGE = PREF/PREF(kmax+1) ! this should eventually change (Comment in Moist code)
ICMIN = max(1,count(PREF < PMIN_DET))

!DERIVE SOME MODEL VARIABLES
DO k = 1,kmax+1
   P(:,:,k) = ak(k)*.01 + bk(k)*Ps
endDO
Ph(:,:,:) = 0.5*(P(:,:,1:kmax) + P(:,:,2:kmax+1))
!Exner Pressures
Pi(:,:,:) = (P(:,:,:)/1000.0)**Rkap
Pih(:,:,:) = (Ph(:,:,:)/1000.0)**Rkap

TEMP(:,:,:) = TH(:,:,:)*Pih(:,:,:)

DO i = 1,imax
   lons(i) = -180.0+(i-1)*(360.0/imax)
endDO
Do j = 1,jmax
   lats(j) = -90.0 +(j-1)*(180.0/jmax)
endDO
IRAS1 = nint(LONS)
JRAS1 = nint(LATS)

!Neither Ts or Ps
QSSFC = 0.0

DO jras = jcol_min,jcol_max

   DO iras = icol_min,icol_max
   
      IRAS_ras(IDIM) = IRAS1(iras)
      JRAS_ras(IDIM) = JRAS1(jras)

      MXDIAMx(IDIM) = 0.0

      !Pick out sounding for RAS
      Q_ras(IDIM,:) = Q(iras,jras,:)
      T_ras(IDIM,:) = TEMP(iras,jras,:)
      TH_ras(IDIM,:) = TH(iras,jras,:)
      U_ras(IDIM,:) = U(iras,jras,:)
      V_ras(IDIM,:) = V(iras,jras,:)

      Ps_ras(IDIM) = Ps(iras,jras)
      P_ras(IDIM,:) = ak(:)*.01 + bk(:)*Ps_ras(IDIM)
      Ph_ras(IDIM,1:kmax) = 0.5*( P_ras(IDIM,1:kmax) +  P_ras(IDIM,2:kmax+1) ) 

      Pi_ras(IDIM,:) = (P_ras(IDIM,:)/1000.0)**Rkap
      Pih_ras(IDIM,:) = (Ph_ras(IDIM,:)/1000.0)**Rkap

      !COMPUTE QSS AND DQS - Depend on perturbed quantities.
      call MRQSAT(T_ras(IDIM,:), Ph_ras(IDIM,:), QSS_ras(IDIM,:), DQS_ras(IDIM,:),1,1,kmax)

      !SEEDRAS TERMS
      SEEDRAS(1) = 1000000 * ( 100*T_ras(IDIM,kmax)   - INT( 100*T_ras(IDIM,kmax) ) )
      SEEDRAS(2) = 1000000 * ( 100*T_ras(IDIM,kmax-1) - INT( 100*T_ras(IDIM,kmax-1) ) )
      
      !BOUNDARY LAYER TERMS ETC
      ZLE_ras(IDIM,:) = 0.0
      DO k=kmax,1,-1
         ZLE_ras(IDIM,k-1) = TH_ras(IDIM,k) * (1.+MAPL_VIREPS*Q_ras(IDIM,k))
         ZLO_ras(IDIM,k  ) = ZLE_ras(IDIM,k) + (MAPL_CP/MAPL_GRAV)*( Pi_ras(IDIM,k)-Pih_ras (IDIM,k  ) ) * ZLE_ras(IDIM,k-1)
         ZLE_ras(IDIM,k-1) = ZLO_ras(IDIM,k) + (MAPL_CP/MAPL_GRAV)*( Pih_ras (IDIM,k)-Pi_ras(IDIM,k-1) ) * ZLE_ras(IDIM,k-1)
      endDO

      KCBL_ras(IDIM) = KCBL(iras,jras)

      WGT0 = 0.0 
      WGT1 = 0.0

      !Case 6
      WGT0(IDIM, KCBL_ras(IDIM):kmax ) = 1.0 
      WGT1(IDIM, KCBL_ras(IDIM):kmax ) = 1.0       

      ZCBLx_ras = ZLE_ras(IDIM, KCBL_ras(IDIM)-1 )

      ! Cheat by adding a kick to CB temp and q
      QSSFC_ras(IDIM) = QSSFC(imax,jmax)

      TPERT_ras(IDIM)  = CBL_TPERT * ( Ts(iras,jras) - ( T_ras(IDIM,kmax)+ MAPL_GRAV*ZLO_ras(IDIM,kmax)/MAPL_CP )  ) 
      QPERT_ras(IDIM)  = CBL_QPERT * ( QSSFC_ras(IDIM) - Q_ras(IDIM,kmax) )          
      TPERT_ras(IDIM)  = MAX( TPERT_ras(IDIM) , 0.0 )
      QPERT_ras(IDIM)  = MAX( QPERT_ras(IDIM) , 0.0 )

      if (FRLAND(iras,jras) < 0.1) then
         TPERT_ras(IDIM)  = MIN( TPERT_ras(IDIM) , CBL_TPERT_MXOCN ) ! ocean
      else
         TPERT_ras(IDIM)  = MIN( TPERT_ras(IDIM) , CBL_TPERT_MXLND ) ! land
      endif
      
      !RASPARAMS
      RASPARAMS( 1) = 1.000
      RASPARAMS( 2) = 0.05
      if (FRLAND(I,J)<0.1) then
         RASPARAMS( 3) = 2.5e-3   ! ocean value
      else
         RASPARAMS( 3) = 2.5e-3  ! land value
      endif
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
      if (imax==144) then
         RASPARAMS(23) = 4000.0   !Low Res
      elseif (imax==576) then
         RASPARAMS(23) = 1000.0   !High Res
      endif
      RASPARAMS(24) = 0.5
      RASPARAMS(25) = 0.65

      N_DTL = KCBL(iras,jras) - ICMIN

      !IF ( (iras == 109) .AND. (jras == 41) ) THEN

         !write(*,*) 'I = ', iras
         !write(*,*) 'J = ', jras

         !write(*,*) 'IDIM = ', IDIM
         !write(*,*) 'IRIN = ', IRUN
         !write(*,*) 'kmax = ', kmax
         !write(*,*) 'ICMIN = ', ICMIN
         !write(*,*) 'DT = ', DT
         !write(*,*) 'CP = ', MAPL_CP
         !write(*,*) 'ALPH = ', MAPL_ALHL
         !write(*,*) 'ALPHS = ', MAPL_ALHS
         !write(*,*) 'TICE = ', MAPL_TICE, TEMV4_ALL
         !write(*,*) 'GRAV = ', MAPL_GRAV
         !write(*,*) 'SEEDRAS = ', SEEDRAS
         !write(*,*) 'IRAS1 = ', IRAS_ras
         !write(*,*) 'JRAS1 = ', JRAS_ras
         !write(*,*) 'SIGE = ', SIGE

         !write(*,*) 'KCBL_ras = ', KCBL_ras
         !write(*,*) 'WGT0 = ', WGT0
         !write(*,*) 'WGT1 = ', WGT1
         !write(*,*) 'ZCBLx_ras = ', ZCBLx_ras
         !write(*,*) 'MXDIAMx = ', MXDIAMx
         !write(*,*) 'TPERT_ras = ', TPERT_ras
         !write(*,*) 'QPERT_ras = ', QPERT_ras
    
         !write(*,*) 'TH = ', TH_ras
         !write(*,*) 'Q = ', Q_ras
         !write(*,*) 'U = ', U_ras
         !write(*,*) 'V = ', V_ras
         !write(*,*) 'QSS = ', QSS_ras
         !write(*,*) 'DQS = ', DQS_ras

         !write(*,*) 'P_ras = ', P_ras
         !write(*,*) 'Pi_ras = ', Pi_ras

         !write(*,*) 'RASPARAMS = ', RASPARAMS
         !write(*,*) 'ITRCR = ', ITRCR

      !endIF 

       
!call RASE(IDIM, IRUN, kmax, ICMIN, DT, MAPL_CP, MAPL_ALHL, MAPL_ALHS, MAPL_TICE, MAPL_GRAV, SEEDRAS, IRAS_ras, JRAS_ras, SIGE, KCBL_ras, WGT0, WGT1, ZCBLx_ras, MXDIAMx, TPERT_ras, QPERT_ras, TH_ras, Q_ras, U_ras, V_ras, QSS_ras, DQS_ras, Pih_ras, Ph_ras, GZLO, GZLE, QLCN, QICN, P_ras, Pi_ras, CNV_DQLDT, CNV_MF0, CNV_MFD, CNV_MFC, CNV_PRC3, CNV_UPDF, CNV_CVW, CNV_QC, CLCN, HHO, HSO, RASPRCP, RASPARAMS, ITRCR)

call rase( IDIM,KMAX,ICMIN,DT,MAPL_CP,MAPL_ALHL,MAPL_GRAV,SEEDRAS,SIGE,N_DTL,WGT0(IDIM,:),WGT1(IDIM,:),TPERT_ras(IDIM),QPERT_ras(IDIM),TH_ras(IDIM,:),Q_ras(IDIM,:),U_ras(IDIM,:),V_ras(IDIM,:),QSS_ras(IDIM,:), DQS_ras(IDIM,:),P_ras(IDIM,:), Pi_ras(IDIM,:), CNV_DQLDT(IDIM,:), CNV_MF0(IDIM,:), CNV_MFD(IDIM,:), CNV_MFC(IDIM,:),CNV_PRC3(IDIM,:), CNV_UPDF(IDIM,:), CNV_CVW(IDIM,:), CNV_QC(IDIM,:),HHO(IDIM,:), HSO(IDIM,:), RASPARAMS )

          !call rase( IDIM, KMAX, ICMIN, DT, MAPL_CP, MAPL_ALHL, MAPL_GRAV, SEEDRAS, SIGE, N_DTL, &
          !           WGT0(IDIM,:), WGT1(IDIM,:), TPERT_ras(IDIM), QPERT_ras(IDIM), &
          !           TH_ras(IDIM,:), Q_ras(IDIM,:), U_ras(IDIM,:), V_ras(IDIM,:), &
          !           QSS_ras(IDIM,:), DQS_ras(IDIM,:),&
          !           P_ras(IDIM,:), Pi_ras(IDIM,:), &
          !           CNV_DQLDT(IDIM,:), CNV_MF0(IDIM,:), CNV_MFD(IDIM,:), CNV_MFC(IDIM,:), &
          !           CNV_PRC3(IDIM,:), CNV_UPDF(IDIM,:), CNV_CVW(IDIM,:), CNV_QC(IDIM,:), &
          !           HHO(IDIM,:), HSO(IDIM,:), RASPARAMS )


      T_ras = TH_ras*Pih_ras

      U_rasout(iras,jras,:) = U_ras(IDIM,:)
      V_rasout(iras,jras,:) = V_ras(IDIM,:)
      TH_rasout(iras,jras,:) = TH_ras(IDIM,:)
      Q_rasout(iras,jras,:) = Q_ras(IDIM,:)
      
      
   endDO !iras loop

endDO !jras loop
!Update temperature
!U(:,:,:) = U_rasout(:,:,:)
!V(:,:,:) = V_rasout(:,:,:)
!TH(:,:,:) = TH_rasout(:,:,:)
!Q(:,:,:) = Q_rasout(:,:,:)
!TEMP(:,:,:) = TH_rasout(:,:,:)*Pih(:,:,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  LARGE SCALE CLOUD LOOP !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO jlsc = jcol_min,jcol_max
   DO ilsc = icol_min,icol_max
!DO jlsc = 83,83!jcol_min,jcol_max
!   DO ilsc = 128,128!

      TH_lsc(:) = TH(ilsc,jlsc,:)
      Q_lsc(:) = Q(ilsc,jlsc,:)
      
      Ps_lsc = Ps(ilsc,jlsc)
      P_lsc(:) = ak(:)*.01 + bk(:)*Ps_lsc
      Ph_lsc(1:kmax) = 0.5*( P_lsc(1:kmax) +  P_lsc(2:kmax+1) ) 
      Pih_lsc(:) = (Ph_lsc(:)/1000.0)**Rkap
      
      call large_cloud(kmax,TH_lsc,Q_lsc,P_lsc,Ph_lsc,Pih_lsc,lsprecip_lsc)

      TH_lscout(ilsc,jlsc,:) = TH_lsc(:)
      Q_lscout(ilsc,jlsc,:) = Q_lsc(:)
      ls_precip(ilsc,jlsc) = lsprecip_lsc

   endDO
endDO

TH(:,:,:) = TH_lscout(:,:,:)
Q(:,:,:) = Q_lscout(:,:,:)
!Q(ilsc-1,jlsc-1,:) = Q_lscout(ilsc-1,jlsc-1,:)
TEMP(:,:,:) = TH_lscout(:,:,:)*Pih(:,:,:)


!UPDATE _OUT TO SEND BACK
U_out = U
V_out = V
TH_out = TH
Q_out = Q

end subroutine global_driver



subroutine large_cloud(kmax,th,q,P,Ph,Pih,precip_surf)

IMPLICIT NONE

!INPUT/OUTPUT
INTEGER :: kmax
REAL, DIMENSION(kmax) :: th,q,ph,pih
REAL, DIMENSION(kmax+1) :: P

!LOCALS
INTEGER :: k, kk

REAL, PARAMETER :: ALHL = 2.4665E6, em = 0.887, grav = 9.80665, cp = 1005.7

REAL, DIMENSION(kmax) :: qs, cpm, dqsdT, T
REAL, DIMENSION(kmax) :: deltaq, deltaT

REAL :: precip_surf, qsatvp(10000)
REAL, DIMENSION(kmax) :: precipz


precipz = 0.0
deltaq = 0.0
deltaT = 0.0

!Constants
cpm = cp*(1+em*q)

!Compute temperature for the column
T = th*Pih

!COMPUTE qs AND dqsdT - Depend on perturbed quantities.
call MRQSAT(T, Ph, qs, DqsdT,1,1,kmax)


!dqsdt = 0.0
!call qtable(qsatvp)
!call MRQSAT_RON(T, Ph, qs, DqsdT, qsatvp)

DO k = 40,kmax

   if (q(k) - qs(k) .le. 0) THEN

      deltaq(k) = 0.0

   elseif (q(k) - qs(k) .gt. 0) THEN

      deltaq(k) = (q(k)-qs(k))/(1+dqsdT(k)*ALHL/cpm(k))

   endif

endDo

deltaT = (ALHL*deltaq)/cpm

!Update temeprature and specific humidity
T = T + deltaT
q = q - deltaQ

!Convect back to potential temp.
th = T/Pih

!COLUMN RAIN
DO k = 1,kmax
   precipz(k) = 100*deltaq(k)*(P(k+1)-P(k))/grav
endDO

!SURFACE RAIN
precip_surf = 0.0
DO k = 1,kmax
   precip_surf = precip_surf + precipz(k)
endDO

end subroutine large_cloud


      subroutine qtable (qsatvp)

!     Construct tables for values for the saturation vapor pressure 

      real :: svpt0=273.15, svp1=.6112, svp2=17.67, svp3=29.65
      integer, parameter :: tmin=130., tmax=330., nqsatvp=10000
      real qsatvp(nqsatvp)

      qsatvp(1)=tmax
      qsatvp(2)=tmin
      tr=(tmax-tmin)/real(nqsatvp-4)
      qsatvp(3)=tr
      do n=4,nqsatvp 
         t=tmin+(n-4)*tr
         a=svp2*(t-svpt0)/(t-svp3)
         qsatvp(n)=svp1*exp(a)
      enddo 
      
      end subroutine qtable


SUBROUTINE MRQSAT_RON(TT,P,Q,DQDT,qsatvp)

!  compute saturation specific humidity and its derivative
!  with repect to temperature. A modification of RAS software to
!  use MAMS1 table lookup.
!implicit none

      parameter ( airmw  = 28.97      )
      parameter ( h2omw  = 18.01      )
      parameter ( one    = 1.0        )
      parameter (esfac = h2omw/airmw        )
      parameter (erfac = (one-esfac)/esfac  )
      real tt(72), p(72), q(72), dqdt(72), tx(72), D(72), v(72), dqdta(72)
      integer ip(72)
      parameter (ep2=0.622)
      integer, parameter :: nqsatvp=10000
      real qsatvp(nqsatvp)
      
      rtr=1./qsatvp(3)

!      do i=1,lenc
         tx=(tt-qsatvp(2))*rtr  
         ip=int(tx)
         v=qsatvp(ip+4)+(tx-ip)*(qsatvp(ip+5)-qsatvp(ip+4))
         d=1./(.1*p-ep2*erfac*v)

!  tx is t-tmin in units of qsatvp(2)
!  factor of 0.1 before p changes from mb to cb (v in units of cb)
   
         q=ep2*v*d

         !if (ldqdt)  then
           dqdta = (qsatvp(ip+5)-qsatvp(ip+4))*rtr
           dqdt = ep2*d*dqdta*(1. + erfac*q)
         !endif
!      enddo


end subroutine MRQSAT_RON

SUBROUTINE MRQSAT(TT,P,Qs,DQsDT,imax,jmax,kmax)

      implicit none

      !Inputs
      integer :: imax, jmax, kmax
      real, dimension(imax,jmax,kmax) :: TT, P, Qs, DQsDT
      logical :: LDQDT

      !Constants
      real, parameter :: svpt0=273.15, svp1=.6112, svp2=17.67, svp3=29.65, airmw = 28.97, h2omw = 18.01, one = 1.0, ep2 = 0.622
      real, parameter :: esfac = h2omw/airmw, erfac = (one-esfac)/esfac

      !Internals
      real :: dqdta, qsatvp, d
      integer :: ii, jj, kk      

      do ii = 1,imax
         do jj = 1,jmax
            do kk=1,kmax

            qsatvp = svp1*exp(svp2*(TT(ii,jj,kk) - svpt0)/(TT(ii,jj,kk) - svp3))
            d = 1./(.1*p(ii,jj,kk)-ep2*erfac*qsatvp)
            Qs(ii,jj,kk)=ep2*qsatvp*d

            dqdta =  (svp2*(svpt0-svp3)/(svp3-TT(ii,jj,kk))**2 )*qsatvp
            dQsdt(ii,jj,kk) = ep2*d*dqdta*(1. + erfac*Qs(ii,jj,kk))
                 
            enddo
         enddo
      enddo

end subroutine MRQSAT

SUBROUTINE RASE(IDIM,K0,ICMIN,DT,CPO,ALHLO,GRAVO,SEEDRAS,SIGE,N_DTL,WGT0,WGT1,TPERT,QPERT,&
THO,QHO,UHO,VHO,QSS,DQS,PLE,PKE,CLW,FLX,FLXD,FLXC,CNV_PRC3,CNV_UPDFRC,CNV_CVW,CNV_QC,HHO,HSO,RASPARAMS)

IMPLICIT NONE

 !ARGUMENTS
 INTEGER                     ::  IDIM, K0, ICMIN
 REAL, PARAMETER                            ::  MAPL_UNDEF = 1.0e15
 REAL, DIMENSION (IDIM,K0  ) ::  THO, QHO, UHO, VHO
 REAL, DIMENSION (IDIM,K0+1) ::  PLE, PKE
 REAL, DIMENSION (IDIM,K0  ) ::  QSS, DQS
 REAL, DIMENSION (     K0+1) ::  SIGE
 REAL, DIMENSION (IDIM,K0  ) ::  FLX, CLW , FLXD, FLXC
 REAL, DIMENSION (IDIM,K0  ) ::  CNV_PRC3
 REAL, DIMENSION (IDIM,K0  ) ::  CNV_UPDFRC, CNV_QC, CNV_CVW
 REAL, DIMENSION (IDIM,K0  ) ::  HHO, HSO
 REAL                        ::  DT,  CPO, ALHLO, GRAVO

 INTEGER, DIMENSION ( 2 )    ::  SEEDRAS
 INTEGER                     ::  N_DTL !KCBL - 31, Around 20-40
 REAL, DIMENSION (IDIM     ) ::  TPERT,QPERT

 REAL, DIMENSION (K0  )      ::  WGT0,WGT1
 REAL, DIMENSION(25)          ::  RASPARAMS
 
 !LOCALS
 REAL, DIMENSION (IDIM     ) ::  MXDIAM
 REAL,  DIMENSION(K0) :: POI_SV, QOI_SV, UOI_SV, VOI_SV
 REAL,  DIMENSION(K0) :: POI, QOI, UOI, VOI,  DQQ, BET, GAM
 REAL,  DIMENSION(K0) :: POI_c, QOI_c
 REAL,  DIMENSION(K0) :: PRH,  PRI,  GHT, DPT, DPB, PKI, DISSK0,DISSK1
 REAL,  DIMENSION(K0) :: TCU, QCU, CLN, POL,DM
 REAL,  DIMENSION(K0) :: QST, SSL,  RMF, RN1
 REAL,  DIMENSION(K0) :: ETA, EHT,  GM1
 REAL,  DIMENSION(K0) :: HOL, HST, QOL, ZOL
 REAL,  DIMENSION(K0) :: WSP, LAMBDSV
 REAL,  DIMENSION(K0) :: MTKWI
 REAL,  DIMENSION(K0) :: WGHT, WGHT0

 REAL,  DIMENSION(K0+1) :: PRJ, PRS, SHT ,ZET, XYD, XYD0
 REAL, DIMENSION(IDIM,K0) :: LAMBDSV2

 INTEGER,  DIMENSION(K0-1) :: RC
 INTEGER :: K,MY_PE

 REAL UHT, VHT, AKM, ACR, TTH, QQH, SHTRG, WSPBL, DQX
 REAL WFN, TRGEXP, EVP, WLQ, MTKW_MAX
 REAL SHTRG_FAC, SIGE_MINHOL

 INTEGER, PARAMETER :: I = 1
 INTEGER IC, L, ICL , ITR 
 INTEGER NDTLEXPON

 !RASE GLOBAL CONSTANTS
 REAL GRAV, CP, ALHL, CPBG, ALHI, CPI, GRAVI, DDT, LBCP, OBG, AFC 
   
!!!!!!!!!
 REAL CO_AUTO, FRICFAC,DPTH_BL,WUPDRFT,PBLFRAC,AUTORAMPB,CO_ZDEP
 REAL RASAL1, RASAL2, CO_T, RASNCL,FRICLAMBDA,SDQVT1,SDQV2 
 REAL LAMBDA_FAC,STRAPPING,ACRITFAC,HMINTRIGGER,LLDISAGGXP
 REAL LAMBMX_FAC, DIAMMN_MIN,RDTLEXPON, CLI_CRIT,SDQV3, MAXDALLOWED
 REAL RHMN, RHMX
 INTEGER KSTRAP

 Real cld_radius, areal_frac, spect_mflx, cvw_cbase
!!!!!!!!!

 REAL, PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
 REAL, PARAMETER :: RHMAX   = 0.9999

 !LAMBDA LIMITS
 REAL            :: LAMBDA_MIN
 REAL            :: LAMBDA_MAX

 !MAPL CONSTANTS
 real, parameter :: MAPL_PI_R8  = 3.14159265358979323846
 real, parameter :: MAPL_PI     = MAPL_PI_R8
 real, parameter :: MAPL_GRAV   = 9.80                   ! m^2/s
 real, parameter :: MAPL_RADIUS = 6376.0E3               ! m
 real, parameter :: MAPL_OMEGA  = 2.0*MAPL_PI/86164.0    ! 1/s
 real, parameter :: MAPL_ALHL   = 2.4665E6               ! J/kg @15C
 real, parameter :: MAPL_ALHF   = 3.3370E5               ! J/kg
 real, parameter :: MAPL_ALHS   = MAPL_ALHL+MAPL_ALHF    ! J/kg
 real, parameter :: MAPL_STFBOL = 5.6734E-8              ! W/(m^2 K^4)
 real, parameter :: MAPL_AIRMW  = 28.97                  ! kg/Kmole
 real, parameter :: MAPL_H2OMW  = 18.01                  ! kg/Kmole
 real, parameter :: MAPL_O3MW   = 47.9982                ! kg/Kmole
 real, parameter :: MAPL_RUNIV  = 8314.3                 ! J/(Kmole K)
 real, parameter :: MAPL_KAPPA  = 2.0/7.0                ! --
 real, parameter :: MAPL_RVAP   = MAPL_RUNIV/MAPL_H2OMW  ! J/(kg K)
 real, parameter :: MAPL_RGAS   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
 real, parameter :: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K)
 real, parameter :: MAPL_P00    = 100000.0               ! Pa
 real, parameter :: MAPL_CAPICE = 2000.                  ! J/(K kg)
 real, parameter :: MAPL_CAPWTR = 4218.                  ! J/(K kg)
 real, parameter :: MAPL_RHOWTR = 1000.                  ! kg/m^3
 real, parameter :: MAPL_NUAIR  = 1.533E-5               ! m^2/S (@ 18C)
 real, parameter :: MAPL_TICE   = 273.16                 ! K
 real, parameter :: MAPL_SRFPRS = 98470                  ! Pa
 real, parameter :: MAPL_KARMAN = 0.40                   ! --
 real, parameter :: MAPL_USMIN  = 1.00                   ! m/s
 real, parameter :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0 
 real, parameter :: MAPL_AVOGAD = 6.023E26
 real, parameter :: PMIN_DET    = 3000.0
 real, parameter :: RKAP        = MAPL_RGAS/MAPL_CP
 real, parameter :: PMIN_CBL    = 50000.0                !Constants for BL bit

 !*********************************************************************
  
 CNV_PRC3  =0.
 CNV_UPDFRC=0. 
 CNV_CVW   =0. 
 CNV_QC    =0.

 FRICFAC      = RASPARAMS(1)     !  ---  1 (= 1.0)
 SHTRG_FAC    = RASPARAMS(2)     !  ---  2
 CO_AUTO      = RASPARAMS(3)     !  ---  3
 CLI_CRIT     = RASPARAMS(4)     !  ---  4
 RASAL1       = RASPARAMS(5)     !  ---  5
 RASAL2       = RASPARAMS(6)     !  ---  6
 RASNCL       = RASPARAMS(7)     !  ---  7
 LAMBDA_FAC   = RASPARAMS(8)     !  ---  8
 LAMBMX_FAC   = RASPARAMS(9)     !  ---  9
 DIAMMN_MIN   = RASPARAMS(10)    !  --- 10
 FRICLAMBDA   = RASPARAMS(11)    !  --- 11
 RDTLEXPON    = RASPARAMS(12)    !  --- 12
 STRAPPING    = RASPARAMS(13)    !  --- 13 (= -1.0)
 SDQV2        = RASPARAMS(14)    !  --- 14
 SDQV3        = RASPARAMS(15)    !  --- 15
 SDQVT1       = RASPARAMS(16)    !  --- 16
 ACRITFAC     = RASPARAMS(17)    !  --- 17
 HMINTRIGGER  = RASPARAMS(18)    !  --- 18
 LLDISAGGXP   = RASPARAMS(19)    !  --- 19
 PBLFRAC      = RASPARAMS(20)    !  --- 20
 AUTORAMPB    = RASPARAMS(21)    !  --- 21
 CO_ZDEP      = RASPARAMS(22)    !  --- 22
 MAXDALLOWED  = RASPARAMS(23)    !  --- 23
 RHMN         = RASPARAMS(24)    !  --- 24
 RHMX         = RASPARAMS(25)    !  --- 25
   
 WUPDRFT = 2.500

 GRAV  = GRAVO
 ALHL  = ALHLO
 CP    = CPO
 CPI   = 1.0/CP      
 ALHI  = 1.0/ALHL
 GRAVI = 1.0/GRAV
 CPBG  = CP*GRAVI
 DDT   = DAYLEN/DT
 AFC   = -1.04E-4*SQRT(DT*113.84)
 LBCP  = ALHL*CPI
 OBG   = 100.*GRAVI

 HHO = MAPL_UNDEF
 HSO = MAPL_UNDEF
 
 !FINDBASE
 K = N_DTL + ICMIN !Same as KCBL

 call  STRAP0(IDIM,PLE,PKE,THO,QHO,UHO,VHO,QSS,DQS,WGT0,MAPL_RGAS,MAPL_CP,ICMIN,I,K,K0,&
ONEPKAP,SEEDRAS,LBCP,ALHL,MAXDALLOWED,PRS,PRJ,POI,QOI,UOI,VOI,QST,DQQ,POL,PRH,PKI,&
DPT,DPB,PRI,MXDIAM,GHT,GM1,BET,GAM,POI_SV,QOI_SV,UOI_SV,VOI_SV)

 RMF  = 0.0

 LAMBDSV = 0.0
 DISSK0  = 0.0
         

 call HTEST(K,K0,ICMIN,CP,QST,PRJ,PRH,POI,QOI,GRAV,ALHL,CPBG,RHMAX,HOL,HST,SSL,QOL,SHT,ZET,ZOL)

 HHO(I,:) = HOL
 HSO(I,:) = HST

 ! CLOUD LOOP
 call CLOUDE(N_DTL,ICMIN,IDIM,I,K,K0,DT,DDT,DAYLEN,RHMN,RHMX,AUTORAMPB,CO_ZDEP,CP,&
RHMAX,GRAV,ALHL,ALHI,CPBG,LBCP,RASAL1,RASAL2,ACRITFAC,SDQVT1,SDQV2,SDQV3,CO_AUTO,CLI_CRIT,&
PBLFRAC,GRAVI,CPI,FRICFAC,FRICLAMBDA,TPERT,QPERT,MXDIAM,PRJ,PRH,PRS,SIGE,&
DPB,DPT,DQQ,POL,GM1,PKI,BET,GHT,POI,QOI,UOI,VOI,QST,GAM,PRI,RMF)

 IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

    call STRAP1(IDIM, I, K, K0, ICMIN, DDT, DAYLEN, POI, QOI, UOI, VOI, POI_SV, QOI_SV, UOI_SV, VOI_SV, &
                  WGT1, PLE, PRS, THO, QHO, UHO, VHO )
          
 ENDIF 


END SUBROUTINE RASE


!=====================================================================================
SUBROUTINE HTEST(K, K0, ICMIN, CP, QST, PRJ, PRH, POI, QOI, GRAV, ALHL, CPBG, RHMAX, &
                 HOL, HST, SSL, QOL, SHT, ZET, ZOL)
!=====================================================================================

IMPLICIT NONE

!!INPUTS
!Grid constants
INTEGER :: K, K0, ICMIN

!Variables
REAL, DIMENSION(k0) :: QST, POI, QOI
REAL :: PRJ(K0+1), PRH(K0)

!Real constants
REAL :: CP, GRAV, ALHL, CPBG, RHMAX

!!OUTPUTS!!
REAL, DIMENSION(K0) :: HOL, HST, SSL, QOL, ZOL
REAL, DIMENSION(K0+1) :: SHT, ZET

REAL :: TEM

! NB: TEM, SSL, QOL, SHT, ZET, ZOL RECOMPED BY CLOUDE STRAIGHT AWAY

!!LOCALS!!
INTEGER :: L

hol=0. ! HOL initialized here in order not to confuse Valgrind debugger

ZET(K+1) = 0
SHT(K+1) = CP*POI(K)*PRJ(K+1)
DO L=K,ICMIN,-1
   QOL(L) = AMIN1(QST(L)*RHMAX,QOI(L))
   QOL(L) = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
   SSL(L) = CP*PRJ(L+1)*POI(L) + GRAV*ZET(L+1)
   HOL(L) = SSL(L) + QOL(L)*ALHL
   HST(L) = SSL(L) + QST(L)*ALHL
   TEM    = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
   ZET(L) = ZET(L+1) + TEM
   ZOL(L) = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
ENDDO


end subroutine HTEST
!===================================================================================



!======================================================================================================================
SUBROUTINE STRAP0(IDIM, PLE, PKE, THO, QHO, UHO, VHO, QSS, DQS, WGT0, MAPL_RGAS, MAPL_CP, ICMIN, I, K, K0, ONEPKAP, &
                  SEEDRAS, LBCP, ALHL, MAXDALLOWED, &
                  PRS, PRJ, POI, QOI, UOI, VOI, QST, DQQ, POL, PRH, PKI, DPT, DPB, PRI, MXDIAM, GHT, GM1, BET, GAM, &
                  POI_SV, QOI_SV, UOI_SV, VOI_SV &
                 )
!======================================================================================================================

IMPLICIT NONE

!!INPUTS!!
!Grid constants
INTEGER :: IDIM, I, K, K0, ICMIN

!Prognostic Variables, nenamed and updated here.
REAL, DIMENSION(IDIM,K0+1) :: PLE, PKE
REAL, DIMENSION(IDIM,K0)   :: THO, QHO, UHO, VHO, QSS, DQS
REAL, DIMENSION(K0) :: WGT0
!Seedras, used to determine mxdiam.
INTEGER ::  SEEDRAS(2)

!Real constants
REAL :: ALHL, LBCP, ONEPKAP, MAPL_RGAS, MAPL_CP, MAXDALLOWED

!!OUTPUTS!!
REAL, DIMENSION(K0)   :: POI, QOI, UOI, VOI, QST, DQQ, POL, PRH, PKI, DPT, DPB, PRI, GHT, GM1, BET, GAM
REAL, DIMENSION(K0+1) :: PRS, PRJ
REAL, DIMENSION(k0)   :: POI_SV, QOI_SV, UOI_SV, VOI_SV

REAL :: MXDIAM(IDIM)

!!LOCALS!!
REAL, DIMENSION(K0) :: WGHT, MASSF
INTEGER, PARAMETER :: nrands=1

REAL :: WGHT0, PRCBL, rndu(nrands)
INTEGER :: KK, L

do kk=icmin,k+1
   PRJ(kk) = PKE(I,kk)
enddo

poi=0. ! These initialized here in order not to confuse Valgrind debugger
qoi=0. ! Do not believe it actually makes any difference.
uoi=0.
voi=0.     

PRS(ICMIN:K0+1) = PLE(I,ICMIN:K0+1)
POI(ICMIN:K)   = THO(I,ICMIN:K)
QOI(ICMIN:K)   = QHO(I,ICMIN:K)
UOI(ICMIN:K)   = UHO(I,ICMIN:K)
VOI(ICMIN:K)   = VHO(I,ICMIN:K)

QST(ICMIN:K) = QSS(I,ICMIN:K)
DQQ(ICMIN:K) = DQS(I,ICMIN:K)

MASSF(:) = WGT0(:)

!!! RESET PRESSURE at bottom edge of CBL 
PRCBL = PRS(K)
do l= K,K0
   PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
end do

PRS(K+1) = PRCBL
PRJ(K+1) = (PRS(K+1)/1000.)**(MAPL_RGAS/MAPL_CP)

DO L=K,ICMIN,-1
   POL(L)  = 0.5*(PRS(L)+PRS(L+1))
   PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) / (ONEPKAP*(PRS(L+1)-PRS(L)))
   PKI(L)  = 1.0 / PRH(L)
   DPT(L)  = PRH(L  ) - PRJ(L)
   DPB(L)  = PRJ(L+1) - PRH(L)
   PRI(L)  = .01 / (PRS(L+1)-PRS(L))
ENDDO


!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
if ( K <= K0) then
   POI(K) = 0.
   QOI(K) = 0.
   UOI(K) = 0.
   VOI(K) = 0.

   !! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
   WGHT = 0.
   DO L=K,K0
      WGHT(L)   = MASSF(L) * ( PLE(I,L+1)-PLE(I,L) )/( PRS(K+1)-PRS(K)  )
   END DO      

   DO L = K,K0
      POI(K) = POI(K) + WGHT(L)*THO(I,L)
      QOI(K) = QOI(K) + WGHT(L)*QHO(I,L)
      UOI(K) = UOI(K) + WGHT(L)*UHO(I,L)
      VOI(K) = VOI(K) + WGHT(L)*VHO(I,L)
   ENDDO

   call MRQSAT_scalar(POI(K)*PRH(K),POL(K),QST(K),DQQ(K))

endif

!SET MAX DIAM OF CONVECTION, GOES INTO LAMBDA CALCULATION
rndu(:) = max( seedras(1)/1000000., 1e-6 )
MXDIAM(I) = maxdallowed*( rndu(1)**(-1./2.) )


DO L=K,ICMIN,-1
   BET(L)  = DQQ(L)*PKI(L)
   GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L))
   IF (L<K) THEN
      GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
      GM1(L+1) = 0.5*LBCP*(DQQ(L  )/(ALHL*(1.0+LBCP*DQQ(L  ))) + DQQ(L+1)/(ALHL*(1.0+LBCP*DQQ(L+1))) )
   ENDIF
ENDDO


!USED IN RESTRAP (STRAP = 1) CALCULATION
POI_SV = POI
QOI_SV = QOI
UOI_SV = UOI
VOI_SV = VOI

END SUBROUTINE STRAP0
!======================================================================================================================


!=======================================================================================================
subroutine STRAP1(IDIM, I, K, K0, ICMIN, DDT, DAYLEN, POI, QOI, UOI, VOI, POI_SV, QOI_SV, UOI_SV, VOI_SV, &
                  WGT1, PLE, PRS, THO, QHO, UHO, VHO )
!=======================================================================================================

IMPLICIT NONE

!!INPUTS!!
INTEGER :: IDIM, I, K, K0, ICMIN
REAL :: DDT, DAYLEN

REAL, DIMENSION(K0) :: POI, QOI, UOI, VOI, POI_SV, QOI_SV, UOI_SV, VOI_SV, WGT1
REAL, DIMENSION(IDIM, K0+1) :: PLE
REAL, DIMENSION(K0+1) :: PRS

!!OUTPUTS!!
REAL, DIMENSION(IDIM,K0) :: THO, QHO, UHO, VHO 

!!LOCALS!!
REAL, DIMENSION(K0) :: WGHT, WGHT0
INTEGER :: L

THO(I,ICMIN:K-1) = POI(ICMIN:K-1)
QHO(I,ICMIN:K-1) = QOI(ICMIN:K-1)
UHO(I,ICMIN:K-1) = UOI(ICMIN:K-1)
VHO(I,ICMIN:K-1) = VOI(ICMIN:K-1)
   
WGHT   = WGT1
wght0 = 0.
DO L=K,K0 
   wght0 = wght0 + WGHT(L)* ( PLE(I,L+1) - PLE(I,L) )
END DO
wght0 = ( PRS(K+1)   - PRS(K)  )/wght0
WGHT  = wght0 * WGHT

DO L=K,K0 
   THO(I,L) =  THO(I,L) + WGHT(L)*(POI(K) - POI_SV(K))
   QHO(I,L) =  QHO(I,L) + WGHT(L)*(QOI(K) - QOI_SV(K))
   UHO(I,L) =  UHO(I,L) + WGHT(L)*(UOI(K) - UOI_SV(K))
   VHO(I,L) =  VHO(I,L) + WGHT(L)*(VOI(K) - VOI_SV(K))
END DO

end subroutine STRAP1
!=======================================================================================================

!============================================================================================================================================
SUBROUTINE CLOUDE(N_DTL,ICMIN,IDIM,I,K,K0,DT,DDT,DAYLEN,RHMN,RHMX,AUTORAMPB,CO_ZDEP,CP,&
RHMAX,GRAV,ALHL,ALHI,CPBG,LBCP,RASAL1,RASAL2,ACRITFAC,SDQVT1,SDQV2,SDQV3,CO_AUTO,CLI_CRIT,&
PBLFRAC,GRAVI,CPI,FRICFAC,FRICLAMBDA,TPERT,QPERT,MXDIAM,PRJ,PRH,PRS,SIGE,&
DPB,DPT,DQQ,POL,GM1,PKI,BET,GHT,POI,QOI,UOI,VOI,QST,GAM,PRI,RMF)
!============================================================================================================================================

IMPLICIT NONE

!!INPUTS!!
INTEGER :: N_DTL, IDIM, I, K, K0, ICMIN
REAL :: DT, DDT, DAYLEN, RHMN, RHMX, AUTORAMPB, CO_ZDEP, CP, RHMAX, GRAV, ALHL, ALHI, CPBG, LBCP
REAL :: RASAL1, RASAL2, ACRITFAC, SDQVT1, SDQV2, SDQV3
REAL :: CO_AUTO, CLI_CRIT, PBLFRAC, GRAVI, CPI, FRICFAC, FRICLAMBDA

REAL, DIMENSION(IDIM) ::  TPERT, QPERT, MXDIAM

REAL, DIMENSION(K0+1) :: PRJ, PRH, PRS, SIGE
REAL, DIMENSION(K0) :: DPB, DPT, DQQ, POL, GM1, PKI, BET, GHT

!!INOUTS!!
REAL, DIMENSION(K0) :: POI, QOI, UOI, VOI

REAL, DIMENSION(K0) :: QST, GAM, PRI
REAL, DIMENSION(K0) :: RMF

!!LOCALS!!

!NOT COLUMN DEP
REAL :: LAMBDA_MIN, LAMBDA_MAX, F3, JUMP1
REAL, PARAMETER :: TE0 = 273.0, TE2 = 200.0


!SAVED - ALL STEPS
INTEGER :: L, ICL_C
INTEGER, DIMENSION(N_DTL) :: RC_ALL

REAL, DIMENSION(N_DTL) :: TRG_ALL, ALM_ALL, F4_ALL, RASAL_ALL
REAL, DIMENSION(N_DTL) :: TOKI_ALL, WFN_ALL, ACR_ALL, TEM_ALL

REAL, DIMENSION(K0,N_DTL) :: POIc, QOIc, HCC_ALL
REAL, DIMENSION(K0,N_DTL) :: ETA_ALL, EHT_ALL, HCLD_ALL, HCLD1_ALL, CVW_ALL
REAL, DIMENSION(K0,N_DTL) :: BK2_ALL, QHT_ALL, TEMV1_ALL, TEMV2_ALL, TEMV3_ALL, TEMV4_ALL

!NONLOCAL BUT ARRAY SAVED
REAL, DIMENSION(K0,N_DTL) :: QOL_ALL, SSL_ALL, HOL_ALL, HST_ALL, ZOL_ALL
REAL, DIMENSION(K0+1,N_DTL) :: ZET_ALL, SHT_ALL


!SAVED - LOOP4 ONLY
INTEGER, DIMENSION(N_DTL) :: L4_MUL

REAL, DIMENSION(N_DTL) :: WLQ1_AL4

REAL,  DIMENSION(K0,N_DTL) :: F2_AL4, TX2a_AL4, CLL0_AL4, TEMV_AL4, UHT_AL4, VHT_AL4, C00_X_AL4
REAL, DIMENSION(K0,N_DTL) :: TE_A_AL4, QCC_AL4, GMS_AL4, GMH_AL4, WLQ_AL4, WLQTEM_AL4, TX2_AL4, CLI_AL4
REAL, DIMENSION(K0,N_DTL) :: CLI_CRIT_X_AL4, RATE_AL4, RATE_AL4H, CVW_X_AL4, DT_LYR_AL4, CLOSS_AL4, AKM_AL4

!SAVED - LOOP5 ONLY
INTEGER, DIMENSION(N_DTL) :: L5_MUL

REAL, DIMENSION(N_DTL) :: TEM_AL5, WFN1_AL5, WFN2_AL5, WFN2_AL5H

REAL, DIMENSION(K0,N_DTL) :: GMS_AL5, GMH_AL5, TEMV_AL5
REAL, DIMENSION(K0,N_DTL+1) :: UCU_AL5, VCU_AL5

!SAVED - MODEL VARIABLES
REAL, DIMENSION(K0,N_DTL+1) :: POI_ALL, QOI_ALL, QST_ALL, UOI_ALL, VOI_ALL


!ZERO OUT LOCAL ARRAYS
RC_ALL = 0
L4_MUL = 1
!L5_MUL = 1
TRG_ALL = 0.0
ALM_ALL = 0.0
F4_ALL = 0.0
RASAL_ALL = 0.0
TOKI_ALL = 0.0
WFN_ALL = 0.0 
ACR_ALL = 0.0
POIc = 0.0
QOIc = 0.0
HCC_ALL = 0.0
ETA_ALL = 0.0
EHT_ALL = 0.0
HCLD_ALL = 0.0
HCLD1_ALL = 0.0
CVW_ALL = 0.0
BK2_ALL = 0.0
QHT_ALL = 0.0
TEM_ALL = 0.0
TEMV1_ALL = 0.0
TEMV2_ALL = 0.0
TEMV3_ALL = 0.0
TEMV4_ALL = 0.0
QOL_ALL = 0.0
SSL_ALL = 0.0
HOL_ALL = 0.0
HST_ALL = 0.0
ZOL_ALL = 0.0
ZET_ALL = 0.0
SHT_ALL = 0.0
CLOSS_AL4 = 0.0
CVW_X_AL4 = 0.0
RATE_AL4 = 0.0
RATE_AL4H = 0.0
DT_LYR_AL4 = 0.0
CLI_CRIT_X_AL4 = 0.0
C00_X_AL4 = 0.0
CLI_AL4 = 0.0
TEMV_AL4 = 0.0
QCC_AL4 = 0.0
TX2_AL4 = 0.0
TE_A_AL4 = 0.0
UHT_AL4 = 0.0
VHT_AL4 = 0.0
WLQ_AL4 = 0.0
WLQTEM_AL4 = 0.0
WLQ1_AL4 = 0.0
AKM_AL4 = 0.0
F2_AL4 = 0.0
TX2a_AL4 = 0.0
CLL0_AL4 = 0.0
GMS_AL4 = 0.0
GMH_AL4 = 0.0
TEM_AL5 = 0.0
TEMV_AL5 = 0.0
WFN1_AL5 = 0.0
WFN2_AL5 = 0.0
WFN2_AL5H = 0.0
UCU_AL5 = 0.0
VCU_AL5 = 0.0
GMS_AL5 = 0.0
GMH_AL5 = 0.0
POI_ALL = 0.0
QOI_ALL = 0.0
QST_ALL = 0.0
UOI_ALL = 0.0
VOI_ALL = 0.0


!SAVING OF UPDATED BUT NONLOCAL
POI_ALL(:,1) = POI
QOI_ALL(:,1) = QOI
UOI_ALL(:,1) = UOI
VOI_ALL(:,1) = VOI
QST_ALL(:,1) = QST


JUMP1 =  (SDQV2-1.0) / ( ( TE0-SDQVT1 )**0.333 )

!CALCULATE LAMBDA
LAMBDA_MIN = .2/MXDIAM(I)
LAMBDA_MAX = .2/  200. 

F3 = 1.0

DO ICL_C = 1,N_DTL

   RC_ALL(ICL_C) = 10
 
   ALM_ALL(ICL_C)   = 0.
   TRG_ALL(ICL_C)   = AMIN1(1.,(QOI_ALL(K,ICL_C)/QST_ALL(K,ICL_C)-RHMN)/(RHMX-RHMN))
   
   F4_ALL(ICL_C)  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(K-ICL_C))/0.2 )  )
   
   !SAVE TH AND Q AND KICK
   POIc(:,ICL_C) = POI_ALL(:,ICL_C)
   QOIc(:,ICL_C) = QOI_ALL(:,ICL_C)

   POIc(K,ICL_C) =  POIc(K,ICL_C) + TPERT(I)
   QOIc(K,ICL_C) =  QOIc(K,ICL_C) + QPERT(I)
   
   ZET_ALL(K+1,ICL_C) = 0.
   SHT_ALL(K+1,ICL_C) = CP*POIc(K,ICL_C)*PRJ(K+1)
   DO L=K,K-ICL_C,-1
      QOL_ALL(L,ICL_C)  = AMIN1(QST_ALL(L,ICL_C)*RHMAX,QOIc(L,ICL_C))
      QOL_ALL(L,ICL_C)  = AMAX1( 0.000, QOL_ALL(L,ICL_C) )     ! GUARDIAN/NEG.PREC.
      SSL_ALL(L,ICL_C)  = CP*PRJ(L+1)*POIc(L,ICL_C) + GRAV*ZET_ALL(L+1,ICL_C)
      HOL_ALL(L,ICL_C)  = SSL_ALL(L,ICL_C) + QOL_ALL(L,ICL_C)*ALHL
      HST_ALL(L,ICL_C)  = SSL_ALL(L,ICL_C) + QST_ALL(L,ICL_C)*ALHL
      TEMV1_ALL(L,ICL_C)     = POIc(L,ICL_C)*(PRJ(L+1)-PRJ(L))*CPBG
      ZET_ALL(L,ICL_C)  = ZET_ALL(L+1,ICL_C) + TEMV1_ALL(L,ICL_C)
      ZOL_ALL(L,ICL_C)  = ZET_ALL(L+1,ICL_C) + (PRJ(L+1)-PRH(L))*POIc(L,ICL_C)*CPBG
   ENDDO
   
   DO L=K-ICL_C+1,K
      TEMV2_ALL(L,ICL_C)  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
      SHT_ALL(L,ICL_C)  = SSL_ALL(L-1,ICL_C) + TEMV2_ALL(L,ICL_C)*(SSL_ALL(L,ICL_C)-SSL_ALL(L-1,ICL_C)) 
      QHT_ALL(L,ICL_C)  = .5*(QOL_ALL(L,ICL_C)+QOL_ALL(L-1,ICL_C))
   ENDDO 
      
   IF (HOL_ALL(K,ICL_C) <= HST_ALL(K-ICL_C,ICL_C)) THEN   !CANNOT REACH IC LEVEL
      RC_ALL(ICL_C) = 1
      !====EXIT====>
   endIF
      
   !LAMBDA CALCULATION: MS-A18
   TEM_ALL(ICL_C) = (HST_ALL(K-ICL_C,ICL_C)-HOL_ALL(K-ICL_C,ICL_C))*(ZOL_ALL(K-ICL_C,ICL_C)-ZET_ALL(K-ICL_C+1,ICL_C)) 
   DO L=K-ICL_C+1,K-1
      TEM_ALL(ICL_C) = TEM_ALL(ICL_C) + (HST_ALL(K-ICL_C,ICL_C)-HOL_ALL(L,ICL_C))*(ZET_ALL(L ,ICL_C)-ZET_ALL(L +1,ICL_C))
   ENDDO
   
   IF (TEM_ALL(ICL_C) <= 0.0) THEN !NO VALID LAMBDA     
      RC_ALL(ICL_C) = 2
      !====EXIT====>
   endIF
         
   ALM_ALL(ICL_C) = (HOL_ALL(K,ICL_C)-HST_ALL(K-ICL_C,ICL_C)) / TEM_ALL(ICL_C)
         
   IF (ALM_ALL(ICL_C) > LAMBDA_MAX) THEN
      RC_ALL(ICL_C) = 3
      !====EXIT====>
   endif   
   
   TOKI_ALL(ICL_C) = 1.0
   
   IF (ALM_ALL(ICL_C) < LAMBDA_MIN) THEN
      TOKI_ALL(ICL_C) = ( ALM_ALL(ICL_C)/LAMBDA_MIN )**2
   ENDIF
   
   !ETA CALCULATION
   DO L=K-ICL_C+1,K
      ETA_ALL(L,ICL_C) = 1.0 + ALM_ALL(ICL_C) * (ZET_ALL(L,ICL_C )-ZET_ALL(K,ICL_C))
   ENDDO
   
   ETA_ALL(K-ICL_C,ICL_C) = 1.0 + ALM_ALL(ICL_C) * (ZOL_ALL(K-ICL_C,ICL_C)-ZET_ALL(K,ICL_C))
         
   !WORKFUNCTION CALCULATION:  MS-A22
   WFN_ALL(ICL_C)     = 0.0
   HCC_ALL(K,ICL_C)  = HOL_ALL(K,ICL_C)
   DO L=K-1,K-ICL_C+1,-1
      HCC_ALL(L,ICL_C) = HCC_ALL(L+1,ICL_C) + ( ETA_ALL(L,ICL_C) - ETA_ALL(L+1,ICL_C) )*HOL_ALL(L,ICL_C)
      TEMV3_ALL(L,ICL_C)    = HCC_ALL(L+1,ICL_C)*DPB(L) + HCC_ALL(L,ICL_C)*DPT(L)
      EHT_ALL(L,ICL_C) = ETA_ALL(L+1,ICL_C)*DPB(L) + ETA_ALL(L,ICL_C)*DPT(L)
      WFN_ALL(ICL_C)    = WFN_ALL(ICL_C) + (TEMV3_ALL(L,ICL_C) - EHT_ALL(L,ICL_C)*HST_ALL(L,ICL_C))*GAM(L)
   ENDDO
   HCC_ALL(K-ICL_C,ICL_C) = HST_ALL(K-ICL_C,ICL_C)*ETA_ALL(K-ICL_C,ICL_C)
   WFN_ALL(ICL_C)     = WFN_ALL(ICL_C) + (HCC_ALL(K-ICL_C+1,ICL_C)&
                    -HST_ALL(K-ICL_C,ICL_C)*ETA_ALL(K-ICL_C+1,ICL_C))*GAM(K-ICL_C)*DPB(K-ICL_C)
   
   !VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
   HCLD_ALL(K,ICL_C) = HOL_ALL(K,ICL_C)
   BK2_ALL(K,ICL_C) = 0.0
   DO L=K-1,K-ICL_C,-1
      HCLD1_ALL(L,ICL_C) = ( ETA_ALL(L+1,ICL_C)*HCLD_ALL(L+1,ICL_C)  +  (ETA_ALL(L,ICL_C) &
                           - ETA_ALL(L+1,ICL_C))*HOL_ALL(L,ICL_C) ) / ETA_ALL(L,ICL_C)
      TEMV4_ALL(L,ICL_C)     = (HCLD1_ALL(L,ICL_C)-HST_ALL(L,ICL_C) )*(ZET_ALL(L,ICL_C)-ZET_ALL(L+1,ICL_C))/ (1.0+LBCP*DQQ(L))
      BK2_ALL(L,ICL_C)  = BK2_ALL(L+1,ICL_C) + GRAV * TEMV4_ALL(L,ICL_C) / ( CP*PRJ(L+1)*POI_ALL(L,ICL_C) )
      CVW_ALL(L,ICL_C) = SQRT(  2.0* MAX( BK2_ALL(L,ICL_C) , 1.0 )  )
   ENDDO
   

   !ALPHA CALCULATION 
   IF ( ZET_ALL(K-ICL_C,ICL_C) <  2000. ) THEN
      RASAL_ALL(ICL_C) = DT / RASAL1
   ELSE
      RASAL_ALL(ICL_C) = DT / (RASAL1 + (RASAL2-RASAL1)*(ZET_ALL(K-ICL_C,ICL_C) - 2000.)/8000.)
   ENDIF
   
   CVW_ALL(K-ICL_C:K,ICL_C) = MAX( CVW_ALL(K-ICL_C:K,ICL_C), 1.00 )
   
   !TEST FOR CRITICAL WORK FUNCTION
   CALL ACRITN(POL(K-ICL_C),PRS(K),ACR_ALL(ICL_C),ACRITFAC)

   IF (WFN_ALL(ICL_C) <= ACR_ALL(ICL_C)) THEN   ! SUB-CRITICAL WORK FUNCTION
      RC_ALL(ICL_C) = 4
      !====EXIT====>
   ENDIF
      
   !IF (RC_ALL(ICL_C) == 10) THEN

   !CLOUD TOP WATER AND MOMENTUM
   WLQ_AL4(K,ICL_C)   = QOL_ALL(K,ICL_C)
   UHT_AL4(K,ICL_C) = UOI_ALL(K,ICL_C)
   VHT_AL4(K,ICL_C) = VOI_ALL(K,ICL_C)
   CLL0_AL4(K,ICL_C) = 0.
      
   DO L=K-1,K-ICL_C,-1
      TEMV_AL4(L,ICL_C)  = ETA_ALL(L,ICL_C) - ETA_ALL(L+1,ICL_C)
      WLQ_AL4(L,ICL_C)   = WLQ_AL4(L+1,ICL_C) + TEMV_AL4(L,ICL_C) * QOL_ALL(L,ICL_C)
      UHT_AL4(L,ICL_C)   = UHT_AL4(L+1,ICL_C) + TEMV_AL4(L,ICL_C) * UOI_ALL(L,ICL_C)
      VHT_AL4(L,ICL_C)   = VHT_AL4(L+1,ICL_C) + TEMV_AL4(L,ICL_C) * VOI_ALL(L,ICL_C)
      
      WLQTEM_AL4(L,ICL_C) = WLQ_AL4(L,ICL_C)

      ! How much condensate (CLI_AL4(ICL_C)) is present here? 
      IF (L>K-ICL_C) THEN
         TX2_AL4(L,ICL_C) = 0.5*(QST_ALL(L,ICL_C)+QST_ALL(L-1,ICL_C))*ETA_ALL(L,ICL_C)
         QCC_AL4(L,ICL_C) = TX2_AL4(L,ICL_C) + GM1(L)*(HCC_ALL(L,ICL_C)&
                 -0.5*(HST_ALL(L,ICL_C)+HST_ALL(L-1,ICL_C))*ETA_ALL(L,ICL_C))
         CLL0_AL4(L,ICL_C) = (WLQTEM_AL4(L,ICL_C)-QCC_AL4(L,ICL_C))
      ELSE
         CLL0_AL4(L,ICL_C)   = (WLQTEM_AL4(L,ICL_C)-QST_ALL(K-ICL_C,ICL_C)*ETA_ALL(K-ICL_C,ICL_C))
      ENDIF
      CLL0_AL4(L,ICL_C)    = MAX(CLL0_AL4(L,ICL_C), 0.00)
      
      CLI_AL4(L,ICL_C)  = CLL0_AL4(L,ICL_C) / ETA_ALL(L,ICL_C)  !Condensate (kg/kg)
      TE_A_AL4(L,ICL_C) = POI_ALL(L,ICL_C)*PRH(L)     !Temperature (K)
      
      !!!!!!!!!!!!!!!!!! SUNDQ3 CALL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      F2_AL4(L,ICL_C)   = 1.0
      IF ( ( TE_A_AL4(L,ICL_C) .GE. SDQVT1 ) .AND. ( TE_A_AL4(L,ICL_C) .LT. TE0 ) )THEN 
         F2_AL4(L,ICL_C)   = 1.0 + JUMP1 * (( TE0 - TE_A_AL4(L,ICL_C) )**0.3333)
      ENDIF
      IF ( TE_A_AL4(L,ICL_C) .LT. SDQVT1 ) THEN 
         F2_AL4(L,ICL_C)   = SDQV2 !+ (SDQV3-SDQV2)*(SDQVT1-TE_A_AL4(L,ICL_C))/(SDQVT1-TE2)
      ENDIF
      !!!!!!!!!!!!!!!!!! END SUNDQ3 CALL !!!!!!!!!!!!!!!!!!!!!!!!
   
      C00_X_AL4(L,ICL_C)  =  CO_AUTO * F2_AL4(L,ICL_C) * F3  * F4_ALL(ICL_C)
      CLI_CRIT_X_AL4(L,ICL_C) =  CLI_CRIT / ( F2_AL4(L,ICL_C)*F3 )

      RATE_AL4H(L,ICL_C) = EXP( -(CLI_AL4(L,ICL_C))**2 / CLI_CRIT_X_AL4(L,ICL_C)**2 )
      RATE_AL4(L,ICL_C) = C00_X_AL4(L,ICL_C) * ( 1.0 - RATE_AL4H(L,ICL_C) )
      
      CVW_X_AL4(L,ICL_C) = MAX( CVW_ALL(L,ICL_C), 1.00 ) !Floor convective velocity, don't trust at low values.
      
      DT_LYR_AL4(L,ICL_C) = ( ZET_ALL(L,ICL_C)-ZET_ALL(L+1,ICL_C) )/CVW_X_AL4(L,ICL_C)
      
      CLOSS_AL4(L,ICL_C) = CLL0_AL4(L,ICL_C) * RATE_AL4(L,ICL_C) * DT_LYR_AL4(L,ICL_C)
      CLOSS_AL4(L,ICL_C) = MIN( CLOSS_AL4(L,ICL_C) , CLL0_AL4(L,ICL_C) )
      
      CLL0_AL4(L,ICL_C) = CLL0_AL4(L,ICL_C) - CLOSS_AL4(L,ICL_C)
      
      IF (CLOSS_AL4(L,ICL_C) > 0.) then
         WLQ_AL4(L,ICL_C) = WLQ_AL4(L,ICL_C) - CLOSS_AL4(L,ICL_C)
      ENDIF
      
   ENDDO

   !endif
   
   IF (RC_ALL(ICL_C) == 10) THEN
      L4_MUL(ICL_C) = 1
   ELSE 
      L4_MUL(ICL_C) = 0
   ENDIF

   WLQ1_AL4(ICL_C) = L4_MUL(ICL_C)*( WLQ_AL4(K-ICL_C,ICL_C) - QST_ALL(K-ICL_C,ICL_C)*ETA_ALL(K-ICL_C,ICL_C) )

   !CALCULATE GAMMAS AND KERNEL
   GMS_AL4(K,ICL_C) = L4_MUL(ICL_C)*( (SHT_ALL(K,ICL_C)-SSL_ALL(K,ICL_C))*PRI(K)   )
   GMH_AL4(K,ICL_C) = L4_MUL(ICL_C)*( GMS_AL4(K,ICL_C) + (QHT_ALL(K,ICL_C)-QOL_ALL(K,ICL_C))*PRI(K)*ALHL   )
   TX2a_AL4(K,ICL_C)  = L4_MUL(ICL_C)*( GMH_AL4(K,ICL_C) )
   AKM_AL4(K,ICL_C) = L4_MUL(ICL_C)*( GMH_AL4(K,ICL_C)*GAM(K-1)*DPB(K-1)      )
        
   GMS_AL4(K-ICL_C+1:K-1,ICL_C) = L4_MUL(ICL_C)*( ( ETA_ALL(K-ICL_C+1:K-1,ICL_C)*(SHT_ALL(K-ICL_C+1:K-1,ICL_C)&
                                  -SSL_ALL(K-ICL_C+1:K-1,ICL_C)) + ETA_ALL(K-ICL_C+2:K,ICL_C)*&
                                  (SSL_ALL(K-ICL_C+1:K-1,ICL_C)-SHT_ALL(K-ICL_C+2:K,ICL_C)) )*PRI(K-ICL_C+1:K-1)  ) 
   GMH_AL4(K-ICL_C+1:K-1,ICL_C) = L4_MUL(ICL_C)*( GMS_AL4(K-ICL_C+1:K-1,ICL_C) &
                                  + ( ETA_ALL(K-ICL_C+1:K-1,ICL_C)*(QHT_ALL(K-ICL_C+1:K-1,ICL_C) &
                                  - QOL_ALL(K-ICL_C+1:K-1,ICL_C)) + ETA_ALL(K-ICL_C+2:K,ICL_C)*&
                                  (QOL_ALL(K-ICL_C+1:K-1,ICL_C) - QHT_ALL(K-ICL_C+2:K,ICL_C)) )*ALHL*PRI(K-ICL_C+1:K-1)   )
   DO L=K-1,K-ICL_C+1,-1
      TX2a_AL4(L,ICL_C) = L4_MUL(ICL_C)*( TX2a_AL4(L+1,ICL_C) + (ETA_ALL(L,ICL_C) - ETA_ALL(L+1,ICL_C)) * GMH_AL4(L,ICL_C)  )
      AKM_AL4(L,ICL_C) = L4_MUL(ICL_C)*( AKM_AL4(L+1,ICL_C) - &
                            GMS_AL4(L,ICL_C)*EHT_ALL(L,ICL_C)*PKI(L) + TX2a_AL4(L,ICL_C)*GHT(L)  )
   ENDDO
      
   GMS_AL4(K-ICL_C,ICL_C) = L4_MUL(ICL_C)*( ETA_ALL(K-ICL_C+1,ICL_C)*(SSL_ALL(K-ICL_C,ICL_C)&
                                 -SHT_ALL(K-ICL_C+1,ICL_C))*PRI(K-ICL_C)  )
   AKM_AL4(K-ICL_C,ICL_C) = L4_MUL(ICL_C)*( AKM_AL4(K-ICL_C+1,ICL_C)-&
                               GMS_AL4(K-ICL_C,ICL_C)*ETA_ALL(K-ICL_C+1,ICL_C)*DPB(K-ICL_C)*PKI(K-ICL_C)  )
      
   GMH_AL4(K-ICL_C,ICL_C) = L4_MUL(ICL_C)*( GMS_AL4(K-ICL_C,ICL_C) &
                            + ( ETA_ALL(K-ICL_C+1,ICL_C)*(QOL_ALL(K-ICL_C,ICL_C)-QHT_ALL(K-ICL_C+1,ICL_C))*ALHL &
                           + ETA_ALL(K-ICL_C,ICL_C)*(HST_ALL(K-ICL_C,ICL_C)-HOL_ALL(K-ICL_C,ICL_C)) )*PRI(K-ICL_C)  )
      
   

   !CLOUD BASE MASS FLUX
   IF  ((AKM_AL4(K-ICL_C,ICL_C) >= 0.0 .OR. WLQ1_AL4(ICL_C) < 0.0)) THEN
      RC_ALL(ICL_C) = 5
      !====EXIT====>
   endIF

   IF (RC_ALL(ICL_C) == 10) THEN
      L5_MUL(ICL_C) = 1
   ELSE 
      L5_MUL(ICL_C) = 0
   ENDIF
       
   IF (RC_ALL(ICL_C) == 10) THEN

   WFN1_AL5(ICL_C) = L5_MUL(ICL_C)*( MIN((RASAL_ALL(ICL_C)*TRG_ALL(ICL_C)*TOKI_ALL(ICL_C))*-(WFN_ALL(ICL_C)-ACR_ALL(ICL_C))&
                          /AKM_AL4(K-ICL_C,ICL_C),(PRS(K+1)-PRS(K))*(100.*PBLFRAC) )  )
  
   !CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
   TEM_AL5(ICL_C) = L5_MUL(ICL_C)*( WFN1_AL5(ICL_C)*GRAVI )
   RMF (K-ICL_C) = L5_MUL(ICL_C)*( RMF (K-ICL_C) + TEM_AL5(ICL_C) )
          
   !THETA AND Q CHANGE DUE TO CLOUD TYPE IC
   GMH_AL5(K-ICL_C:K,ICL_C) = L5_MUL(ICL_C)*( GMH_AL4(K-ICL_C:K,ICL_C) * WFN1_AL5(ICL_C) )
   GMS_AL5(K-ICL_C:K,ICL_C) = L5_MUL(ICL_C)*( GMS_AL4(K-ICL_C:K,ICL_C) * WFN1_AL5(ICL_C) )
 
   !CUMULUS FRICTION
   WFN2_AL5H(ICL_C) = EXP( -ALM_ALL(ICL_C) / FRICLAMBDA )
   WFN2_AL5(ICL_C) = 0.5*WFN1_AL5(ICL_C)*FRICFAC*WFN2_AL5H(ICL_C)

   TEMV_AL5(K-ICL_C:K,ICL_C) = WFN2_AL5(ICL_C)*PRI(K-ICL_C:K)
      
   ENDIF   

   UCU_AL5(K,ICL_C)  = UCU_AL5(K,ICL_C) + L5_MUL(ICL_C)*( TEMV_AL5(K,ICL_C) * (UOI_ALL(K-1,ICL_C) - UOI_ALL(K,ICL_C)) )
   VCU_AL5(K,ICL_C)  = VCU_AL5(K,ICL_C) + L5_MUL(ICL_C)*( TEMV_AL5(K,ICL_C) * (VOI_ALL(K-1,ICL_C) - VOI_ALL(K,ICL_C)) )
      
   DO L=K-1,K-ICL_C+1,-1
      UCU_AL5(L,ICL_C) = UCU_AL5(L,ICL_C) + L5_MUL(ICL_C)*( TEMV_AL5(L,ICL_C) * ( (UOI_ALL(L-1,ICL_C) &
                         - UOI_ALL(L,ICL_C)) * ETA_ALL(L,ICL_C) &
                         + (UOI_ALL(L,ICL_C) - UOI_ALL(L+1,ICL_C)) * ETA_ALL(L+1,ICL_C) ) )
      VCU_AL5(L,ICL_C) = VCU_AL5(L,ICL_C) + L5_MUL(ICL_C)*( TEMV_AL5(L,ICL_C) * ( (VOI_ALL(L-1,ICL_C) &
                         - VOI_ALL(L,ICL_C)) * ETA_ALL(L,ICL_C) &
                         + (VOI_ALL(L,ICL_C) - VOI_ALL(L+1,ICL_C)) * ETA_ALL(L+1,ICL_C) ) )
   ENDDO

   UCU_AL5(K-ICL_C,ICL_C) = UCU_AL5(K-ICL_C,ICL_C) + L5_MUL(ICL_C)*( (2.*(UHT_AL4(K-ICL_C,ICL_C) &
                          - UOI_ALL(K-ICL_C,ICL_C)*(ETA_ALL(K-ICL_C,ICL_C)-ETA_ALL(K-ICL_C+1,ICL_C))) &
                          - (UOI_ALL(K-ICL_C,ICL_C)+UOI_ALL(K-ICL_C+1,ICL_C))&
                          *ETA_ALL(K-ICL_C+1,ICL_C))*TEMV_AL5(K-ICL_C,ICL_C) )
   VCU_AL5(K-ICL_C,ICL_C) = VCU_AL5(K-ICL_C,ICL_C) + L5_MUL(ICL_C)*( (2.*(VHT_AL4(K-ICL_C,ICL_C) &
                          - VOI_ALL(K-ICL_C,ICL_C)*(ETA_ALL(K-ICL_C,ICL_C)-ETA_ALL(K-ICL_C+1,ICL_C))) &
                          - (VOI_ALL(K-ICL_C,ICL_C)+VOI_ALL(K-ICL_C+1,ICL_C))&
                          *ETA_ALL(K-ICL_C+1,ICL_C))*TEMV_AL5(K-ICL_C,ICL_C) )

   
   QOI_ALL(:,ICL_C+1) = QOI_ALL(:,ICL_C)
   POI_ALL(:,ICL_C+1) = POI_ALL(:,ICL_C)
   QST_ALL(:,ICL_C+1) = QST_ALL(:,ICL_C)

   UOI_ALL(:,ICL_C+1) = UOI_ALL(:,ICL_C)
   VOI_ALL(:,ICL_C+1) = VOI_ALL(:,ICL_C)

   DO L=K-ICL_C,K
      QOI_ALL(L,ICL_C+1) = QOI_ALL(L,ICL_C) + (GMH_AL5(L,ICL_C) - GMS_AL5(L,ICL_C)) * ALHI
      POI_ALL(L,ICL_C+1) = POI_ALL(L,ICL_C) +  GMS_AL5(L,ICL_C) * PKI(L)*CPI
      QST_ALL(L,ICL_C+1) = QST_ALL(L,ICL_C) +  GMS_AL5(L,ICL_C) * BET(L)*CPI

      UOI_ALL(L,ICL_C+1) = UOI_ALL(L,ICL_C) + UCU_AL5(L,ICL_C)
      VOI_ALL(L,ICL_C+1) = VOI_ALL(L,ICL_C) + VCU_AL5(L,ICL_C)
   ENDDO
      
ENDDO

QOI = QOI_ALL(:,N_DTL+1)
POI = POI_ALL(:,N_DTL+1)
QST = QST_ALL(:,N_DTL+1)
UOI = UOI_ALL(:,N_DTL+1)
VOI = VOI_ALL(:,N_DTL+1)

END SUBROUTINE CLOUDE
!============================================================================================================================================

!=========================================
SUBROUTINE ACRITN( PL, PLB, ACR, ACRITFAC)
!=========================================
      
      IMPLICIT NONE

      REAL :: PL, PLB, ACRITFAC
      REAL :: ACR
      INTEGER IWK
      
      !!REAL, PARAMETER :: FACM=0.5

      REAL, PARAMETER :: &
         PH(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, &
                 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)

      REAL, PARAMETER :: & 
         A(15)=(/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677, &
                 0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
                 0.0553, 0.0445, 0.0633     /)   !!*FACM

      IWK = PL * 0.02 - 0.999999999

      IF (IWK .GT. 1 .AND. IWK .LE. 15) THEN
       ACR = A(IWK-1) + (PL-PH(IWK-1))*.02*(A(IWK)-A(IWK-1))
      ELSEIF(IWK > 15) THEN
       ACR = A(15)
      ELSE
       ACR = A(1)
      ENDIF

      ACR = ACRITFAC  * ACR * (PLB - PL)


   END SUBROUTINE ACRITN
!========================================


!======================================
SUBROUTINE MRQSAT_scalar(TT,P,Qs,DQsDT)
!======================================

      implicit none

      !Inputs
      real :: TT, P, Qs, DQsDT
      logical :: LDQDT

      !Constants
      real, parameter :: svpt0=273.15, svp1=.6112, svp2=17.67, svp3=29.65, airmw = 28.97, h2omw = 18.01, one = 1.0, ep2 = 0.622
      real, parameter :: esfac = h2omw/airmw, erfac = (one-esfac)/esfac

      !Internals
      real :: dqdta, qsatvp, d

      
      qsatvp = svp1*exp(svp2*(TT - svpt0)/(TT - svp3))
      d = 1./(.1*p-ep2*erfac*qsatvp)
      Qs=ep2*qsatvp*d

      dqdta =  (svp2*(svpt0-svp3)/(svp3-TT)**2 )*qsatvp
      dQsdt = ep2*d*dqdta*(1. + erfac*Qs)          

end subroutine MRQSAT_scalar
!======================================














