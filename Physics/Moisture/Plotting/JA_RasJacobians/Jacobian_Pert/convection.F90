subroutine convection_driver(kmax,DT,KCBL,Ts,FRLAND,PREF,U,V,TH,Q,P,PRC3_RAS,PRC3_LSC,RASPARAMS,TH_fixed)

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
 REAL, PARAMETER :: MAPL_EM     = 0.887

!!INPUTS!!
INTEGER :: kmax
REAL :: DT

INTEGER :: KCBL
REAL :: Ts, FRLAND

REAL, DIMENSION(kmax+1) :: PREF

!!INPUT/OUTPUT!!
REAL, DIMENSION(kmax) :: U, V, TH, Q, TH_fixed
REAL, DIMENSION(kmax+1) :: P

REAL :: RASPARAMS(25)

!!OUTPUTS!!
REAL, DIMENSION(kmax) :: PRC3_RAS, PRC3_LSC

!!LOCALS!!
INTEGER :: k

INTEGER, PARAMETER :: IDIM = 1, IRUN = 1


!!!!!!!!!!! RAS !!!!!!!!!!!!

INTEGER :: N_DTL, ICMIN, SEEDRAS(2), point_convective
REAL :: SIGE(kmax+1)

!MODEL VARIABLES READY FOR RAS - PICKS OUT COLUMN AND ALLOWS FOR TRUNCATION ABOVE CONVECTION THROUGH kmax
REAL, DIMENSION(IDIM,kmax) :: U_ras, V_ras, TH_ras, Q_ras, T_ras, TH_ras_fixed, T_fixed
REAL, DIMENSION(IDIM,kmax+1) :: P_ras

REAL, DIMENSION(IDIM,kmax) :: Ph_ras, Pih_ras
REAL, DIMENSION(IDIM,kmax+1) :: Pi_ras
REAL :: ZLE_ras(IDIM,0:kmax), ZLO_ras(IDIM,kmax), ZCBLx_ras(IDIM)
REAL :: QSS_ras(IDIM,kmax), DQS_ras(IDIM,kmax)

!BOUNDARY LAYER TERMS
REAL :: WGT0(IDIM,kmax), WGT1(IDIM,kmax), MXDIAMx(IDIM), TPERT_ras(IDIM), QPERT_ras(IDIM)
REAL, PARAMETER :: CBL_TPERT = 1.0, CBL_QPERT = 0.0, CBL_TPERT_MXOCN = 2.0, CBL_TPERT_MXLND = 4.0
REAL :: QSSFC_ras(IDIM)

INTEGER :: KCBL_ras(IDIM)

!COLUMN OUTPUTS - AS RETURNED FROM RAS, MOSTLY UNUSED
REAL :: CNV_DQLDT(IDIM,kmax), CNV_MF0(IDIM,kmax), CNV_MFD(IDIM,kmax), CNV_MFC(IDIM,kmax)
REAL :: CNV_PRC3(IDIM,kmax), CNV_UPDF(IDIM,kmax), CNV_CVW(IDIM,kmax), CNV_QC(IDIM,kmax)
REAL :: CLCN(IDIM,kmax), HHO(IDIM,kmax), HSO(IDIM,kmax), RASPRCP(IDIM)

!OTHER PARAMETERS
INTEGER, PARAMETER :: ITRCR = 1

REAL, DIMENSION(kmax) :: U_rasout, V_rasout, TH_rasout, Q_rasout


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  LARGE SCALE CLOUD LOOP !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BEGIN RAS PART OF THE CALULATION
SIGE = PREF/PREF(kmax+1) ! this should eventually change (Comment in Moist code)
ICMIN = max(1,count(PREF < PMIN_DET))

MXDIAMx(IDIM) = 0.0

!Pick out sounding for RAS
 Q_ras(IDIM,:) = Q
 TH_ras(IDIM,:) = TH
 U_ras(IDIM,:) = U
 V_ras(IDIM,:) = V
 P_ras(IDIM,:) = P

 TH_ras_fixed(IDIM,:) = TH_fixed

!Compute other required variables
 Ph_ras(IDIM,1:kmax) = 0.5*( P_ras(IDIM,1:kmax) +  P_ras(IDIM,2:kmax+1) ) 

 Pi_ras(IDIM,:) = (P_ras(IDIM,:)/1000.0)**Rkap
 Pih_ras(IDIM,:) = (Ph_ras(IDIM,:)/1000.0)**Rkap

 T_ras(IDIM,:) = TH_ras(IDIM,:)*Pih_ras(IDIM,:)
 T_fixed(IDIM,:) = TH_ras_fixed(IDIM,:)*Pih_ras(IDIM,:)

 !COMPUTE QSS AND DQS - Depend on perturbed quantities.
 call MRQSAT(T_ras(IDIM,:), Ph_ras(IDIM,:), QSS_ras(IDIM,:), DQS_ras(IDIM,:),kmax)

 !SEEDRAS TERMS
 SEEDRAS(1) = 1000000 * ( 100*T_fixed(IDIM,kmax)   - INT( 100*T_fixed(IDIM,kmax) ) )
 SEEDRAS(2) = 1000000 * ( 100*T_fixed(IDIM,kmax-1) - INT( 100*T_fixed(IDIM,kmax-1) ) )
      
 !BOUNDARY LAYER TERMS ETC
 ZLE_ras(IDIM,:) = 0.0
 DO k=kmax,1,-1
    ZLE_ras(IDIM,k-1) = TH_ras(IDIM,k) * (1.+MAPL_VIREPS*Q_ras(IDIM,k))
    ZLO_ras(IDIM,k  ) = ZLE_ras(IDIM,k) + (MAPL_CP/MAPL_GRAV)*( Pi_ras(IDIM,k)-Pih_ras (IDIM,k  ) ) * ZLE_ras(IDIM,k-1)
    ZLE_ras(IDIM,k-1) = ZLO_ras(IDIM,k) + (MAPL_CP/MAPL_GRAV)*( Pih_ras (IDIM,k)-Pi_ras(IDIM,k  ) ) * ZLE_ras(IDIM,k-1)
 endDO

 KCBL_ras(IDIM) = KCBL

 WGT0 = 0.0 
 WGT1 = 0.0

 !Case 6
 WGT0(IDIM, KCBL_ras(IDIM):kmax ) = 1.0 
 WGT1(IDIM, KCBL_ras(IDIM):kmax ) = 1.0       

 ZCBLx_ras = ZLE_ras(IDIM, KCBL_ras(IDIM)-1 )

 ! Cheat by adding a kick to CB temp and q
 QSSFC_ras(IDIM) = 0.0

 TPERT_ras(IDIM)  = CBL_TPERT * ( Ts - ( T_ras(IDIM,kmax)+ MAPL_GRAV*ZLO_ras(IDIM,kmax)/MAPL_CP )  ) 
 QPERT_ras(IDIM)  = CBL_QPERT * ( QSSFC_ras(IDIM) - Q_ras(IDIM,kmax) )          
 TPERT_ras(IDIM)  = MAX( TPERT_ras(IDIM) , 0.0 )
 QPERT_ras(IDIM)  = MAX( QPERT_ras(IDIM) , 0.0 )

 if (FRLAND < 0.1) then
    TPERT_ras(IDIM)  = MIN( TPERT_ras(IDIM) , CBL_TPERT_MXOCN ) ! ocean
 else
    TPERT_ras(IDIM)  = MIN( TPERT_ras(IDIM) , CBL_TPERT_MXLND ) ! land
 endif
   
 N_DTL = KCBL - ICMIN
 
call rase( IDIM, KMAX, ICMIN, DT, MAPL_CP, MAPL_ALHL, MAPL_GRAV, SEEDRAS, SIGE, N_DTL, &
            WGT0(IDIM,:), WGT1(IDIM,:), TPERT_ras(IDIM), QPERT_ras(IDIM), &
            TH_ras(IDIM,:), Q_ras(IDIM,:), U_ras(IDIM,:), V_ras(IDIM,:), &
            QSS_ras(IDIM,:), DQS_ras(IDIM,:),&
            P_ras(IDIM,:), Pi_ras(IDIM,:), &
            CNV_DQLDT(IDIM,:), CNV_MF0(IDIM,:), CNV_MFD(IDIM,:), CNV_MFC(IDIM,:), &
            CNV_PRC3(IDIM,:), CNV_UPDF(IDIM,:), CNV_CVW(IDIM,:), CNV_QC(IDIM,:), &
            HHO(IDIM,:), HSO(IDIM,:), RASPARAMS, POINT_CONVECTIVE )

 U_rasout = U_ras(IDIM,:)
 V_rasout = V_ras(IDIM,:)
 TH_rasout = TH_ras(IDIM,:)
 Q_rasout = Q_ras(IDIM,:)

 PRC3_ras = 100*CNV_PRC3(IDIM,:)*(P_ras(IDIM,2:kmax+1)-P_ras(IDIM,1:kmax)) / DT
      

!SAVE TO INPUT/OUTPUT ARRAYS 
U = U_rasout
V = V_rasout
TH = TH_rasout
Q = Q_rasout

return

end subroutine convection_driver

SUBROUTINE MRQSAT(TT,P,Qs,DQsDT,kmax)

      implicit none

      !Inputs
      integer :: kmax
      real, dimension(kmax) :: TT, P, Qs, DQsDT
      logical :: LDQDT

      !Constants
      real, parameter :: svpt0=273.15, svp1=.6112, svp2=17.67, svp3=29.65, airmw = 28.97, h2omw = 18.01, one = 1.0, ep2 = 0.622
      real, parameter :: esfac = h2omw/airmw, erfac = (one-esfac)/esfac

      !Internals
      real :: dqdta, qsatvp, qsatvp_exc, d
      integer :: kk      

      do kk=30,kmax
            
         qsatvp_exc = svp2*(TT(kk) - svpt0)/(TT(kk) - svp3)
         qsatvp_exc = min(qsatvp_exc, 10.0)

         qsatvp = svp1*exp(qsatvp_exc)

         d = 1./(.1*p(kk)-ep2*erfac*qsatvp)
         Qs(kk)=ep2*qsatvp*d

         dqdta =  (svp2*(svpt0-svp3)/(svp3-TT(kk))**2 )*qsatvp
         dQsdt(kk) = ep2*d*dqdta*(1. + erfac*Qs(kk))
      enddo

end subroutine MRQSAT

SUBROUTINE MRQSAT_SCALAR(TT,P,Qs,DQsDT)

      implicit none

      !Inputs
      real :: TT, P, Qs, DQsDT
      logical :: LDQDT

      !Constants
      real, parameter :: svpt0=273.15, svp1=.6112, svp2=17.67, svp3=29.65, airmw = 28.97, h2omw = 18.01, one = 1.0, ep2 = 0.622
      real, parameter :: esfac = h2omw/airmw, erfac = (one-esfac)/esfac

      !Internals
      real :: dqdta, qsatvp, qsatvp_exc, d     

      qsatvp_exc = svp2*(TT - svpt0)/(TT - svp3)
      qsatvp_exc = min(qsatvp_exc, 10.0)

      qsatvp = svp1*exp(qsatvp_exc)
      d = 1./(.1*p-ep2*erfac*qsatvp)
      Qs=ep2*qsatvp*d

      dqdta =  (svp2*(svpt0-svp3)/(svp3-TT)**2 )*qsatvp
      dQsdt = ep2*d*dqdta*(1. + erfac*Qs)


end subroutine MRQSAT_SCALAR

SUBROUTINE RASE(IDIM,K0,ICMIN,DT,CPO,ALHLO,GRAVO,SEEDRAS,SIGE,N_DTL,WGT0,WGT1,TPERT,QPERT,&
THO,QHO,UHO,VHO,QSS,DQS,PLE,PKE,CLW,FLX,FLXD,FLXC,CNV_PRC3,CNV_UPDFRC,CNV_CVW,CNV_QC,HHO,HSO,RASPARAMS)

IMPLICIT NONE

 !!ARGUMENTS!!
 INTEGER                     ::  IDIM, K0, ICMIN
 REAL                        ::  DT,  CPO, ALHLO, GRAVO
 INTEGER, DIMENSION ( 2 )    ::  SEEDRAS
 REAL, DIMENSION (K0+1)      ::  SIGE
 INTEGER                     ::  N_DTL
 REAL, DIMENSION (K0  )      ::  WGT0,WGT1
 REAL, DIMENSION (IDIM     ) ::  TPERT,QPERT
 REAL, DIMENSION (IDIM,K0  ) ::  THO, QHO, UHO, VHO
 REAL, DIMENSION (IDIM,K0  ) ::  QSS, DQS
 REAL, DIMENSION (IDIM,K0+1) ::  PLE, PKE
 REAL, DIMENSION (IDIM,K0  ) ::  CLW ,FLX, FLXD, FLXC
 REAL, DIMENSION (IDIM,K0  ) ::  CNV_PRC3, CNV_UPDFRC, CNV_CVW, CNV_QC
 REAL, DIMENSION (IDIM,K0  ) ::  HHO, HSO
 REAL, DIMENSION(25)         ::  RASPARAMS
 
 !!LOCALS!!
 INTEGER :: K, L

 REAL, DIMENSION (IDIM) :: MXDIAM
 REAL,  DIMENSION(K0)   :: POI_SV, QOI_SV, UOI_SV, VOI_SV
 REAL,  DIMENSION(K0)   :: POI, QOI, UOI, VOI, DQQ, BET, GAM
 REAL,  DIMENSION(K0)   :: PRH, PRI, GHT, DPT, DPB, PKI, GM1, RNS
 REAL,  DIMENSION(K0)   :: POL, QST, SSL, RMF, HST, QOL, ZOL, HOL
 REAL,  DIMENSION(K0)   :: TEM_PRC3
 REAL,  DIMENSION(K0+1) :: PRJ, PRS, SHT, ZET
 
 !OTHER CONSTANTS
 REAL                   :: WUPDRFT,GRAV,ALHL,CP,CPI,ALHI,GRAVI,CPBG,DDT,AFC,LBCP,OBG
 !RASPARAMS
 REAL                   :: FRICFAC,SHTRG_FAC,CO_AUTO,CLI_CRIT,RASAL1
 REAL                   :: RASAL2,RASNCL,LAMBDA_FAC,LAMBMX_FAC,DIAMMN_MIN
 REAL                   :: FRICLAMBDA,RDTLEXPON,STRAPPING,SDQV2,SDQV3
 REAL                   :: SDQVT1,ACRITFAC,HMINTRIGGER,LLDISAGGXP,PBLFRAC
 REAL                   :: AUTORAMPB,CO_ZDEP,MAXDALLOWED,RHMN,RHMX

 !PARAMETERS
 INTEGER, PARAMETER     :: I = 1

 REAL, PARAMETER        :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
 REAL, PARAMETER        :: RHMAX   = 0.9999

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
 real, parameter :: MAPL_UNDEF  = 1.0e15

 MXDIAM = 0.
 POI_SV = 0.
 QOI_SV = 0.
 UOI_SV = 0.
 VOI_SV = 0.
 POI = 0.
 QOI = 0.
 UOI = 0.
 VOI = 0.
 DQQ = 0.
 BET = 0.
 GAM = 0.
 PRH = 0.
 PRI = 0.
 GHT = 0.
 DPT = 0.
 DPB = 0.
 PKI = 0.
 GM1 = 0.
 RNS = 0.
 POL = 0.
 QST = 0.
 SSL = 0.
 RMF = 0.
 HST = 0.
 QOL = 0.
 ZOL = 0.
 HOL = 0.
 TEM_PRC3 = 0.
 PRJ = 0.
 PRS = 0.
 SHT = 0.
 ZET = 0.
  
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

 RMF = 0.0
 RNS = 0.0

 call HTEST(K,K0,ICMIN,CP,QST,PRJ,PRH,POI,QOI,GRAV,ALHL,CPBG,RHMAX,HOL,HST,SSL,QOL,SHT,ZET,ZOL)

 HHO(I,:) = HOL
 HSO(I,:) = HST

 ! CLOUD LOOP
 call CLOUDE(N_DTL,ICMIN,IDIM,I,K,K0,DT,DDT,DAYLEN,RHMN,RHMX,AUTORAMPB,CO_ZDEP,CP,&
             RHMAX,GRAV,ALHL,ALHI,CPBG,LBCP,RASAL1,RASAL2,ACRITFAC,SDQVT1,SDQV2,SDQV3,CO_AUTO,CLI_CRIT,&
             PBLFRAC,GRAVI,CPI,FRICFAC,FRICLAMBDA,TPERT,QPERT,MXDIAM,PRJ,PRH,PRS,SIGE,&
             DPB,DPT,DQQ,POL,GM1,PKI,BET,GHT,POI,QOI,UOI,VOI,QST,GAM,PRI,RMF,RNS)

 IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN


    DO L=ICMIN,K
       TEM_PRC3(L) = PRI(L)*GRAV
       CNV_PRC3(IDIM,L) = RNS(L)*TEM_PRC3(L)
    ENDDO

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

   call MRQSAT_SCALAR(POI(K)*PRH(K),POL(K),QST(K),DQQ(K))

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
DPB,DPT,DQQ,POL,GM1,PKI,BET,GHT,POI,QOI,UOI,VOI,QST,GAM,PRI,RMF,RNS)
!============================================================================================================================================

IMPLICIT NONE

!!INPUTS!!
INTEGER :: N_DTL, ICL_V(N_DTL), IDIM, I, K, K0, ICMIN
REAL :: DT, DDT, DAYLEN, RHMN, RHMX, AUTORAMPB, CO_ZDEP, CP, RHMAX, GRAV, ALHL, ALHI, CPBG, LBCP
REAL :: RASAL1, RASAL2, ACRITFAC, SDQVT1, SDQV2, SDQV3
REAL :: CO_AUTO, CLI_CRIT, PBLFRAC, GRAVI, CPI, FRICFAC, FRICLAMBDA

REAL, DIMENSION(IDIM) ::  TPERT, QPERT, MXDIAM

REAL, DIMENSION(K0+1) :: PRJ, PRH, PRS, SIGE
REAL, DIMENSION(K0) :: DPB, DPT, DQQ, POL, GM1, PKI, BET, GHT

!!INOUTS!!
REAL, DIMENSION(K0) :: POI, QOI, UOI, VOI
REAL, DIMENSION(K0) :: QST, QOL, SSL, HOL, HST, ZOL, GAM, PRI
REAL, DIMENSION(K0) :: RASAL_1, RASAL_2, RNS
REAL, DIMENSION(K0) :: RMF
REAL, DIMENSION(K0+1) :: ZET, SHT


!!LOCALS!!
REAL :: TRG, WFN, WFN1, WFN2, WFN3, TX2, TX3, QCC, AKM, WLQ
REAL :: LAMBDA_MIN, LAMBDA_MAX, DEEP_FACT, CU_DIAM, WSCALE, WLQ1, WLQ2
REAL :: CLI, C00_X, CLI_CRIT_X, TOKI, DT_LYR, RATE, CVW_X, CLOSS, F3, F4
REAL,  DIMENSION(K0) :: POI_c, QOI_c, QHT, ETA, EHT, TX2a

REAL, DIMENSION(K0) :: CVW, RASAL, HCLD, HCLD1, GMS, GMS1, UCU, VCU
REAL, DIMENSION(K0) :: GMH, RNN, CLL, CLL0, CLLI, CLLB, DLL0, RMFC, RMFD, RMFP, DLLX
REAL :: TEM, TEM1, ALM
REAL, DIMENSION(IDIM,k0) :: ENTLAM


REAL :: UHT, VHT
REAL, DIMENSION(K0) :: F2, BK2, HCC
REAL :: ACR, TE_A

INTEGER :: RC
INTEGER :: L, IC, ICL_C

REAL, PARAMETER :: TE0 = 273.0, TE2 = 200.0

REAL :: JUMP1

!ZERO OUT LOCALS
TRG = 0.0
WFN = 0.0
WFN1 = 0.0
WFN2 = 0.0
WFN3 = 0.0
TX2 = 0.0
TX3 = 0.0
QCC = 0.0
AKM = 0.0
WLQ = 0.0
LAMBDA_MIN = 0.0
LAMBDA_MAX = 0.0
DEEP_FACT = 0.0
CU_DIAM = 0.0
WSCALE = 0.0
WLQ1 = 0.0
WLQ2 = 0.0
CLI = 0.0
C00_X = 0.0
CLI_CRIT_X = 0.0
TOKI = 0.0
DT_LYR = 0.0
RATE = 0.0
CVW_X = 0.0
CLOSS = 0.0
F3 = 0.0
F4 = 0.0
POI_c = 0.0
QOI_c = 0.0
QHT = 0.0
ETA = 0.0
EHT = 0.0
TX2a = 0.0
UHT = 0.0
VHT = 0.0
F2 = 0.0
BK2 = 0.0
HCC = 0.0
ACR = 0.0
TE_A = 0.0
JUMP1 = 0.0
CVW = 0.0
RASAL = 0.0
HCLD = 0.0
HCLD1 = 0.0
GMS = 0.0
GMS1 = 0.0
UCU = 0.0
VCU = 0.0
GMH = 0.0
RNN = 0.0
CLL = 0.0
CLL = 0.0
CLLI = 0.0
CLLB = 0.0
DLL0 = 0.0
RMFC = 0.0
RMFD = 0.0
RMFP = 0.0
DLLX = 0.0
TEM = 0.0
TEM1 = 0.0
ALM = 0.0
ENTLAM = 0.0

 !!!!!!!!!!!!!!! FINDDTLS !!!!!!!!!!!!!!
 do L=1,N_DTL
    ICL_V(L) = K - L 
 enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

JUMP1 =  (SDQV2-1.0) / ( ( TE0-SDQVT1 )**0.333 )

DO ICL_C = 1,N_DTL          

   IC   = ICL_V( ICL_C )
   RC = 10          
 
   UCU(ICMIN:) = 0.    ! This change makes cumulus friction
   VCU(ICMIN:) = 0.    ! correct.

   ALM   = 0.
   TRG   = AMIN1(1.,(QOI(K)/qst(K)-RHMN)/(RHMX-RHMN))
   
   F4  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(IC))/0.2 )  )  ! F4 should ramp from 0 at SIG=AUTORAMPB
         
   !RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
   POI_c = poi
   QOI_c = qoi
   POI_c(K) =  POI_c(K) + TPERT(I)
   QOI_c(K) =  QOI_c(K) + QPERT(I)
   
   ZET(K+1) = 0.
   SHT(K+1) = CP*POI_c(K)*PRJ(K+1)
   DO L=K,IC,-1
      QOL(L)  = AMIN1(qst(L)*RHMAX,QOI_c(L))
      QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
      SSL(L)  = CP*PRJ(L+1)*POI_c(L) + GRAV*ZET(L+1)
      HOL(L)  = SSL(L) + QOL(L)*ALHL
      HST(L)  = SSL(L) + qst(L)*ALHL
      TEM     = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
      ZET(L)  = ZET(L+1) + TEM
      ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
   ENDDO
   
   DO L=IC+1,K
      TEM  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
      SHT(L)  = SSL(L-1) + TEM*(SSL(L)-SSL(L-1)) 
      QHT(L)  = .5*(QOL(L)+QOL(L-1))
   ENDDO
   
   !CALCULATE LAMBDA, ETA, AND WORKFUNCTION
   LAMBDA_MIN = .2/MXDIAM(I)
   LAMBDA_MAX = .2/  200. 
   
   IF (HOL(K) <= HST(IC)) THEN   !CANNOT REACH IC LEVEL
      RC = 1
      cycle
      !====EXIT====>
   endIF
      
   !LAMBDA CALCULATION: MS-A18
   TEM = (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
   DO L=IC+1,K-1
      TEM = TEM + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
   ENDDO
   
   IF (TEM <= 0.0) THEN !NO VALID LAMBDA     
      RC = 2
      cycle
      !====EXIT====>
   endIF
         
   TEM1 = TEM   
   ALM = (HOL(K)-HST(IC)) / TEM
         
   IF (ALM > LAMBDA_MAX) THEN
      RC = 3
      cycle
      !====EXIT====>
   endif   
   
   TOKI = 1.0
   
   IF (ALM < LAMBDA_MIN) THEN
      TOKI = ( ALM/LAMBDA_MIN )**2
   ENDIF
   
   !ETA CALCULATION: MS-A2
   DO L=IC+1,K
      ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
   ENDDO
   
   ETA(IC) = 1.0 + ALM * (ZOL(IC)-ZET(K))
         
   !WORKFUNCTION CALCULATION:  MS-A22
   WFN     = 0.0
   HCC(K)  = HOL(K)
   DO L=K-1,IC+1,-1
      HCC(L) = HCC(L+1) + ( ETA(L) - ETA(L+1) )*HOL(L)
      TEM1    = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
      EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
      WFN  = WFN + (TEM1 - EHT(L)*HST(L))*GAM(L)
   ENDDO
   HCC(IC) = HST(IC)*ETA(IC)
   WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)
   
   !VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
   BK2(K)  = 0.0
   HCLD(K) = HOL(K)
   DO L=K-1,IC,-1
      HCLD1(L) = ( ETA(L+1)*HCLD(L+1)  +  (ETA(L) - ETA(L+1))*HOL(L) ) / ETA(L)
      TEM1     = (HCLD1(L)-HST(L) )*(ZET(L)-ZET(L+1))/ (1.0+LBCP*DQQ(L))
      BK2(L)  = BK2(L+1) + GRAV * MAX(TEM1,0.0) / ( CP*PRJ(L+1)*poi(L) )
      CVW(L) = SQRT(  2.0* MAX( BK2(L) , 0.0 )  )
   ENDDO
   CU_DIAM   = 1000.0
   
   !ALPHA CALCULATION 
   IF ( ZET(IC) <  2000. ) THEN
      RASAL(IC) = RASAL1
   ENDIF

   IF ( ZET(IC)  >= 2000. ) THEN 
      RASAL(IC) = RASAL1 + (RASAL2-RASAL1)*(ZET(IC) - 2000.)/8000.
   ENDIF
   RASAL_1(IC) = MIN( RASAL(IC) , 1.0e5 )
       
   RASAL_2(IC) = DT / RASAL_1(IC)
   
   CVW(IC:K) = MAX( CVW(IC:K), 1.00 )

   !TEST FOR CRITICAL WORK FUNCTION
   CALL ACRITN(POL(IC),PRS(K),ACR,ACRITFAC)
   
   IF (WFN <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION
      RC = 4
      cycle
      !====EXIT====>
   endIF
   
   IF (RC == 10) THEN

      !CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
      WLQ      = QOL(K)
      UHT = uoi(K)
      VHT = voi(K)
      RNN(K)  = 0.
      CLL0(K) = 0.
      
      
      DO L=K-1,IC,-1
         TEM1  = ETA(L) - ETA(L+1)
         WLQ   = WLQ + TEM1 * QOL(L)
         UHT   = UHT + TEM1 * uoi(L)
         VHT   = VHT + TEM1 * voi(L)
         
         ! How much condensate (CLI) is present here? 
         IF (L>IC) THEN
            TX2 = 0.5*(qst(L)+qst(L-1))*ETA(L)
            TX3 = 0.5*(HST(L)+HST(L-1))*ETA(L)
            QCC = TX2 + GM1(L)*(HCC(L)-TX3)
            CLL0(L) = (WLQ-QCC)
         ELSE
            CLL0(L)   = (WLQ-qst(IC)*ETA(IC))
         ENDIF
         CLL0(L)    = MAX(CLL0(L), 0.00)
         
         CLI  = CLL0(L) / ETA(L)  !Condensate (kg/kg)
         TE_A = poi(L)*PRH(L)     !Temperature (K)
         
         !!!!!!!!!!!!!!!!!! SUNDQ3 CALL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( TE_A .GE. TE0 ) THEN 
            F2(L)   = 1.0
         ENDIF
         IF ( ( TE_A .GE. SDQVT1 ) .AND. ( TE_A .LT. TE0 ) )THEN 
            F2(L)   = 1.0 + JUMP1 * (( TE0 - TE_A )**0.3333)
         ENDIF
         IF ( TE_A .LT. SDQVT1 ) THEN 
            F2(L)   = SDQV2 + (SDQV3-SDQV2)*(SDQVT1-TE_A)/(SDQVT1-TE2)
         ENDIF
         IF ( F2(L) .GT. 27.0 ) THEN
            F2(L) = 27.0
         ENDIF
         F3 = 1.0
         !!!!!!!!!!!!!!!!!! END SUNDQ3 CALL !!!!!!!!!!!!!!!!!!!!!!!!
   
         C00_X  =  CO_AUTO * F2(L) * F3  * F4  
         CLI_CRIT_X =  CLI_CRIT / ( F2(L)*F3 )
      
         RATE = C00_X * ( 1.0 - EXP( -(CLI)**2 / CLI_CRIT_X**2 ) )
         
         CVW_X = MAX( CVW(L), 1.00 ) !Floor convective velocity since we don't really trust it at low values.

         DT_LYR = ( ZET(L)-ZET(L+1) )/CVW_X  !L.H.S. DT_LYR => time in layer (L,L+1)
         
         CLOSS = CLL0(L) * RATE * DT_LYR
         CLOSS = MIN( CLOSS , CLL0(L) )

         CLL0(L) = CLL0(L) - CLOSS
         DLL0(L) = CLOSS
         
         IF (CLOSS>0.) then
            WLQ = WLQ - CLOSS
            RNN(L) = CLOSS 
         else
            RNN(L) = 0.
         ENDIF
         
      ENDDO
      
      WLQ = WLQ - qst(IC)*ETA(IC)

      !CALCULATE GAMMAS AND KERNEL
      GMS(K) = (SHT(K)-SSL(K))*PRI(K)                ! MS-A30 (W/O GRAV)
      GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*ALHL  ! MS-A31 (W/O GRAV)
      AKM    = GMH(K)*GAM(K-1)*DPB(K-1)              ! MS-A37 (W/O GRAV)
      
      TX2a(K-1)     = GMH(K)
      DO L=K-1,IC+1,-1
         GMS(L) = ( ETA(L)*(SHT(L)-SSL(L  )) + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
         GMH(L) = GMS(L) + ( ETA(L)*(QHT(L) - QOL(L)) + ETA(L+1)*(QOL(L) - QHT(L+1)) )*ALHL*PRI(L)
         TX2a(L-1) = TX2a(L) + (ETA(L) - ETA(L+1)) * GMH(L)
         !AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2a(L-1)*GHT(L)
      ENDDO
      
      DO L=K-1,IC+1,-1
         AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2a(L-1)*GHT(L)
      ENDDO

      GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
      AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)
      
      GMH(IC) = GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL + ETA(IC)*(HST(IC)-HOL(IC)) )*PRI(IC)
      
   endIF

   !CLOUD BASE MASS FLUX
   IF (AKM >= 0.0 .OR. WLQ < 0.0) THEN
      RC = 5
      cycle
      !====EXIT====>
   endIF
   
   IF (RC == 10) THEN

      WFN1 = -(WFN-ACR)/AKM ! MS-A39 MASS-FLUX IN Pa/step
      WFN1 = MIN( ( RASAL_2(IC)*TRG*TOKI )*WFN1, (PRS(K+1)-PRS(K))*(100.*PBLFRAC) )
     
      !CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
      
      TEM1     = WFN1*GRAVI
      CLL (IC) = CLL (IC) + WLQ*TEM1       ! (kg/m^2/step)
      RMF (IC) = RMF (IC) + TEM1           ! (kg/m^2/step)

          
      !THETA AND Q CHANGE DUE TO CLOUD TYPE IC
      !DO L=IC,K
         RNS(IC:K) = RNS(IC:K) + RNN(IC:K)*TEM1 ! (kg/m^2/step)
         GMH(IC:K) = GMH(IC:K) * WFN1
         GMS1(IC:K) = GMS(IC:K) * WFN1
         qoi(IC:K) = qoi(IC:K) + (GMH(IC:K) - GMS1(IC:K)) * ALHI
         poi(IC:K) = poi(IC:K) + GMS1(IC:K)*PKI(IC:K)*CPI
         qst(IC:K) = qst(IC:K) + GMS1(IC:K)*BET(IC:K)*CPI
      !ENDDO
      
            
      !CUMULUS FRICTION
      WFN2 = 0.5*WFN1*FRICFAC*EXP( -ALM / FRICLAMBDA )
      TEM1 = WFN2*PRI(K)
      
      UCU(K)  = UCU(K) + TEM1 * (uoi(K-1) - uoi(K))      
      VCU(K)  = VCU(K) + TEM1 * (voi(K-1) - voi(K))
      
      DO L=K-1,IC+1,-1
         TEM1    = WFN2*PRI(L)
         UCU(L) = UCU(L) + TEM1 * ( (uoi(L-1) - uoi(L)) * ETA(L) &
                                 + (uoi(L) - uoi(L+1)) * ETA(L+1) )
         VCU(L) = VCU(L) + TEM1 * ( (voi(L-1) - voi(L)) * ETA(L) &
                                 + (voi(L) - voi(L+1)) * ETA(L+1) )
      ENDDO

      TEM1  = WFN2*PRI(IC)

      UCU(IC) = UCU(IC) + (2.*(UHT - uoi(IC)*(ETA(IC)-ETA(IC+1))) &
              - (uoi(IC)+uoi(IC+1))*ETA(IC+1))*TEM1
      VCU(IC) = VCU(IC) + (2.*(VHT - voi(IC)*(ETA(IC)-ETA(IC+1))) &
              - (voi(IC)+voi(IC+1))*ETA(IC+1))*TEM1


      uoi = uoi
      voi = voi

      DO L=IC,K
         uoi(L) = uoi(L) + UCU(L)
         voi(L) = voi(L) + VCU(L)
      ENDDO

      RC = 0


   ENDIF
                  
   ENTLAM(I,IC) = ALM

ENDDO

POI = poi
QOI = qoi
UOI = uoi
VOI = voi
QST = qst

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





















