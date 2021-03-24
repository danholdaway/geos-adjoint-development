subroutine convection_driver( IM, JM, LM, PREF, DT_MOIST, FRLAND, KCBL, Ts, RASPARAMS, U, V, TH, Q, PLE, &
                        CNV_DQLDT,CNV_MFD,CNV_PRC3,CNV_UPDF,ESTBLX )


IMPLICIT NONE

!!INPUTS!!
INTEGER, INTENT(IN) :: IM, JM, LM
real(8), INTENT(IN) :: DT_MOIST

real(8), DIMENSION(0:LM), INTENT(IN) :: PREF

INTEGER, DIMENSION(IM,JM), INTENT(IN) :: KCBL
real(8), DIMENSION(IM,JM), INTENT(IN) :: FRLAND

real(8), INTENT(IN) :: RASPARAMS(25)

real(8), DIMENSION(IM,JM), INTENT(IN) :: Ts
real(8), DIMENSION(IM,JM,0:LM), INTENT(IN) :: PLE

!!INOUTS!!
real(8), DIMENSION(IM,JM,LM), INTENT(INOUT) :: U, V, TH, Q

!!OUTPUTS!!
real(8), DIMENSION(IM,JM,LM), INTENT(OUT) :: CNV_DQLDT,CNV_MFD,CNV_PRC3,CNV_UPDF


!!LOCALS!!
real(8), parameter :: CONS_RUNIV  = 8314.3, CONS_KAPPA = 2.0/7.0, CONS_AIRMW  = 28.97, CONS_H2OMW  = 18.01, CONS_GRAV   = 9.80
real(8), parameter :: CONS_RGAS = CONS_RUNIV/CONS_AIRMW, CONS_CP = CONS_RGAS/CONS_KAPPA, CONS_VIREPS = CONS_AIRMW/CONS_H2OMW-1.0

integer :: i,j,k,l,ii,jj, IDIM, K0, ICMIN

real(8), parameter :: PMIN_DET = 3000.0, CBL_TPERT = 1.0, CBL_QPERT = 0.0
real(8), parameter :: AUTOC_CN_OCN = 2.5e-3, AUTOC_CN_LAND = 2.5e-3
real(8), parameter :: CBL_TPERT_MXOCN = 2.0, CBL_TPERT_MXLND = 4.0

real(8), dimension(0:LM) :: SIGE

integer, dimension(im,jm)  :: SEEDRAS

real(8), dimension(im,jm) :: qssfc, TPERT, QPERT, CO_AUTO

real(8), dimension(im,jm,0:lm) :: cnv_ple, pke, zle

real(8), dimension(im,jm,lm) :: plo, pk, temp, dqs, qss, zlo, wgt0, wgt1

!QSATVP TABLE
real(8) :: ESTBLX(:)

!ZERO OUTPUTS
CNV_DQLDT = 0.0
CNV_MFD = 0.0
CNV_PRC3 = 0.0
CNV_UPDF = 0.0

SIGE = 0.0
SEEDRAS = 0
qssfc = 0.0
TPERT = 0.0
QPERT = 0.0
cnv_ple = 0.0
pke = 0.0
zle = 0.0
plo = 0.0
pk = 0.0
temp = 0.0
dqs = 0.0
qss = 0.0
zlo = 0.0
wgt0 = 0.0
wgt1 = 0.0

IDIM = 1

K0 = LM
ICMIN = 0
DO L=0,LM
   IF (PREF(L) < PMIN_DET) THEN
      ICMIN = ICMIN + 1
   endif
endDO

CNV_PLE  = PLE*.01

PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
PKE      = (CNV_PLE/1000.)**(CONS_RGAS/CONS_CP)
PK       = (PLO/1000.)**(CONS_RGAS/CONS_CP)
     
TEMP     = TH*PK

call DQSAT_sub(DQS,QSS,TEMP,PLO,im,jm,lm,ESTBLX)

QSSFC    = 0.0 !GEOS_QSAT( TS , PLE(:,:,LM) ) - Not needed currently as CBL_qpert = 0

ZLE(:,:,LM) = 0.
do L=LM,1,-1
   ZLE(:,:,L-1) = TH (:,:,L) * (1.+CONS_VIREPS*Q(:,:,L))
   ZLO(:,:,L  ) = ZLE(:,:,L) + (CONS_CP/CONS_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
   ZLE(:,:,L-1) = ZLO(:,:,L) + (CONS_CP/CONS_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
end do

TPERT  = CBL_TPERT * ( TS - ( TEMP(:,:,LM)+ CONS_GRAV*ZLO(:,:,LM)/CONS_CP )  ) 
QPERT  = CBL_QPERT * ( QSSFC - Q(:,:,LM) )          
TPERT  = MAX( TPERT , 0.0 )

where (FRLAND<0.1) 
   TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
elsewhere
   TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
end where

where (FRLAND<0.1) 
   CO_AUTO = AUTOC_CN_OCN   ! ocean value
elsewhere
   CO_AUTO = AUTOC_CN_LAND  ! land value
end where

SEEDRAS(:,:) = INT(1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) ))

SIGE = PREF/PREF(LM)

DO I = 1,IM
   DO J = 1,JM

      WGT0(I,J,:)            = 0.
      WGT0(I,J,KCBL(I,J):K0) = 1.0 
      WGT1(I,J,:)            = 0.
      WGT1(I,J,KCBL(I,J):K0) = 1.0 

   endDO
endDO

call RASE     (  IDIM                 , &
                 K0                   , &
                 ICMIN                , &
                 DT_MOIST             , &  !!where? see below.
                 SEEDRAS              , &
                 SIGE                 , &
! inputs for CBL
                 KCBL                 , &
                 WGT0                 , &
                 WGT1                 , &
                 TPERT                , &
                 QPERT                , &
! inputs
                 TH                  , &
                 Q                   , & 
                 U                   , &
                 V                   , &
                 QSS                  , & 
                 DQS                  , &
! Pass in CO_AUTO
                 CO_AUTO              , &
!
                 CNV_PLE              , &
                 PKE                  , &
! outputs
                 CNV_DQLDT            , &   ! -> progno_cloud
                 CNV_MFD              , &   ! -> progno_cloud
                 CNV_PRC3             , &   ! -> progno_cloud 
                 CNV_UPDF             , &   ! -> progno_cloud
! params
                 RASPARAMS            , &
                 ESTBLX                )


end subroutine convection_driver


SUBROUTINE RASE(IDIM, K0, ICMIN, DT ,                      &
                SEEDRAS,SIGE,                          &
                KCBL,WGT0,WGT1,CTPERT,CQPERT,          &
                THO, QHO, UHO, VHO,                              & 
                CQSS, CDQS,                                        &
                CCO_AUTO,                                         &
                CPLE, CPKE, &
                CLW, FLXD, CCNV_PRC3,CNV_UPDFRC,                                      &
                RASPARAMS,ESTBLX                                        )

 IMPLICIT NONE

 !ARGUMENTS
 INTEGER,                     INTENT(IN   ) ::  IDIM, K0, ICMIN
 real(8), DIMENSION (IDIM,K0  ), INTENT(INOUT) ::  THO, QHO, UHO, VHO
 real(8), DIMENSION (IDIM,K0+1), INTENT(IN   ) ::  CPLE, CPKE
 real(8), DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  CQSS, CDQS
 real(8), DIMENSION (     K0+1), INTENT(IN   ) ::  SIGE
 real(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CLW , FLXD
 real(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CCNV_PRC3
 real(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CNV_UPDFRC 
 real(8),                        INTENT(IN   ) ::  DT
 INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  SEEDRAS
 INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  KCBL
 real(8), DIMENSION (IDIM     ), INTENT(IN   ) ::  CTPERT,CQPERT
 real(8), DIMENSION (IDIM     ), INTENT(IN   ) ::  CCO_AUTO
 real(8), DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  WGT0,WGT1
 real(8), DIMENSION(:),          INTENT(IN   ) ::  RASPARAMS
 real(8), DIMENSION(:),          INTENT(IN   ) ::  ESTBLX

 !RASPARAMS REDEFINED
 real(8) :: FRICFAC,CLI_CRIT,RASAL1,RASAL2,FRICLAMBDA,SDQV2,SDQV3,SDQVT1, MXDIAM
 real(8) :: ACRITFAC,HMINTRIGGER,PBLFRAC,AUTORAMPB,CO_ZDEP,MAXDALLOWED,RHMN,RHMX

 !GLOBAL CONSTANTS
 real(8), PARAMETER :: CONS_RUNIV  = 8314.3, CONS_GRAV = 9.80, CONS_KAPPA = 2.0/7.0, CONS_AIRMW  = 28.97, CONS_ALHL = 2.4665E6
 real(8), PARAMETER :: CONS_RGAS = CONS_RUNIV/CONS_AIRMW, CONS_CP = CONS_RGAS/CONS_KAPPA

 real(8), PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0, RHMAX = 0.9999
 real(8) GRAV, CP, ALHL, CPBG, ALHI, CPI, GRAVI, DDT, LBCP 

 !LOCALS
 INTEGER I, K, IC, L

 real(8), DIMENSION(K0) :: POI_SV, QOI_SV, UOI_SV, VOI_SV
 real(8), DIMENSION(K0) :: POI, QOI, UOI, VOI,  DQQ, BET, GAM, CLL
 real(8), DIMENSION(K0) :: POI_c, QOI_c
 real(8), DIMENSION(K0) :: PRH,  PRI,  GHT, DPT, DPB, PKI
 real(8), DIMENSION(K0) :: TCU, QCU, RNS, POL
 real(8), DIMENSION(K0) :: QST, SSL,  RMF, RMFC, RMFP
 real(8), DIMENSION(K0) :: GM1, RMFD
 real(8), DIMENSION(K0) :: HOL, HST, QOL, ZOL, CLL0,CLLI,CLLB
 real(8), DIMENSION(K0) :: CVW, UPDFRC
 real(8), DIMENSION(K0) :: UPDFRP,DLL0,DLLX

 real(8), DIMENSION(K0+1) :: PRJ, PRS, SHT ,ZET

 real(8) :: TEMa, TEMb, TEMc

 !QSATVP LOOKUP TABLE
 integer, parameter :: DEGSUBS    =  100
 real(8),    parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
       
 !ZERO OUTPUTS
 CLW = 0.0
 FLXD = 0.0
 CCNV_PRC3 = 0.0
 CNV_UPDFRC = 0.0
      
 !Zero locals, taf linearisation safety measure
 POI_SV = 0.0
 QOI_SV = 0.0
 UOI_SV = 0.0
 VOI_SV = 0.0
 POI = 0.0
 QOI = 0.0
 UOI = 0.0
 VOI = 0.0
 DQQ = 0.0
 BET = 0.0
 GAM = 0.0
 CLL = 0.0
 POI_c = 0.0
 QOI_c = 0.0
 PRH = 0.0
 PRI = 0.0
 GHT = 0.0
 DPT = 0.0
 DPB = 0.0
 PKI = 0.0
 TCU = 0.0
 QCU = 0.0
 RNS = 0.0
 POL = 0.0
 QST = 0.0
 SSL = 0.0
 RMF = 0.0
 RMFC = 0.0
 RMFP = 0.0
 GM1 = 0.0
 RMFD = 0.0
 HOL = 0.0
 HST = 0.0
 QOL = 0.0
 ZOL = 0.0
 CLL0 = 0.0
 CLLI = 0.0
 CLLB = 0.0
 CVW = 0.0
 UPDFRC = 0.0
 UPDFRP = 0.0
 DLL0 = 0.0
 DLLX = 0.0
 PRJ  = 0.0
 PRS = 0.0
 SHT = 0.0
 ZET = 0.0 
 TEMa = 0.0
 TEMb = 0.0
 TEMc = 0.0

 !ASSIGN RASPARAMS NEEDED
 FRICFAC      = RASPARAMS(1)     !  ---  1
 CLI_CRIT     = RASPARAMS(4)     !  ---  4
 RASAL1       = RASPARAMS(5)     !  ---  5
 RASAL2       = RASPARAMS(6)     !  ---  6
 FRICLAMBDA   = RASPARAMS(11)    !  --- 11
 SDQV2        = RASPARAMS(14)    !  --- 14
 SDQV3        = RASPARAMS(15)    !  --- 15
 SDQVT1       = RASPARAMS(16)    !  --- 16
 ACRITFAC     = RASPARAMS(17)    !  --- 17
 HMINTRIGGER  = RASPARAMS(18)    !  --- 18
 PBLFRAC      = RASPARAMS(20)    !  --- 20
 AUTORAMPB    = RASPARAMS(21)    !  --- 21
 CO_ZDEP      = RASPARAMS(22)    !  --- 22
 MAXDALLOWED  = RASPARAMS(23)    !  --- 23
 RHMN         = RASPARAMS(24)    !  --- 24
 RHMX         = RASPARAMS(25)    !  --- 25

 !ASIGN RASE CONSTANTS
 GRAV  = CONS_GRAV
 ALHL  = CONS_ALHL
 CP    = CONS_CP
 CPI   = 1.0/CP      
 ALHI  = 1.0/ALHL
 GRAVI = 1.0/GRAV
 CPBG  = CP*GRAVI
 DDT   = DAYLEN/DT
 LBCP  = ALHL*CPI

 I = 1

 !!CALL FINDBASE
 K = KCBL(I)
 
 IF ( K > 0 ) THEN 

    CALL STRAP0(I,K,K0,IDIM,ICMIN,CPKE,CPLE,UHO,VHO,THO,QHO,CDQS,CQSS,ESTBLX,SEEDRAS,WGT0,MXDIAM,&
                POI,QOI,UOI,VOI,PRJ,PRS,QST,DQQ,CONS_RGAS,CONS_CP,ONEPKAP,&
                POL,PRH,PKI,DPT,DPB,PRI,MAXDALLOWED,BET,GAM,LBCP,GHT,GM1,ALHL,TCU,QCU,&
                RNS,CLL,RMF,RMFD,RMFC,RMFP,CLL0,DLL0,DLLX,CLLI,CLLB,CVW,UPDFRC,UPDFRP,&
                POI_SV,QOI_SV,UOI_SV,VOI_SV)

    !HTEST ROUTINE
    HOL = 0.0
    ZET(K+1) = 0
    SHT(K+1) = CP*POI(K)*PRJ(K+1)
    DO L = K,ICMIN,-1
       QOL(L)  = MIN(QST(L)*RHMAX,QOI(L))
       QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
       SSL(L)  = CP*PRJ(L+1)*POI(L) + GRAV*ZET(L+1)
       HOL(L)  = SSL(L) + QOL(L)*ALHL
       HST(L)  = SSL(L) + QST(L)*ALHL
       TEMa    = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
       ZET(L)  = ZET(L+1) + TEMa
       ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
    ENDDO
    !END HTEST ROUTINE


    !MAIN CLOUD LOOP
    CALL CLOUDE(I,K,K0,IDIM,ICMIN,DT,DDT,DAYLEN,RHMN,RHMX,AUTORAMPB,     &
                CO_ZDEP,CP,CPI,RHMAX,GRAV,GRAVI,ALHI,ALHL,CPBG,LBCP,RASAL1,RASAL2,   &
                ACRITFAC,SDQV2,SDQV3,SDQVT1,CLI_CRIT,PBLFRAC,FRICFAC,FRICLAMBDA,     &
                QOI,QST,POI,UOI,VOI,SIGE,POI_C,QOI_C,CTPERT,CQPERT,MXDIAM,CCO_AUTO,  &
                ZET,SHT,QOL,SSL,HOL,HST,ZOL,PRJ,PRH,DPB,DPT,BET,GAM,CVW,DQQ,POL,PRS, &
                CLL0,DLL0,GM1,PRI,GHT,PKI,CLL,RMF,RMFD,RMFP,RMFC,DLLX,UPDFRP,UPDFRC, &
                CLLI,CLLB,RNS)
      

    IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

       !RNEVP ROUTINE
       ZET(K+1) = 0
       DO L = K,ICMIN,-1
          TEMb     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
          ZET(L)  = ZET(L+1) + TEMb
       ENDDO

       DO L = ICMIN,K
          TEMc    = PRI(L)*GRAV
          CCNV_PRC3(I,L) = RNS(L)*TEMc     
       ENDDO
       !END RNEVP ROUTINE

       CALL STRAP1(I,K,K0,IDIM,ICMIN,&
                   THO,QHO,UHO,VHO,CNV_UPDFRC,CPLE,FLXD,CLW,WGT1,&
                   POI,QOI,UOI,VOI,UPDFRC,CVW,CLLB,PRS,POI_SV,QOI_SV,UOI_SV,VOI_SV,&
                   RMF,RMFD,RMFC,CLL,DDT,DAYLEN )

    ELSE
   
       FLXD(I,:) = 0.
       CLW (I,:) = 0.

    ENDIF 

 ELSE 
     
    FLXD(I,:) = 0.
    CLW (I,:) = 0.

 ENDIF

 RETURN

END SUBROUTINE RASE


SUBROUTINE CLOUDE(I,K,K0,IDIM,ICMIN,DT,DDT,DAYLEN,RHMN,RHMX,AUTORAMPB,     &
                  CO_ZDEP,CP,CPI,RHMAX,GRAV,GRAVI,ALHI,ALHL,CPBG,LBCP,RASAL1,RASAL2,   &
                  ACRITFAC,SDQV2,SDQV3,SDQVT1,CLI_CRIT,PBLFRAC,FRICFAC,FRICLAMBDA,     &
                  QOI,QST,POI,UOI,VOI,SIGE,POI_C,QOI_C,CTPERT,CQPERT,MXDIAM,CCO_AUTO,  &
                  ZET,SHT,QOL,SSL,HOL,HST,ZOL,PRJ,PRH,DPB,DPT,BET,GAM,CVW,DQQ,POL,PRS, &
                  CLL0,DLL0,GM1,PRI,GHT,PKI,CLL,RMF,RMFD,RMFP,RMFC,DLLX,UPDFRP,UPDFRC, &
                  CLLI,CLLB,RNS)

 IMPLICIT NONE

 !!GLOBALS!!
 INTEGER :: I, K, K0, IDIM, ICMIN

 real(8) :: DT, DDT, DAYLEN, RHMN, RHMX, AUTORAMPB, CO_ZDEP, CP, CPI, RHMAX, GRAV, GRAVI, ALHI, ALHL, CPBG, LBCP
 real(8) :: RASAL1, RASAL2, ACRITFAC,SDQV2,SDQV3,SDQVT1,CLI_CRIT,PBLFRAC,FRICFAC,FRICLAMBDA, MXDIAM

 real(8), DIMENSION(IDIM) :: CTPERT, CQPERT, CCO_AUTO
 real(8), DIMENSION(K0) :: QOI, QST, POI, UOI, VOI, POI_C, QOI_C
 real(8), DIMENSION(K0) :: QOL, SSL, HOL, HST, ZOL, PRH, DPB, DPT, BET, GAM, CVW
 real(8), DIMENSION(K0) :: DQQ, POL, CLL0, DLL0, GM1, PRI, GHT, PKI
 real(8), DIMENSION(K0) :: CLL,RMF,RMFD,RMFP,RMFC,DLLX,UPDFRP,UPDFRC,CLLI,CLLB,RNS
 real(8), DIMENSION(K0+1) :: SIGE, ZET, SHT, PRJ, PRS

 !!LOCALS!!
 INTEGER :: IC, L, RC, ICL, ICL_C, N_DTL

 real(8) :: TE0, TE1, TE2, JUMP1

 real(8) :: TEM1, TEM2, TEM3, TEM4, TEM5, TEM6, TEM7, TEM8, TEM9, TEM10
 real(8) :: ALM, TRG, LAMBDA_MIN, LAMBDA_MAX
 real(8) :: WFN, WFN1, WFN2, ACR, AKM, WLQ, UHT, VHT, TX2, TX3, QCC

 real(8) :: CLI , TE_A, C00_X, CLI_CRIT_X, TOKI
 real(8) :: DT_LYR, RATE, CVW_X, CLOSS, F2, F4, F5

 real(8), DIMENSION(K0) :: ETA, HCC, EHT, BK2, HCLD, RASAL, RNN, UCU, VCU
 real(8), DIMENSION(K0) :: GMS, GMS1, GMH
 real(8), DIMENSION(K0+1) :: QHT

 !LOCALS TO ZERO
 CLI = 0.0
 TE_A = 0.0
 C00_X = 0.0
 CLI_CRIT_X = 0.0
 TOKI = 0.0
 DT_LYR = 0.0
 RATE = 0.0
 CVW_X = 0.0
 CLOSS = 0.0
 F2 = 0.0
 F4 = 0.0
 F5 = 0.0
 ALM = 0.0
 TRG = 0.0
 TEM1 = 0.0
 TEM2 = 0.0
 TEM3 = 0.0
 TEM4 = 0.0
 TEM5 = 0.0
 TEM6 = 0.0
 TEM7 = 0.0
 TEM8 = 0.0
 TEM9 = 0.0
 TEM10 = 0.0
 LAMBDA_MIN = 0.0
 LAMBDA_MAX = 0.0
 WFN = 0.0
 WFN1 = 0.0
 WFN2 = 0.0
 ACR = 0.0
 AKM = 0.0
 WLQ = 0.0
 UHT = 0.0
 VHT = 0.0
 TX2 = 0.0
 TX3 = 0.0
 QCC = 0.0
 ETA = 0.0
 HCC = 0.0
 EHT = 0.0
 BK2 = 0.0
 HCLD = 0.0
 RASAL = 0.0
 RNN = 0.0
 UCU = 0.0
 VCU = 0.0
 QHT = 0.0
 GMS = 0.0
 GMS1 = 0.0
 GMH = 0.0

 TE0=273.0
 TE1=SDQVT1
 TE2=200.
 JUMP1=  (SDQV2-1.0) / ( ( TE0-TE1 )**0.333 ) 

 N_DTL = K - ICMIN
 
 !CLOUD LOOP
 DO ICL_C = 1,N_DTL

    ICL   = K - ICL_C !ICL_V( ICL_C )
 
    UCU(ICMIN:) = 0.    ! This change makes cumulus friction
    VCU(ICMIN:) = 0.    ! correct.

    IF ( ICL > ICMIN ) THEN

       IC = ICL

       RC = 10
      
       ALM   = 0.
       if (abs(QST(K)) .gt. 0.0) then
          TRG   = MIN(1.,(QOI(K)/QST(K)-RHMN)/(RHMX-RHMN))
       endif
      
       F4  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(IC))/0.2 )  )  ! F4 should ramp from 0 at SIG=AUTORAMPB
                                                                    ! to 1 at SIG=AUTORAMPB-0.2
       !if ( SIGE(IC) >= 0.5 ) then
       !   F5 = 1.0
       !else
       !   F5 = 1.0 - 2.*CO_ZDEP *( 0.5 - SIGE(IC) )
       !   F5 = MAX( F5 , 0.0 )
       !endif
      
       IF (TRG <= 1.0E-5) THEN    ! TRIGGER  =========>>      
          RC = 1
          !RETURN
       ENDIF

       IF (RC == 10) THEN

          !RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
          POI_c = POI
          QOI_c = QOI
          POI_c(K) =  POI_c(K) + CTPERT(I)
          QOI_c(K) =  QOI_c(K) + CQPERT(I)
   
          ZET(K+1) = 0.
          SHT(K+1) = CP*POI_c(K)*PRJ(K+1)
          DO L = K,IC,-1   
             QOL(L)  = MIN(QST(L)*RHMAX,QOI_c(L))
             QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
             SSL(L)  = CP*PRJ(L+1)*POI_c(L) + GRAV*ZET(L+1)
             HOL(L)  = SSL(L) + QOL(L)*ALHL
             HST(L)  = SSL(L) + QST(L)*ALHL
             TEM1    = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
             ZET(L)  = ZET(L+1) + TEM1
             ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG   
             ENDDO
   
          DO L = IC+1,K
             TEM2  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
             SHT(L)  = SSL(L-1) + TEM2*(SSL(L)-SSL(L-1)) 
             QHT(L)  = .5*(QOL(L)+QOL(L-1))
          ENDDO
      
          !CALCULATE LAMBDA, ETA, AND WORKFUNCTION
          LAMBDA_MIN = .2/MXDIAM
          LAMBDA_MAX = .2/  200. 
         
          IF (HOL(K) <= HST(IC)) THEN   ! CANNOT REACH IC LEVEL  ======>>
             RC = 2
             !RETURN
             !====EXIT====>
          ENDIF

          IF (RC == 10) THEN

             !LAMBDA CALCULATION: MS-A18
      
             TEM3 = (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
             DO L = IC+1,K-1
                TEM3 = TEM3 + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
             ENDDO

             IF (TEM3 <= 0.0) THEN         ! NO VALID LAMBDA  ============>>
                RC = 3
                !RETURN
                !====EXIT====>
             ENDIF

             IF (RC == 10) THEN

                if (abs(TEM3) .gt. 0.0) then !Linearisation security
                   ALM  = (HOL(K)-HST(IC)) / TEM3
                endif

                IF (ALM > LAMBDA_MAX) THEN
                   RC = 4
                   !RETURN
                   !====EXIT====>
                ENDIF

                IF (RC == 10) THEN

                   TOKI=1.0
      
                   IF (ALM < LAMBDA_MIN) THEN      
                      TOKI = ( ALM/LAMBDA_MIN )**2
                   ENDIF
      
                   !ETA CALCULATION: MS-A2
                   DO L = IC+1,K
                      ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
                   ENDDO
                   ETA(IC) = 1.0 + ALM * (ZOL(IC)-ZET(K))
      
                   !WORKFUNCTION CALCULATION:  MS-A22
                   WFN     = 0.0
                   HCC(K)  = HOL(K)
                   DO L = K-1,IC+1,-1
                      HCC(L) = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
                      TEM4   = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
                      EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
                      WFN    = WFN + (TEM4 - EHT(L)*HST(L))*GAM(L)
                   ENDDO
                   HCC(IC) = HST(IC)*ETA(IC)
                   WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)
      
                   !VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)


                   HCLD(K)  = HOL(K)
                   DO L = K-1,IC,-1
      
                      if (abs(ETA(L)) .gt. 0.0) then !Linearisation security
                         HCLD(L) = ( ETA(L+1)*HCLD(L+1) + (ETA(L) - ETA(L+1))*HOL(L) ) / (ETA(L) + 1e-15)
                      endif

                   ENDDO

                   BK2(K)   = 0.0
                   DO L = K-1,IC,-1

                      if (abs(1.0+LBCP*DQQ(L)) .gt. 0.0) then !Linearisation security
                         TEM5     = (HCLD(L)-HST(L) )*(ZET(L)-ZET(L+1))/ (1.0+LBCP*DQQ(L)+1e-15)
                         if ( abs( CP*PRJ(L+1)*POI(L) ) .gt. 0.0) then !Linearisation security
                            BK2(L)  = BK2(L+1) + GRAV * MAX(TEM5,0.0) / ( CP*PRJ(L+1)*POI(L) )
                         endif
                      endif

                      CVW(L) = SQRT(  2.0* MAX( BK2(L) , 0.0 )  )

                   ENDDO

                   !ALPHA CALCULATION 
                   RASAL(IC) = DT / RASAL1
                   IF ( ZET(IC)  >= 2000. ) THEN 
                       RASAL(IC) = DT / (RASAL1 + (RASAL2-RASAL1)*(ZET(IC) - 2000.)/8000.)
                       !RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
                       !RASAL(IC) = DT / RASAL(IC)
                   ENDIF
      
                   CVW(IC:K) = MAX(  CVW(IC:K) , 1.00 )
                         
                   !TEST FOR CRITICAL WORK FUNCTION
                   CALL ACRITN(POL(IC), PRS(K), ACRITFAC, ACR)
      
                   IF (WFN <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION ======>>
                      RC = 5
                      !RETURN
                      !====EXIT====>
                   ENDIF
      
                   !CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
                   IF (RC == 10) THEN

                      WLQ     = QOL(K)
                      UHT     = UOI(K)
                      VHT     = VOI(K)
                      RNN(K)  = 0.
                      CLL0(K) = 0.
     
                      DO L=K-1,IC,-1
                         TEM6   = ETA(L) - ETA(L+1)
                         WLQ   = WLQ + TEM6 * QOL(L)
                         UHT   = UHT + TEM6 * UOI(L)
                         VHT   = VHT + TEM6 * VOI(L)
         
                         !!!! How much condensate (CLI) is present here? 
                         IF (L>IC) THEN
                            TX2   = 0.5*(QST(L)+QST(L-1))*ETA(L)
                            TX3   = 0.5*(HST(L)+HST(L-1))*ETA(L)
                            QCC   = TX2 + GM1(L)*(HCC(L)-TX3)
                            CLL0(L)  = (WLQ-QCC)
                         ELSE
                            CLL0(L)   = (WLQ-QST(IC)*ETA(IC))
                         ENDIF

                         CLL0(L)    = MAX( CLL0(L) , 0.00 )
    
                         
                         TE_A = POI(L)*PRH(L)     ! Temperature (K)
     
                         if (abs(ETA(L)) .gt. 0.0) then !Linearisation security

                            !SUNDQ3_ICE
                            F2   = 1.0
                            IF ( ( TE_A .GE. TE1 ) .AND. ( TE_A .LT. TE0 ) )THEN 
                               F2   = 1.0 + JUMP1 * (( TE0 - TE_A )**0.3333)
                            ELSEIF   ( TE_A .LT. TE1 ) THEN 
                               F2   = SDQV2 + (SDQV3-SDQV2)*(TE1-TE_A)/(TE1-TE2)
                            ENDIF
                            !end sundq3_ice
                                        
                            C00_X  =  CCO_AUTO(I) * F2 * F4  ! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                            CLI_CRIT_X =  CLI_CRIT / F2 
                            CLI  = CLL0(L) / ETA(L)  ! condensate (kg/kg)  
                            RATE = C00_X * ( 1.0 - EXP( -(CLI)**2 / (CLI_CRIT_X**2+10e-15) ) )

                            CVW_X = MAX( CVW(L) , 1.00 )  ! Floor conv. vel. since we don't real(8)ly trust it at low values
                   
                            DT_LYR  = ( ZET(L)-ZET(L+1) )/CVW_X           ! l.h.s. DT_LYR => time in layer (L,L+1)

                            CLOSS   = CLL0(L) * RATE * DT_LYR
                            CLOSS   = MIN( CLOSS , CLL0(L) )

                            CLL0(L) = CLL0(L) - CLOSS
                            DLL0(L) = CLOSS

                            IF (CLOSS>0.) then
                               WLQ = WLQ - CLOSS
                               RNN(L) = CLOSS 
                            else
                               RNN(L) = 0.
                            ENDIF

                         endif 

                      ENDDO

                      WLQ = WLQ - QST(IC)*ETA(IC)

                      !CALCULATE GAMMAS AND KERNEL
                      GMS(K) =          (SHT(K)-SSL(K))*PRI(K)          ! MS-A30 (W/O GRAV)
                      GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*ALHL     ! MS-A31 (W/O GRAV)
                      AKM    = GMH(K)*GAM(K-1)*DPB(K-1)                 ! MS-A37 (W/O GRAV)

                      TX2     = GMH(K)
                      DO L=K-1,IC+1,-1
                         GMS(L) =          ( ETA(L  )*(SHT(L)-SSL(L  )) + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
                         GMH(L) = GMS(L) + ( ETA(L  )*(QHT(L)-QOL(L  )) + ETA(L+1)*(QOL(L)-QHT(L+1)) )*ALHL*PRI(L)
                         TX2 = TX2 + (ETA(L) - ETA(L+1)) * GMH(L)
                         AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2*GHT(L)
                      ENDDO

                      GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
                      AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)

                      GMH(IC) = GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL + ETA(IC  )*(HST(IC)-HOL(IC  ))     )*PRI(IC)
     
                      !CLOUD BASE MASS FLUX
      
                      IF (AKM >= 0.0 .OR. WLQ < 0.0)  THEN  !  =========>
                         RC = 6
                         !RETURN
                         !====EXIT====>
                      ENDIF
      
                      IF (RC == 10) then
      
                         if (abs(AKM) .gt. 0.0) then !Linearisation security
                            WFN1 = - (WFN-ACR)/AKM ! MS-A39 MASS-FLUX IN Pa/step
                            WFN1 = MIN( ( RASAL(IC)*TRG*TOKI )*WFN1  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))
                         endif
     
                         !CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                         TEM7      = WFN1*GRAVI
                         CLL (IC) = CLL (IC) + WLQ*TEM7           ! (kg/m^2/step)
                         RMF (IC) = RMF (IC) +     TEM7           ! (kg/m^2/step)
                         RMFD(IC) = RMFD(IC) +     TEM7 * ETA(IC) ! (kg/m^2/step)

                         DO L = IC+1,K
                            RMFP(L) = TEM7 * ETA(L)                 ! (kg/m^2/step)
                            RMFC(L) = RMFC(L)  +  RMFP(L)          ! (kg/m^2/step)
           
                            DLLX(L) = DLLX(L)  +  TEM7*DLL0(L)        

                            IF ( CVW(L) > 0.0 ) THEN
                               if ( abs(CVW(L) * PRS(L)) .gt. 0.0) then !Linearisation security
                                  UPDFRP(L) = rmfp(L) * (DDT/DAYLEN) * 1000. / ( CVW(L) * PRS(L) )
                               endif
                            ELSE 
                               UPDFRP(L) = 0.0       
                            ENDIF

                            if (abs(ETA(L)) .gt. 0.0) then !Linearisation security
                               CLLI(L) = CLL0(L)/ETA(L)  ! current cloud; incloud condensate        
                            endif
                            CLLB(L) = CLLB(L) +  UPDFRP(L) * CLLI(L) !  cumulative grid mean convective condensate        

                            UPDFRC(L) =  UPDFRC(L) +  UPDFRP(L)      
     
                         ENDDO

                         !THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                         DO L = IC,K
                            RNS(L) = RNS(L) + RNN(L)*TEM7 ! (kg/m^2/step)
                            GMH(L) = GMH(L) * WFN1
                            GMS1(L) = GMS(L) * WFN1
                            QOI(L) = QOI(L) + (GMH(L) - GMS1(L)) * ALHI
                            POI(L) = POI(L) + GMS1(L)*PKI(L)*CPI
                            QST(L) = QST(L) + GMS1(L)*BET(L)*CPI
                         ENDDO
      
                         !CUMULUS FRICTION
                         IF(FRICFAC <= 0.0) THEN
                            RC = 6
                            !RETURN  !  NO CUMULUS FRICTION =========>>
                         ENDIF

                         IF (RC == 10) THEN

                            WFN2 = WFN1*0.5 *1.0           !*FRICFAC*0.5

                            WFN2     = WFN2*FRICFAC*EXP( -ALM / FRICLAMBDA )
                            TEM8     = WFN2*PRI(K)

                            UCU(K)  = UCU(K) + TEM8 * (UOI(K-1) - UOI(K))
                            VCU(K)  = VCU(K) + TEM8 * (VOI(K-1) - VOI(K))

                            DO L = K-1,IC+1,-1
                               TEM9    = WFN2*PRI(L)
                               UCU(L) = UCU(L) + TEM9 * ( (UOI(L-1) - UOI(L)) * ETA(L) + (UOI(L) - UOI(L+1)) * ETA(L+1) )
                               VCU(L) = VCU(L) + TEM9 * ( (VOI(L-1) - VOI(L)) * ETA(L) + (VOI(L) - VOI(L+1)) * ETA(L+1) )
                            ENDDO

                            TEM10     = WFN2*PRI(IC)
                            UCU(IC) = UCU(IC) + (2.*(UHT - UOI(IC)*(ETA(IC)-ETA(IC+1))) - (UOI(IC)+UOI(IC+1))*ETA(IC+1))*TEM10
                            VCU(IC) = VCU(IC) + (2.*(VHT - VOI(IC)*(ETA(IC)-ETA(IC+1))) - (VOI(IC)+VOI(IC+1))*ETA(IC+1))*TEM10
    
                            DO L = IC,K
                               UOI(L) = UOI(L) + UCU(L)
                               VOI(L) = VOI(L) + VCU(L)
                            ENDDO

                         ENDIF

                         RC = 0 !DID EVERYTHING (Exc Cumulus Friction)

                      endif !RC = 6
                   endif !RC = 5
                endif !RC = 4
             endif !RC = 3
          endif !RC = 2
       endif !RC = 1

    ENDIF

 ENDDO

 RETURN

END SUBROUTINE CLOUDE


SUBROUTINE STRAP0(I,K,K0,IDIM,ICMIN,CPKE,CPLE,UHO,VHO,THO,QHO,CDQS,CQSS,ESTBLX,SEEDRAS,WGT0,MXDIAM,&
                  POI,QOI,UOI,VOI,PRJ,PRS,QST,DQQ,CONS_RGAS,CONS_CP,ONEPKAP,&
                  POL,PRH,PKI,DPT,DPB,PRI,MAXDALLOWED,BET,GAM,LBCP,GHT,GM1,ALHL,TCU,QCU,&
                  RNS,CLL,RMF,RMFD,RMFC,RMFP,CLL0,DLL0,DLLX,CLLI,CLLB,CVW,UPDFRC,UPDFRP,&
                  POI_SV,QOI_SV,UOI_SV,VOI_SV)

 IMPLICIT NONE

 !!GLOBALS!!
 INTEGER :: I, K, K0, IDIM, ICMIN

 INTEGER, DIMENSION(IDIM) :: SEEDRAS
 real(8), DIMENSION(IDIM,K0) :: UHO, VHO, THO, QHO, CDQS, CQSS, WGT0
 real(8), DIMENSION(IDIM,K0+1) :: CPKE, CPLE

 real(8) :: CONS_RGAS,CONS_CP
 real(8) :: ONEPKAP,MAXDALLOWED, LBCP, ALHL, MXDIAM
 real(8), DIMENSION(K0+1) :: PRJ, PRS
 real(8), DIMENSION(K0) :: POI, QOI, UOI, VOI, QST, DQQ, TEMP1
 real(8), DIMENSION(K0) :: POL,PRH,PKI,DPT,DPB,PRI
 real(8), DIMENSION(K0) :: BET, GAM, GHT, GM1, TCU, QCU
 real(8), DIMENSION(K0) :: POI_SV,QOI_SV,UOI_SV,VOI_SV
 real(8), DIMENSION(K0) :: RNS,CLL,RMF,RMFD,RMFC,RMFP,CLL0,DLL0,DLLX,CLLI,CLLB,CVW,UPDFRC,UPDFRP

 !!LOCALS!!
 INTEGER :: L, KK

 real(8) :: WGHT0, PRCBL, rndu
 real(8), DIMENSION(K0) :: WGHT, MASSF

 !Table lookup constants
 real(8) :: ESTBLX(:)

 !Zero locals
 WGHT0 = 0.0
 PRCBL = 0.0
 rndu = 0.0
 WGHT = 0.0
 MASSF = 0.0

 do kk=icmin,k+1
    PRJ(kk) = CPKE(I,kk)
 enddo

 poi=0.        ! These initialized here in order not to confuse Valgrind debugger
 qoi=0.        ! Do not believe it actually makes any difference.
 uoi=0.
 voi=0.     

 PRS(ICMIN:K0+1) = CPLE(I,ICMIN:K0+1)
 POI(ICMIN:K)   = THO(I,ICMIN:K)
 QOI(ICMIN:K)   = QHO(I,ICMIN:K)
 UOI(ICMIN:K)   = UHO(I,ICMIN:K)
 VOI(ICMIN:K)   = VHO(I,ICMIN:K)

 QST(ICMIN:K) = CQSS(I,ICMIN:K)
 DQQ(ICMIN:K) = CDQS(I,ICMIN:K)
       
 !!! Mass fraction of each layer below cloud base
 !!! contributed to aggregate cloudbase layer (CBL) 
 MASSF(:) = WGT0(I,:)

 !!! RESET PRESSURE at bottom edge of CBL 
 PRCBL = PRS(K)
 do l = K,K0
    PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
 end do
 PRS(K+1) = PRCBL
 PRJ(K+1) = (PRS(K+1)/1000.)**(CONS_RGAS/CONS_CP)

 DO L=K,ICMIN,-1
    POL(L)  = 0.5*(PRS(L)+PRS(L+1))
    PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) &
              / (ONEPKAP*(PRS(L+1)-PRS(L)))
    PKI(L)  = 1.0 / PRH(L)
    DPT(L)  = PRH(L  ) - PRJ(L)
    DPB(L)  = PRJ(L+1) - PRH(L)
    PRI(L)  = .01 / (PRS(L+1)-PRS(L))
 ENDDO

 !!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
 if ( K <= K0 ) then
    POI(K) = 0.
    QOI(K) = 0.
    UOI(K) = 0.
    VOI(K) = 0.

    !! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
    WGHT = 0.
    DO L = K,K0
       WGHT(L)   = MASSF(L) * ( CPLE(I,L+1)-CPLE(I,L) ) / ( PRS(K+1)-PRS(K) )
    END DO      

    DO L = K,K0
       POI(K) = POI(K) + WGHT(L)*THO(I,L)
       QOI(K) = QOI(K) + WGHT(L)*QHO(I,L)
       UOI(K) = UOI(K) + WGHT(L)*UHO(I,L)
       VOI(K) = VOI(K) + WGHT(L)*VHO(I,L)
    ENDDO

    DQQ(K) = 0.0
    QST(K) = 0.0

    TEMP1(K) = POI(K)*PRH(K)

    call DQSAT_sub(DQQ(K),QST(K),TEMP1(K),POL(K),1,1,1,ESTBLX)

 endif

 rndu = max( seedras(I)/1000000., 1e-6 )

 MXDIAM = maxdallowed*( rndu**(-1./2.) )
 DO L = K,ICMIN,-1
    BET(L)  = DQQ(L)*PKI(L)  !*
    GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L)) !*
    IF (L < K) THEN
       GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
       GM1(L+1) = 0.5*LBCP*(DQQ(L  )/(ALHL*(1.0+LBCP*DQQ(L  ))) + DQQ(L+1)/(ALHL*(1.0+LBCP*DQQ(L+1))) )
    ENDIF
 ENDDO

 TCU(ICMIN:K) = -POI(ICMIN:K)*PRH(ICMIN:K)
 QCU(ICMIN:K) = -QOI(ICMIN:K)

 RNS  = 0.
 CLL  = 0.
 RMF  = 0.
 RMFD = 0.
 RMFC = 0.
 RMFP = 0.
 CLL0 = 0.
 DLL0 = 0.
 DLLX = 0.
 CLLI = 0.
 CLLB = 0.

 POI_SV = POI
 QOI_SV = QOI
 UOI_SV = UOI
 VOI_SV = VOI

 CVW     = 0.0
 UPDFRC  = 0.0
 UPDFRP  = 0.0

end SUBROUTINE STRAP0

SUBROUTINE STRAP1(I,K,K0,IDIM,ICMIN,&
                  THO,QHO,UHO,VHO,CNV_UPDFRC,CPLE,FLXD,CLW,WGT1,&
                  POI,QOI,UOI,VOI,UPDFRC,CVW,CLLB,PRS,POI_SV,QOI_SV,UOI_SV,VOI_SV,&
                  RMF,RMFD,RMFC,CLL,DDT,DAYLEN )

 IMPLICIT NONE

 !!GLOBALS!!
 INTEGER :: I, K, K0, IDIM, ICMIN
 real(8) :: DDT, DAYLEN

 real(8), DIMENSION(IDIM,K0) :: THO, QHO, UHO, VHO
 real(8), DIMENSION(IDIM,K0) :: CNV_UPDFRC, FLXD, CLW, WGT1
 real(8), DIMENSION(IDIM,K0+1) :: CPLE

 real(8), DIMENSION(K0) :: POI, QOI, UOI, VOI, POI_SV, QOI_SV, UOI_SV, VOI_SV
 real(8), DIMENSION(K0) :: UPDFRC, CVW, CLLB, PRS, RMF, RMFD, RMFC, CLL
 
 !!LOCALS!! 
 INTEGER :: L
 real(8) :: WGHT0
 real(8), DIMENSION(K0) :: WGHT
 
 !Zero locals
 WGHT0 = 0.0
 WGHT = 0.0
 
 THO(I,ICMIN:K-1) = POI(ICMIN:K-1)
 QHO(I,ICMIN:K-1) = QOI(ICMIN:K-1)
 UHO(I,ICMIN:K-1) = UOI(ICMIN:K-1)
 VHO(I,ICMIN:K-1) = VOI(ICMIN:K-1)
 CNV_UPDFRC(I,ICMIN:K-1)   =  UPDFRC(ICMIN:K-1)

 WGHT   = WGT1(I,:)

 !! Scale properly by layer masses
 wght0 = 0.
 DO L = K,K0 
    wght0 = wght0 + WGHT(L)* ( CPLE(I,L+1) - CPLE(I,L) )
 ENDDO
         
 wght0 = ( PRS(K+1)   - PRS(K)  )/wght0

 WGHT  = wght0 * WGHT


 DO L = K,K0 
    THO(I,L) =  THO(I,L) + WGHT(L)*(POI(K) - POI_SV(K))
    QHO(I,L) =  QHO(I,L) + WGHT(L)*(QOI(K) - QOI_SV(K))
    UHO(I,L) =  UHO(I,L) + WGHT(L)*(UOI(K) - UOI_SV(K))
    VHO(I,L) =  VHO(I,L) + WGHT(L)*(VOI(K) - VOI_SV(K))
 ENDDO

 FLXD(I,ICMIN:K) = RMFD(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
 CLW (I,ICMIN:K) = CLL (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s )

 FLXD(I,1:ICMIN-1) = 0.
 CLW (I,1:ICMIN-1) = 0.

 IF ( K < K0 ) THEN 
    FLXD(I,K:K0) = 0.
    CLW (I,K:K0) = 0.
 END IF

  
END SUBROUTINE STRAP1



SUBROUTINE ACRITN( PL, PLB, ACRITFAC, ACR)
      
      IMPLICIT NONE

      real(8), INTENT(IN ) :: PL, PLB, ACRITFAC
      real(8), INTENT(OUT) :: ACR

      integer :: IWK
      
      !!real(8), PARAMETER :: FACM=0.5

      real(8) :: PH(15), A(15)

      PH(1) = 150.0
      PH(2) = 200.0
      PH(3) = 250.0
      PH(4) = 300.0
      PH(5) = 350.0
      PH(6) = 400.0
      PH(7) = 450.0
      PH(8) = 500.0
      PH(9) = 550.0
      PH(10) = 600.0
      PH(11) = 650.0
      PH(12) = 700.0
      PH(13) = 750.0
      PH(14) = 800.0
      PH(15) = 850.0

      A(1)  = 1.6851
      A(2)  = 1.1686
      A(3)  = 0.7663
      A(4)  = 0.5255
      A(5)  = 0.4100
      A(6)  = 0.3677
      A(7)  = 0.3151
      A(8)  = 0.2216
      A(9)  = 0.1521
      A(10)  = 0.1082
      A(11)  = 0.0750
      A(12)  = 0.0664
      A(13)  = 0.0553
      A(14)  = 0.0445
      A(15)  = 0.0633

      IWK = INT(PL * 0.02 - 0.999999999)

      IF (IWK .GT. 1 .AND. IWK .LE. 15) THEN
       ACR = A(IWK-1) + (PL-PH(IWK-1))*.02*(A(IWK)-A(IWK-1))
      ELSEIF(IWK > 15) THEN
       ACR = A(15)
      ELSE
       ACR = A(1)
      ENDIF

      ACR = ACRITFAC  * ACR * (PLB - PL)

      RETURN

   END SUBROUTINE ACRITN



subroutine DQSAT_sub(DQSi,QSSi,TEMP,PLO,im,jm,lm,ESTBLX)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 integer :: im,jm,lm
 real(8), dimension(im,jm,lm) :: TEMP, PLO
 real(8) :: ESTBLX(:)

 !Outputs
 real(8), dimension(im,jm,lm) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8), parameter :: CONS_H2OMW = 18.01, CONS_AIRMW = 28.97
 real(8),    parameter :: ESFAC = CONS_H2OMW/CONS_AIRMW

 integer :: i, j, k

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 integer, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 do K=1,LM
    do J=1,JM
       do I=1,IM

          TL = TEMP(I,J,K)
          PL = PLO(I,J,K)

          PP = PL*100.0

          if (TL<=TMINTBL) then
             TI = TMINTBL
          elseif(TL>=TMAXTBL-.001) then
             TI = TMAXTBL-.001
          else
             TI = TL
          end if

          TT = (TI - TMINTBL)*DEGSUBS+1
          IT = int(TT)

          DQQ =  ESTBLX(IT+1) - ESTBLX(IT)
          QQ  =  (TT-IT)*DQQ + ESTBLX(IT)

          if (PP <= QQ) then
             QSAT = MAX_MIXING_RATIO
             DQSAT = 0.0
          else
             DD = 1.0/(PP - (1.0-ESFAC)*QQ)
             QSAT = ESFAC*QQ*DD
             DQSAT = (ESFAC*DEGSUBS)*DQQ*PP*(DD*DD)
          end if

          DQSi(I,J,K) = DQSAT
          QSSi(I,J,K) = QSAT

       end do
    end do
 end do

end subroutine DQSAT_sub
