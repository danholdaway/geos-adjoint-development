module CONVECTION

IMPLICIT NONE

PRIVATE
PUBLIC :: RASE, ACRITN, SUNDQ3_ICE, DQSATs_RAS, DQSAT_RAS

CONTAINS

 SUBROUTINE RASE( K0             , &
                  ICMIN          , &
                  DT             , &
                  SEEDRAS        , &
                  SIGE           , &
                  KCBL           , &
                  WGT0           , &
                  WGT1           , &
                  FRLAND         , &
                  TS             , &
                  THO            , &
                  QHO            , &
                  UHO            , &
                  VHO            , &
                  CO_AUTO        , &
                  PLE            , &
                  PKE            , &
                  CLW            , & !CNV_DQLDT
                  FLXD           , & !CNV_MFD
                  CNV_PRC3       , & !CNV_PRC3
                  CNV_UPDFRC     , & !CNV_UPDFRC
                  RASPARAMS      , &
                  CONS_CP        , &
                  CONS_ALHL      , &
                  CONS_ALHS      , &
                  CONS_GRAV      , &
                  CONS_UNDEF     , &
                  CONS_RGAS      , &
                  CONS_H2OMW     , &
                  CONS_AIRMW     , &
                  CONS_VIREPS    , &
                  ESTBLX         , &
                  TAP_DUMMY      , &
                  ITRCR          , &
                  XHO            , &
                  FSCAV          , &
                  DISSKE           )

 IMPLICIT NONE

 !INPUTS
 INTEGER, INTENT(IN) :: K0, ICMIN, KCBL, SEEDRAS
 REAL(8),    INTENT(IN) :: DT, FRLAND, TS, CO_AUTO, TAP_DUMMY

 REAL(8), DIMENSION(K0  ), INTENT(IN) :: WGT0, WGT1
 REAL(8), DIMENSION(K0+1), INTENT(IN) :: PLE, PKE, SIGE
 REAL(8), DIMENSION(:)   , INTENT(IN) :: RASPARAMS, ESTBLX

 REAL(4), INTENT(IN) :: CONS_ALHS, CONS_CP, CONS_ALHL
 REAL(4), INTENT(IN) :: CONS_GRAV, CONS_UNDEF, CONS_RGAS
 REAL(4), INTENT(IN) :: CONS_H2OMW, CONS_AIRMW, CONS_VIREPS

 !OUTPUTS
! REAL(8), DIMENSION (K0+1), INTENT(OUT) ::  FLXC
! REAL(8), DIMENSION (K0  ), INTENT(OUT) ::  FLX, CNV_CVW, CNV_QC, ENTLAM
! REAL(8), DIMENSION (K0+1)              ::  FLXC
! REAL(8), DIMENSION (K0  )              ::  FLX, CNV_CVW, CNV_QC, ENTLAM

 !PROGNOSTIC
 REAL(8), DIMENSION (K0  ), INTENT(INOUT) ::  THO, QHO, UHO, VHO
 REAL(8), DIMENSION (K0  ), INTENT(INOUT) ::  CLW, FLXD, CNV_PRC3, CNV_UPDFRC

 REAL(8), DIMENSION (K0  ) ::  CLWTMP, FLXDTMP, CNV_PRC3TMP, CNV_UPDFRCTMP, DISSKETMP

 !SCAVENGING
 INTEGER,              INTENT(IN)    :: ITRCR
 REAL(8)   , OPTIONAL, INTENT(IN)    :: FSCAV(ITRCR) ! = 0 (no scav), = 1 (full scav)
 REAL(8)   , OPTIONAL, INTENT(INOUT) :: DISSKE(K0)        !OUT
 REAL(8)   , OPTIONAL, INTENT(INOUT) :: XHO(K0,ITRCR)     !INOUT

! REAL(8)   , OPTIONAL, INTENT(INOUT) :: TRIEDLEV_DIAG(K0) !OUT

 !LOCALS
 REAL(8), PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
 REAL(8), PARAMETER :: RHMAX   = 0.9999
 REAL(8), PARAMETER :: CBL_QPERT = 0.0, CBL_TPERT = 1.0
 REAL(8), PARAMETER :: CBL_TPERT_MXOCN = 2.0, CBL_TPERT_MXLND = 4.0

 INTEGER :: K, KK, IC, L, ICL , ITR , ICL_C, N_DTL, RC(K0-1)
 INTEGER, ALLOCATABLE :: ICL_V(:)

 !Constants & parameters
 REAL(8) :: CPBG, ALHI, CPI, GRAVI, DDT, LBCP
 REAL(8) :: FRICFAC, PBLFRAC, AUTORAMPB, CO_ZDEP, RASAL1, RASAL2
 REAL(8) :: FRICLAMBDA, SDQVT1, SDQV2, SDQV3, ACRITFAC, CLI_CRIT
 REAL(8) :: MAXDALLOWED, RHMN, RHMX

 REAL(8) :: LAMBDA_MIN, LAMBDA_MAX
 REAL(8) :: TPERT, QPERT, MXDIAM, SEEDRAS1
 REAL(8) :: TX2, TX3, UHT, VHT, AKM, ACR, ALM, DQX
 REAL(8) :: WFN, TEM, TRG, WLQ, QCC
 REAL(8) :: WFNOG
 REAL(8) :: PRCBL, rndu 
 REAL(8) :: CLI , TE_A, C00_X, CLI_CRIT_X,  TOKI, GMHx, HSTx
 REAL(8) :: DT_LYR, RATE, CVW_X, CLOSS, F2, F3, F4, F5

 REAL(8), DIMENSION(K0+1) :: PRJ,PRS,QHT,SHT,ZET
 REAL(8), DIMENSION(0:K0) :: ZLE
 REAL(8), DIMENSION(K0)   :: QSS,DQS,TEMPf,Pf,PK,ZLO
 REAL(8), DIMENSION(K0)   :: POI_SV,QOI_SV,UOI_SV,VOI_SV
 REAL(8), DIMENSION(K0)   :: POI,QOI,UOI,VOI,DQQ,BET,GAM,CLL
 REAL(8), DIMENSION(K0)   :: POI_c,QOI_c
 REAL(8), DIMENSION(K0)   :: PRH, PRI, GHT, DPT, DPB, PKI
 REAL(8), DIMENSION(K0)   :: UCU, VCU, RNS, POL, DISSK0
 REAL(8), DIMENSION(K0)   :: QST, SSL, RMF, RNN, RMFP, RMFD
 REAL(8), DIMENSION(K0)   :: GMS, ETA, GMH, EHT, GM1, HCC
 REAL(8), DIMENSION(K0)   :: HOL, HST, QOL, ZOL, HCLD, CLL0
 REAL(8), DIMENSION(K0)   :: BKE , CVW, CVW1, UPDFRC
 REAL(8), DIMENSION(K0)   :: RASAL, UPDFRP,BK2,BK3
 REAL(8), DIMENSION(K0)  :: WGHT, WGHT0, WGHT1, MASSF
! REAL(8), DIMENSION(K0)   :: RMFC, CLLB, CLLI

 !SCAVANGING RELATED PARAMETERS
 LOGICAL                       :: DO_TRACERS, SMOOTH_HST
 REAL(8)                       :: DELZKM  ! layer thickness in km
 REAL(8)                       :: FNOSCAV ! fraction of tracer *not* scavenged
 REAL(8)                       :: FSCAV_(ITRCR) ! Fraction scavenged per km
 REAL(8),  DIMENSION(ITRCR)    :: XHT
 REAL(8),  DIMENSION(K0,ITRCR) :: XOI, XCU, XOI_SV


! BEGIN CALCULATIONS
!-------------------

 !Initialize local arrays
 PRJ = 0.0
 PRS = 0.0
 QHT = 0.0
 SHT = 0.0
 ZET = 0.0
 ZLE = 0.0
 QSS = 0.0
 DQS = 0.0
 TEMPf = 0.0
 Pf = 0.0
 PK = 0.0
 ZLO = 0.0
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
 UCU = 0.0
 VCU = 0.0
 RNS = 0.0
 POL = 0.0
 DISSK0 = 0.0
 QST = 0.0
 SSL = 0.0
 RMF = 0.0
 RNN = 0.0
 RMFP = 0.0
 RMFD = 0.0
 GMS = 0.0
 ETA = 0.0
 GMH = 0.0
 EHT = 0.0
 GM1 = 0.0
 HCC = 0.0
 HOL = 0.0
 HST = 0.0
 QOL = 0.0
 ZOL = 0.0
 HCLD = 0.0
 CLL0 = 0.0
 BKE = 0.0
 CVW = 0.0
 CVW1 = 0.0
 UPDFRC = 0.0
 RASAL = 0.0
 UPDFRP = 0.0
 BK2 = 0.0
 BK3 = 0.0
 WGHT = 0.0
 WGHT0 = 0.0
 WGHT1 = 0.0
 MASSF = 0.0

 IF ( PRESENT(FSCAV) ) then
    FSCAV_ = FSCAV
 ELSE
    FSCAV_ = 0.0 ! NO SCAVENGING BY DEFAULT
 ENDIF      

! IF (PRESENT(TRIEDLEV_DIAG)) THEN
!    TRIEDLEV_DIAG=0.
! ENDIF

 !Make outputs prognostic to trick the autodiff
 CLWTMP        = CLW        * TAP_DUMMY
 FLXDTMP       = FLXD       * TAP_DUMMY
 CNV_PRC3TMP   = CNV_PRC3   * TAP_DUMMY
 CNV_UPDFRCTMP = CNV_UPDFRC * TAP_DUMMY
 if ( PRESENT( DISSKE )) THEN
    DISSKETMP     = DISSKE     * TAP_DUMMY
 endif

! CNV_CVW    = 0. 
! CNV_QC     = 0.
! ENTLAM     = 0.

 if ( PRESENT( DISSKE )) THEN
    DISSKE    =0.
 endif
      
 SMOOTH_HST   = .FALSE.  

 FRICFAC      = RASPARAMS(1)     !  ---  1
 CLI_CRIT     = RASPARAMS(4)     !  ---  4
 RASAL1       = RASPARAMS(5)     !  ---  5
 RASAL2       = RASPARAMS(6)     !  ---  6
 FRICLAMBDA   = RASPARAMS(11)    !  --- 11
 SDQV2        = RASPARAMS(14)    !  --- 14
 SDQV3        = RASPARAMS(15)    !  --- 15
 SDQVT1       = RASPARAMS(16)    !  --- 16
 ACRITFAC     = RASPARAMS(17)    !  --- 17
 PBLFRAC      = RASPARAMS(20)    !  --- 20
 AUTORAMPB    = RASPARAMS(21)    !  --- 21
 CO_ZDEP      = RASPARAMS(22)    !  --- 22
 MAXDALLOWED  = RASPARAMS(23)    !  --- 23
 RHMN         = RASPARAMS(24)    !  --- 24
 RHMX         = RASPARAMS(25)    !  --- 25

 DO_TRACERS = present(XHO) .and. ITRCR > 0

 CPI   = 1.0/CONS_CP      
 ALHI  = 1.0/CONS_ALHL
 GRAVI = 1.0/CONS_GRAV
 CPBG  = CONS_CP*GRAVI
 DDT   = DAYLEN/DT
 LBCP  = CONS_ALHL*CPI

 !Get saturation specific humidity and gradient wrt to T
 Pf    = 0.5*(PLE(1:K0) +  PLE(2:K0+1  ) )
 PK   = (Pf/1000.)**(CONS_RGAS/CONS_CP)
 TEMPf = THO(:)*PK

 ZLE(K0+1) = 0.
 do L=K0,1,-1
    ZLE(L) = THO(L) * (1.+CONS_VIREPS*QHO(L))
    ZLO(L) = ZLE(L+1) + (CONS_CP/CONS_GRAV)*( PKE(L+1)-PK (L  ) ) * ZLE(L)
    ZLE(L) = ZLO(L)   + (CONS_CP/CONS_GRAV)*( PK (L)    -PKE(L) ) * ZLE(L)
 end do

 TPERT  = CBL_TPERT * ( TS - ( TEMPf(K0)+ CONS_GRAV*ZLO(K0)/CONS_CP )  ) 
 QPERT  = CBL_QPERT !* ( QSSFC - Q(:,:,K0) ) [CBL_QPERT = 0.0]
 TPERT  = MAX( TPERT , 0.0 )
 QPERT  = MAX( QPERT , 0.0 )

 if (FRLAND < 0.1) then
    TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
 else
    TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
 endif

 call DQSAT_RAS(DQS,QSS,TEMPf,Pf,1,1,K0,ESTBLX,CONS_H2OMW, CONS_AIRMW)

 !CALL FINDBASE
 K = KCBL

 rc(icmin) = 0

 !FIND DTLS
 N_DTL = K - ICMIN 

 if(allocated(ICL_V)) deallocate(ICL_V)
 allocate(ICL_V(N_DTL))

 do L=1,N_DTL
    ICL_V(L) = K - L
 enddo
 !---------


 !STRAP FINAL=0 SUBROUTINE
 !------------------
 do kk=icmin,k+1
    PRJ(kk) = PKE(kk)
 enddo

 poi=0.        ! These initialized here in order not to confuse Valgrind debugger
 qoi=0.        ! Do not believe it actually makes any difference.
 uoi=0.
 voi=0.     

 PRS(ICMIN:K0+1) = PLE(ICMIN:K0+1)
 POI(ICMIN:K)   = THO(ICMIN:K)
 QOI(ICMIN:K)   = QHO(ICMIN:K)
 UOI(ICMIN:K)   = UHO(ICMIN:K)
 VOI(ICMIN:K)   = VHO(ICMIN:K)

 QST(ICMIN:K) = QSS(ICMIN:K)
 DQQ(ICMIN:K) = DQS(ICMIN:K)

 IF (DO_TRACERS) THEN 
     DO ITR=1,ITRCR 
         XOI(ICMIN:K,ITR) = XHO(ICMIN:K,ITR)
     END DO
 END IF

 !Mass fraction of each layer below cloud base contributed to aggregate cloudbase layer 
 MASSF(:) = WGT0(:)

 !RESET PRESSURE at bottom edge of CBL 
 PRCBL = PRS(K)
 do l= K,K0
     PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
 end do
 PRS(K+1) = PRCBL
 PRJ(K+1) = (PRS(K+1)/1000.)**(CONS_RGAS/CONS_CP)

 DO L=K,ICMIN,-1
    POL(L)  = 0.5*(PRS(L)+PRS(L+1))
    PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) / (ONEPKAP*(PRS(L+1)-PRS(L)))
    PKI(L)  = 1.0 / PRH(L)
    DPT(L)  = PRH(L  ) - PRJ(L)
    DPB(L)  = PRJ(L+1) - PRH(L)
    PRI(L)  = .01 / (PRS(L+1)-PRS(L))
 ENDDO

 !RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
 if ( K <= K0) then

    POI(K) = 0.
    QOI(K) = 0.
    UOI(K) = 0.
    VOI(K) = 0.

    !SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
    WGHT = 0.
    DO L=K,K0
       WGHT(L) = MASSF(L) * ( PLE(L+1) - PLE(L) ) / ( PRS(K+1) - PRS(K) )
    END DO      

    DO L=K,K0
       POI(K) = POI(K) + WGHT(L)*THO(L)
       QOI(K) = QOI(K) + WGHT(L)*QHO(L)
       UOI(K) = UOI(K) + WGHT(L)*UHO(L)
       VOI(K) = VOI(K) + WGHT(L)*VHO(L)
    ENDDO

    IF (DO_TRACERS) THEN 
       XOI(K,:)=0.
       DO ITR=1,ITRCR
          DO L=K,K0
             XOI(K,ITR) = XOI(K,ITR) + WGHT(L)*XHO(L,ITR)
          END DO
       END DO
    END IF

    call DQSATs_RAS(DQQ(K), QST(K), POI(K)*PRH(K),POL(K),ESTBLX,CONS_H2OMW, CONS_AIRMW)

 endif
 
 seedras1 = max( seedras/1000000., 1e-6 )
 rndu = seedras1
 MXDIAM = maxdallowed*( rndu**(-1./2.) )

 DO L=K,ICMIN,-1
    BET(L)  = DQQ(L)*PKI(L)  !*
    GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L)) !*
    IF (L < K) THEN
       GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
       GM1(L+1) = 0.5*LBCP*(DQQ(L  )/(CONS_ALHL*(1.0+LBCP*DQQ(L  ))) + &
                            DQQ(L+1)/(CONS_ALHL*(1.0+LBCP*DQQ(L+1))) )
    ENDIF
 ENDDO

 POI_SV = POI
 QOI_SV = QOI
 UOI_SV = UOI
 VOI_SV = VOI

 IF (DO_TRACERS) THEN
    XOI_SV = XOI
 ENDIF

 CVW     = 0.0
 UPDFRC  = 0.0
 UPDFRP  = 0.0
 DISSK0  = 0.0
 !END STRAP (FINAL=0) SUBROUTINE
 !------------------------------

 !HTEST 
 !-----
 HOL = 0.0
 ZET(K+1) = 0
 SHT(K+1) = CONS_CP*POI(K)*PRJ(K+1)
 DO L=K,ICMIN,-1
    QOL(L)  = MIN(QST(L)*RHMAX,QOI(L))
    QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
    SSL(L)  = CONS_CP*PRJ(L+1)*POI(L) + CONS_GRAV*ZET(L+1)
    HOL(L)  = SSL(L) + QOL(L)*CONS_ALHL
    HST(L)  = SSL(L) + QST(L)*CONS_ALHL
    TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
    ZET(L)  = ZET(L+1) + TEM
    ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
 ENDDO
 !HTEST

 DO ICL_C = 1,N_DTL

    ICL   = ICL_V( ICL_C )

    IF (DO_TRACERS) THEN
       XCU(ICMIN:,:) = 0.
    ENDIF

!    IF (PRESENT(TRIEDLEV_DIAG)) THEN
!       TRIEDLEV_DIAG(ICL) = 1.
!    ENDIF

    UCU(ICMIN:) = 0.    ! This change makes cumulus friction
    VCU(ICMIN:) = 0.    ! correct.

    IF ( ICL > ICMIN ) THEN
   
       !--------------------------!
       !!  MAIN CLOUD SUBROUTINE !!
       !--------------------------!
 
       IC = ICL

       RC(IC) = 0

       ALM   = 0.
       TRG   = MIN(1.,(QOI(K)/QST(K)-RHMN)/(RHMX-RHMN))

       !F4 should ramp from 0 at SIG=AUTORAMPB to 1 at SIG=AUTORAMPB-0.2
       F4  = MIN(1.0, MAX( 0.0,(AUTORAMPB-SIGE(IC))/0.2 )) 
     
       if ( SIGE(IC) >= 0.5 ) then
          F5 = 1.0
       else
          F5 = 1.0 - 2.*CO_ZDEP *( 0.5 - SIGE(IC) )
          F5 = MAX( F5 , 0.0 )
       endif

       IF (TRG <= 1.0E-5) THEN    ! TRIGGER  =========>>
          RC(IC) = 7
       ENDIF

       if (RC(IC) == 0) then

          !RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL

          POI_c = POI
          QOI_c = QOI
          POI_c(K) =  POI_c(K) + TPERT
          QOI_c(K) =  QOI_c(K) + QPERT

          ZET(K+1) = 0.
          SHT(K+1) = CONS_CP*POI_c(K)*PRJ(K+1)
          DO L=K,IC,-1
             QOL(L)  = MIN(QST(L)*RHMAX,QOI_c(L))
             QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
             SSL(L)  = CONS_CP*PRJ(L+1)*POI_c(L) + CONS_GRAV*ZET(L+1)
             HOL(L)  = SSL(L) + QOL(L)*CONS_ALHL
             HST(L)  = SSL(L) + QST(L)*CONS_ALHL
             TEM     = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
             ZET(L)  = ZET(L+1) + TEM
             ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
          ENDDO

          DO L=IC+1,K
             TEM  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
             SHT(L)  = SSL(L-1) + TEM*(SSL(L)-SSL(L-1)) 
             QHT(L)  = .5*(QOL(L)+QOL(L-1))
          ENDDO

          !SMOOTH HSTAR W/ 1-2-1 Filter
          if ( SMOOTH_HST ) then
             HSTx = HST(IC) ! save for later
             DO L=K-1,IC+1,-1
                HST(L) = 0.25*(HST(L+1)+HST(L-1))+0.5*HST(L)
             END DO
             DO L=IC,IC
                HST(L) = 0.5*HST(L+1)+0.5*HST(L)
             END DO
          endif
     
          !CALCULATE LAMBDA, ETA, AND WORKFUNCTION
          LAMBDA_MIN = .2/MXDIAM
          LAMBDA_MAX = .2/  200. 

          IF (HOL(K) <= HST(IC)) THEN   ! CANNOT REACH IC LEVEL  ======>>
             RC(IC) = 1
          ENDIF
 
          if (RC(IC) == 0) then

             !LAMBDA CALCULATION: MS-A18
             TEM  =       (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
             DO L=IC+1,K-1
                TEM = TEM + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
             ENDDO

             IF (TEM <= 0.0) THEN         ! NO VALID LAMBDA  ============>>
                RC(IC) = 2
             ENDIF

             if (RC(IC) == 0) then

                ALM     = (HOL(K)-HST(IC)) / TEM

                IF (ALM > LAMBDA_MAX) THEN
                   RC(IC) = 3
                ENDIF
   
                if (RC(IC) == 0) then

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
                      HCC(L) = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
                      TEM    = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
                      EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
                      WFN    = WFN + (TEM - EHT(L)*HST(L))*GAM(L)
                   ENDDO
                   HCC(IC) = HST(IC)*ETA(IC)
                   WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)

                   !VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
                   BK3(K)   = 0.0
                   BK2(K)   = 0.0
                   BKE(K)   = 0.0
                   HCLD(K)  = HOL(K)
                   DO L=K-1,IC,-1
                      HCLD(L) = ( ETA(L+1)*HCLD(L+1)   +      & 
                               (ETA(L) - ETA(L+1))*HOL(L) ) / ETA(L)
                      TEM     = (HCLD(L)-HST(L) )*(ZET(L)-ZET(L+1))/ (1.0+LBCP*DQQ(L))
                      BKE(L)  = BKE(L+1) + CONS_GRAV * TEM / ( CONS_CP*PRJ(L+1)*POI(L) )
                      BK2(L)  = BK2(L+1) + CONS_GRAV * MAX(TEM,0.0) / ( CONS_CP*PRJ(L+1)*POI(L) )
                      BK3(L)  = BK3(L+1) + CONS_GRAV * MIN(TEM,0.0) / ( CONS_CP*PRJ(L+1)*POI(L) )
                      CVW(L) = SQRT(  2.0* MAX( BK2(L) , 0.0 )  )    
                   ENDDO

                   !ALPHA CALCULATION 
                   IF ( ZET(IC) <  2000. ) THEN
                      RASAL(IC) = RASAL1
                   ENDIF
                   IF ( ZET(IC)  >= 2000. ) THEN 
                      RASAL(IC) = RASAL1 + (RASAL2-RASAL1)*(ZET(IC) - 2000.)/8000.
                   ENDIF
                   RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
                   RASAL(IC) = DT / RASAL(IC)

                   CVW1 = 0.0
                   CVW1(IC:K) = MAX(  CVW(IC:K) , 1.00 )

                   !TEST FOR CRITICAL WORK FUNCTION
                   CALL ACRITN(POL(IC), PRS(K), ACR, ACRITFAC)

                   IF (WFN <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION ======>>
                      RC(IC) = 4
                   ENDIF

                   if (RC(IC) == 0) then

                      IF (DO_TRACERS) THEN 
                         DO ITR=1,ITRCR
                            !Scavenging of the below cloud tracer
                            DELZKM = ( ZET(IC) - ZET(K) ) / 1000.
                            FNOSCAV = max( min(exp(- FSCAV_(ITR) * DELZKM),1.), 0.)
                            XHT(ITR) = XOI(K,ITR) * FNOSCAV
                         END DO
                      END IF

                      WLQ     = QOL(K)
                      UHT     = UOI(K)
                      VHT     = VOI(K)
                      RNN(K)  = 0.
                      CLL0(K) = 0.

                      DO L=K-1,IC,-1
                         TEM   = ETA(L) - ETA(L+1)
                         WLQ   = WLQ + TEM * QOL(L)
                         UHT   = UHT + TEM * UOI(L)
                         VHT   = VHT + TEM * VOI(L)

                         IF (DO_TRACERS)  THEN  
                            DO ITR=1,ITRCR
                               !Scavenging of the entrained tracer.  Updates transported tracer mass.
                               DELZKM = ( ZET(IC) - ZET(L+1) ) / 1000.
                               FNOSCAV = max( min(exp(- FSCAV_(ITR) * DELZKM),1.), 0.)
                               XHT(ITR) = XHT(ITR) + TEM * XOI(L,ITR) * FNOSCAV
                            END DO
                         END IF

                         !How much condensate (CLI) is present here? 
                         IF (L>IC) THEN
                            TX2   = 0.5*(QST(L)+QST(L-1))*ETA(L)
                            TX3   = 0.5*(HST(L)+HST(L-1))*ETA(L)
                            QCC   = TX2 + GM1(L)*(HCC(L)-TX3)
                            CLL0(L)  = (WLQ-QCC)
                         ELSE
                            CLL0(L)   = (WLQ-QST(IC)*ETA(IC))
                         ENDIF

                         CLL0(L)    = MAX( CLL0(L) , 0.00 )

                         CLI  = CLL0(L) / ETA(L)  ! condensate (kg/kg)
                         TE_A = POI(L)*PRH(L)     ! Temperature (K)

                         CALL SUNDQ3_ICE( TE_A,  SDQV2, SDQV3, SDQVT1, F2 , F3 )

                         !F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                         C00_X  =  CO_AUTO * F2 * F3  * F4 !* F5
                         CLI_CRIT_X =  CLI_CRIT / ( F2 * F3 )
       
                         RATE = C00_X * ( 1.0 - EXP( -(CLI)**2 / CLI_CRIT_X**2 ) )

                         !Floor conv. vel. since, dont trust it at low values
                         CVW_X     = MAX( CVW1(L) , 1.00 )

                         DT_LYR  = ( ZET(L)-ZET(L+1) )/CVW_X ! l.h.s. DT_LYR => time in layer (L,L+1)

                         CLOSS   = CLL0(L) * RATE * DT_LYR
                         CLOSS   = MIN( CLOSS , CLL0(L) )

                         CLL0(L) = CLL0(L) - CLOSS

                         IF (CLOSS>0.) then
                            WLQ = WLQ - CLOSS
                            RNN(L) = CLOSS 
                         else
                            RNN(L) = 0.
                         ENDIF

                      ENDDO

                      WLQ = WLQ - QST(IC)*ETA(IC)

                      !CALCULATE GAMMAS AND KERNEL
                      GMS(K) =          (SHT(K)-SSL(K))*PRI(K)          ! MS-A30 (W/O CONS_GRAV)
                      GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*CONS_ALHL     ! MS-A31 (W/O CONS_GRAV)
                      AKM    = GMH(K)*GAM(K-1)*DPB(K-1)                 ! MS-A37 (W/O CONS_GRAV)

                      TX2     = GMH(K)
                      DO L=K-1,IC+1,-1
                         GMS(L) = ( ETA(L  )*(SHT(L)-SSL(L  ))   & 
                                  + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
                         GMH(L) = GMS(L) + ( ETA(L  )*(QHT(L)-QOL(L  ))   &
                                         +   ETA(L+1)*(QOL(L)-QHT(L+1)) )*CONS_ALHL*PRI(L)
                         TX2 = TX2 + (ETA(L) - ETA(L+1)) * GMH(L)
                         AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2*GHT(L)

                      ENDDO

                      GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
                      AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)

                      GMH(IC) =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*CONS_ALHL &
                                     + ETA(IC  )*(HST(IC)-HOL(IC  ))     )*PRI(IC)

                      if ( SMOOTH_HST ) then
                         GMHx    =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*CONS_ALHL &
                                     + ETA(IC  )*(HSTx   -HOL(IC  ))     )*PRI(IC)
                      endif

                      !CLOUD BASE MASS FLUX

                      IF (AKM >= 0.0 .OR. WLQ < 0.0)  THEN  !  =========>
                         RC(IC) = 5
                      ENDIF

                      if (RC(IC) == 0) then

                         WFN = - (WFN-ACR)/AKM ! MS-A39 MASS-FLUX IN Pa/step
                         WFN = MIN( ( RASAL(IC)*TRG*TOKI )*WFN  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))
     
                         !CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                         WFNOG    = WFN*GRAVI
                         TEM      = WFN*GRAVI
                         CLL (IC) = CLL (IC) + WLQ*TEM           ! (kg/m^2/step)
                         RMF (IC) = RMF (IC) +     TEM           ! (kg/m^2/step)
                         RMFD(IC) = RMFD(IC) +     TEM * ETA(IC) ! (kg/m^2/step)

                         DO L=IC+1,K
                            RMFP(L) = TEM * ETA(L)                 ! (kg/m^2/step)
!                            RMFC(L) = RMFC(L)  +  RMFP(L)          ! (kg/m^2/step)
        
                            IF ( CVW1(L) > 0.0 ) THEN
                               UPDFRP(L) = rmfp(L) * (DDT/DAYLEN) * 1000. / ( CVW1(L) * PRS(L) )
                            ELSE 
                               UPDFRP(L) = 0.0       
                            ENDIF

!                            CLLI(L) = CLL0(L)/ETA(L)  ! current cloud; incloud condensate        
!                            CLLB(L) = CLLB(L) +  UPDFRP(L) * CLLI(L) !  cumulative grid mean convective condensate        

                            UPDFRC(L) =  UPDFRC(L) +  UPDFRP(L)      
 
                         ENDDO

                         !THETA AND Q CHANGE DUE TO CLOUD TYPE IC

                         DO L=IC,K
                            RNS(L) = RNS(L) + RNN(L)*TEM ! (kg/m^2/step)
                            GMH(L) = GMH(L) * WFN
                            GMS(L) = GMS(L) * WFN
                            QOI(L) = QOI(L) + (GMH(L) - GMS(L)) * ALHI
                            POI(L) = POI(L) + GMS(L)*PKI(L)*CPI
                            QST(L) = QST(L) + GMS(L)*BET(L)*CPI
                         ENDDO
      
                         if ( SMOOTH_HST ) then
                            GMHx    = GMHx * WFN
                            dQx     = (GMHx - GMH(IC)) * ALHI
                            RNS(IC) = RNS(IC) + dQx / ( PRI(IC)*CONS_GRAV )
                         endif

                         IF (DO_TRACERS) THEN
                            WFN     = WFN*0.5 *1.0 !*FRICFAC*0.5
                            TEM     = WFN*PRI(K)

                            DO ITR=1,ITRCR 
                               XCU(K,ITR) =  XCU(K,ITR) + TEM * (XOI(K-1,ITR) - XOI(K,ITR))
                            END DO

                            DO ITR=1,ITRCR 
                               DO L=K-1,IC+1,-1
                                  TEM    = WFN*PRI(L)
                                  XCU(L,ITR) = XCU(L,ITR) + TEM *                        &
                                              ( (XOI(L-1,ITR) - XOI(L,ITR  )) * ETA(L)   &
                                              + (XOI(L,ITR  ) - XOI(L+1,ITR)) * ETA(L+1) )
                               ENDDO
                            ENDDO

                            TEM     = WFN*PRI(IC)

                            DO ITR=1,ITRCR 
                               XCU(IC,ITR) = XCU(IC,ITR) +  &
                                             (2.*(XHT(ITR) - XOI(IC,ITR)*(ETA(IC)-ETA(IC+1))) & 
                                             - (XOI(IC,ITR)+XOI(IC+1,ITR))*ETA(IC+1))*TEM
                            ENDDO

                            DO ITR=1,ITRCR 
                               DO L=IC,K
                                  XOI(L,ITR) = XOI(L,ITR) + XCU(L,ITR)
                               ENDDO
                            ENDDO
                         else 

                            WFN = WFN*0.5 *1.0           !*FRICFAC*0.5
 
                         endif

                         !CUMULUS FRICTION
                         IF (FRICFAC <= 0.0) THEN
                            RC(IC) = 6
                         ENDIF

                         if (RC(IC) == 0) then

                            WFN     = WFN*FRICFAC*EXP( -ALM / FRICLAMBDA )
                            TEM     = WFN*PRI(K)

                            UCU(K)  = UCU(K) + TEM * (UOI(K-1) - UOI(K))
                            VCU(K)  = VCU(K) + TEM * (VOI(K-1) - VOI(K))

                            DO L=K-1,IC+1,-1
                               TEM    = WFN*PRI(L)
                               UCU(L) = UCU(L) + TEM *                        &
                                           ( (UOI(L-1) - UOI(L  )) * ETA(L  ) &
                                           + (UOI(L  ) - UOI(L+1)) * ETA(L+1) )
                               VCU(L) = VCU(L) + TEM *                        &
                                           ( (VOI(L-1) - VOI(L  )) * ETA(L)   &
                                           + (VOI(L  ) - VOI(L+1)) * ETA(L+1) )
                            ENDDO

                            TEM     = WFN*PRI(IC)

                            UCU(IC) = UCU(IC) + (2.*(UHT - UOI(IC)*(ETA(IC)-ETA(IC+1))) & 
                                       - (UOI(IC)+UOI(IC+1))*ETA(IC+1))*TEM
                            VCU(IC) = VCU(IC) + (2.*(VHT - VOI(IC)*(ETA(IC)-ETA(IC+1))) & 
                                       - (VOI(IC)+VOI(IC+1))*ETA(IC+1))*TEM


                            DISSK0(IC) = ETA(IC)* CONS_GRAV * WFNOG * PRI(IC) * 0.5 & 
                                       *( (UHT/ETA(IC)  - UOI(IC) )**2    & 
                                        + (VHT/ETA(IC)  - VOI(IC) )**2 )  
                            DO L=IC,K
                             UOI(L) = UOI(L) + UCU(L)
                             VOI(L) = VOI(L) + VCU(L)
                            ENDDO

                            RC(IC) = 10

                         endif
                      endif
                   endif
                endif
             endif
          endif
       endif
       !--------------------------!
       !!  MAIN CLOUD SUBROUTINE !!
       !--------------------------!


    ENDIF

!    ENTLAM(ICL) = ALM

 ENDDO

 IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

    !RNEVP
    ZET(K+1) = 0
    DO L=K,ICMIN,-1
       TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
       ZET(L)  = ZET(L+1) + TEM
    ENDDO

    DO L=ICMIN,K
       TEM    = PRI(L)*CONS_GRAV
       CNV_PRC3(L) = RNS(L)*TEM + CNV_PRC3TMP(L)
    ENDDO

    !If hst is smoothed then adjusted precips may be negative
    if ( SMOOTH_HST ) then
       DO L=ICMIN,K
          if ( CNV_PRC3(L) < 0. ) then
             QOI(L)        = QOI(L) +  CNV_PRC3(L)
             POI(L)        = POI(L) -  CNV_PRC3(L) * (CONS_ALHL/CONS_CP) / PRJ(L+1)
             CNV_PRC3(L) = 0.
          endif
       END DO
    endif

    !STRAP (FINAL = 1) SUBROUTINE
    !----------------------------
    THO(ICMIN:K-1) = POI(ICMIN:K-1)
    QHO(ICMIN:K-1) = QOI(ICMIN:K-1)
    UHO(ICMIN:K-1) = UOI(ICMIN:K-1)
    VHO(ICMIN:K-1) = VOI(ICMIN:K-1)
    CNV_UPDFRC(ICMIN:K-1)   =  UPDFRC(ICMIN:K-1) + CNV_UPDFRCTMP(ICMIN:K-1)

!    CNV_CVW   (ICMIN:K-1)   =     CVW1(ICMIN:K-1)
!    CNV_QC(ICMIN:K-1)       =  CLLB(ICMIN:K-1)

    !De-strap tendencies from RAS specify weighting "SHAPE"
    WGHT1   = WGT1(:)

    !Scale properly by layer masses
    wght0 = 0.
    DO L=K,K0 
       wght0 = wght0 + WGHT1(L)* ( PLE(L+1) - PLE(L) )
    END DO
         
    wght0 = ( PRS(K+1)   - PRS(K)  )/wght0
    WGHT1  = wght0 * WGHT1

    DO L=K,K0 
       THO(L) =  THO(L) + WGHT1(L)*(POI(K) - POI_SV(K))
       QHO(L) =  QHO(L) + WGHT1(L)*(QOI(K) - QOI_SV(K))
       UHO(L) =  UHO(L) + WGHT1(L)*(UOI(K) - UOI_SV(K))
       VHO(L) =  VHO(L) + WGHT1(L)*(VOI(K) - VOI_SV(K))
    END DO

    IF (DO_TRACERS) THEN 
       XHO(ICMIN:K-1,:) = XOI(ICMIN:K-1,:) 
       DO ITR=1,ITRCR
          DO L=K,K0
             XHO(L,ITR) =  XHO(L,ITR) + WGHT1(L)*(XOI(K,ITR) - XOI_SV(K,ITR))
          END DO
       END DO
    END IF

!    FLX (ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN                     !  (KG/m^2/s @ CLOUD BASE)
    FLXD(ICMIN:K) = RMFD(ICMIN:K) * DDT/DAYLEN + FLXDTMP(ICMIN:K)  !  (KG/m^2/s @ CLOUD TOP)
!    FLXC(ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN                     !  (KG/m^2/s @ CLOUD TOP)
    CLW (ICMIN:K) = CLL (ICMIN:K) * DDT/DAYLEN + CLWTMP(ICMIN:K)   !  (KG/m^2/s )

    if ( PRESENT( DISSKE )) THEN
       DISSKE(ICMIN:K-1) = DISSK0(ICMIN:K-1) * DDT/DAYLEN + DISSKETMP(ICMIN:K-1)
    endif

!    FLX (1:ICMIN-1) = 0.
    FLXD(1:ICMIN-1) = 0.
!    FLXC(1:ICMIN-1) = 0.
    CLW (1:ICMIN-1) = 0.
 
    IF ( K < K0 ) THEN 
!       FLX (K:K0) = 0.
       FLXD(K:K0) = 0.
!       FLXC(K:K0) = 0.
       CLW (K:K0) = 0.
    END IF

    IF (ALLOCATED(ICL_V)) THEN
       DEALLOCATE( ICL_V )
    ENDIF

 ELSE

    !STRAP (FINAL = 2) SUBROUTINE
!    FLX (:) = 0.
    FLXD(:) = 0.
!    FLXC(:) = 0.
    CLW (:) = 0.

    !----------------------------

 ENDIF 

 IF (ALLOCATED(ICL_V)) THEN
    DEALLOCATE(ICL_V)
 ENDIF

      
END SUBROUTINE RASE


! Subroutines
!------------

SUBROUTINE ACRITN( PL, PLB, ACR, ACRITFAC )
      
 IMPLICIT NONE

 REAL(8), INTENT(IN ) :: PL, PLB, ACRITFAC
 REAL(8), INTENT(OUT) :: ACR

 INTEGER :: IWK
      
 REAL(8), PARAMETER :: PH(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, &
                             550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)

 REAL(8), PARAMETER :: A(15)=(/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677, &
                             0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
                             0.0553, 0.0445, 0.0633  /) 

 IWK = nint(PL * 0.02 - 0.999999999)

 IF (IWK .GT. 1 .AND. IWK .LE. 15) THEN
    ACR = A(IWK-1) + (PL-PH(IWK-1))*.02*(A(IWK)-A(IWK-1))
 ELSEIF(IWK > 15) THEN
    ACR = A(15)
 ELSE
    ACR = A(1)
 ENDIF

 ACR = ACRITFAC  * ACR * (PLB - PL)

ENDSUBROUTINE ACRITN

SUBROUTINE SUNDQ3_ICE( TEMP,RATE2,RATE3,TE1, F2, F3)

IMPLICIT NONE

REAL(8), INTENT( IN) :: TEMP,RATE2,RATE3,TE1
REAL(8), INTENT(OUT) :: F2, F3

REAL(8) :: XX, YY,TE0,TE2,JUMP1  !,RATE2,RATE3,TE1

 TE0=273.
 TE2=200.
 JUMP1=  (RATE2-1.0) / ( ( TE0-TE1 )**0.333 ) 

 ! Ice - phase treatment
 IF ( TEMP .GE. TE0 ) THEN 
    F2   = 1.0
    F3   = 1.0
 ENDIF

 IF ( ( TEMP .GE. TE1 ) .AND. ( TEMP .LT. TE0 ) )THEN 
    F2   = 1.0 + JUMP1 * (( TE0 - TEMP )**0.3333)
    F3   = 1.0
 ENDIF

 IF ( TEMP .LT. TE1 ) THEN 
    F2   = RATE2 + (RATE3-RATE2)*(TE1-TEMP)/(TE1-TE2)
    F3   = 1.0
 ENDIF

 IF ( F2 .GT. 27.0 ) THEN
    F2 = 27.0
 ENDIF

endsubroutine sundq3_ice

subroutine DQSAT_RAS(DQSi,QSSi,TEMP,PLO,im,jm,lm,ESTBLX,CONS_H2OMW, CONS_AIRMW)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 integer :: im,jm,lm
 real(8), dimension(im,jm,lm) :: TEMP, PLO
 real(8) :: ESTBLX(:)
 real(4) :: CONS_H2OMW, CONS_AIRMW

 !Outputs
 real(8), dimension(im,jm,lm) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8) :: ESFAC

 integer :: i, j, k

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 integer, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 ESFAC = CONS_H2OMW/CONS_AIRMW

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

end subroutine DQSAT_RAS

subroutine DQSATs_RAS(DQSi,QSSi,TEMP,PLO,ESTBLX,CONS_H2OMW,CONS_AIRMW)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 real(8) :: TEMP, PLO
 real(8) :: ESTBLX(:)
 real(4) :: CONS_H2OMW, CONS_AIRMW

 !Outputs
 real(8) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8) :: ESFAC

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 integer :: IT

 integer, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 ESFAC = CONS_H2OMW/CONS_AIRMW

 TL = TEMP
 PL = PLO

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

 DQSi = DQSAT
 QSSi = QSAT

end subroutine DQSATs_RAS


END MODULE CONVECTION
