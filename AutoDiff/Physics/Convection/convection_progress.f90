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
                  CO_AUTO        , &
                  PLE            , &
                  RASPARAMS      , &
                  CONS_CP        , &
                  CONS_ALHL      , &
                  CONS_GRAV      , &
                  CONS_RGAS      , &
                  CONS_H2OMW     , &
                  CONS_AIRMW     , &
                  CONS_VIREPS    , &
                  ESTBLX         , &
                  TAP_DUMMY      )

 IMPLICIT NONE

 !INPUTS
 INTEGER, INTENT(IN) :: K0, ICMIN, KCBL, SEEDRAS
 real(8), INTENT(IN) :: DT, FRLAND, TS, CO_AUTO, TAP_DUMMY

 real(8), DIMENSION(K0  ), INTENT(IN) :: WGT0, WGT1
 real(8), DIMENSION(K0+1), INTENT(IN) :: PLE, SIGE
 real(8), DIMENSION(:)   , INTENT(IN) :: RASPARAMS, ESTBLX

 REAL*4, INTENT(IN) :: CONS_CP, CONS_ALHL
 REAL*4, INTENT(IN) :: CONS_GRAV, CONS_RGAS
 REAL*4, INTENT(IN) :: CONS_H2OMW, CONS_AIRMW, CONS_VIREPS

 !PROGNOSTIC
 real(8), DIMENSION (K0  ), INTENT(INOUT) ::  THO, QHO

 !LOCALS
 INTEGER :: IC, CDET, K, LL, L, ITR

 real(8), PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
 real(8), PARAMETER :: CBL_QPERT = 0.0, CBL_TPERT = 1.0
 real(8), PARAMETER :: CBL_TPERT_MXOCN = 2.0, CBL_TPERT_MXLND = 4.0

 !Constants
 real(8) :: CPI, DDT, LBCP, RHMX, RHMN

 !Rasparams
 real(8) :: MAXDALLOWED

 !Variables constant for the column - not passed to cloude
 real(8) :: PRCBL, rndu
 real(8), DIMENSION(K0)   :: TEM8
 real(8), DIMENSION(K0+1) :: PKE, ZLE
 real(8), DIMENSION(K0)   :: QSS, DQS, TEMPf, Pf, PK, ZLO
 real(8), DIMENSION(K0)   :: POI_SV, QOI_SV
 real(8), DIMENSION(K0)   :: WGHT, WGHT0, WGHT1, MASSF

 !Variables constant for the column - passed to cloude
 real(8) :: TPERT, QPERT, MXDIAM
 real(8), DIMENSION(K0+1) :: PRJ, PRS
 real(8), DIMENSION(K0)   :: DQQ, BET, GAM, POL, GM1
 real(8), DIMENSION(K0)   :: PRH, PRI, GHT, DPT, DPB, PKI

 !Variables output by cloude
 real(8), DIMENSION(K0)   :: RNS, UPDFRC

 !Variables seen by all parts
 real(8), DIMENSION(K0)   :: POI0, QOI0, QST0
 real(8), DIMENSION(K0)   :: POI1, QOI1
 
 real(8), DIMENSION(K0+2,K0) :: POIC, QOIC, QSTC
 real(8), DIMENSION(K0+2,K0) :: RNSC, UPDFRCC

 real(8), DIMENSION(K0) :: POICa, QOICa, QSTCa
 real(8), DIMENSION(K0) :: RNSCa, UPDFRCCa
 real(8), DIMENSION(K0) :: POICb, QOICb, QSTCb
 real(8), DIMENSION(K0) :: RNSCb, UPDFRCCb

! BEGIN CALCULATIONS
!-------------------

 !Initialize local arrays
 PRJ = 0.0
 PRS = 0.0
 ZLE = 0.0
 QSS = 0.0
 DQS = 0.0
 TEMPf = 0.0
 Pf = 0.0
 PK = 0.0
 ZLO = 0.0
 POI_SV = 0.0
 QOI_SV = 0.0
 DQQ = 0.0
 BET = 0.0
 GAM = 0.0
 PRH = 0.0
 PRI = 0.0
 GHT = 0.0
 DPT = 0.0
 DPB = 0.0
 PKI = 0.0
 RNS = 0.0
 POL = 0.0
 QST0 = 0.0
 GM1 = 0.0
 UPDFRC = 0.0
 WGHT = 0.0
 WGHT1 = 0.0
 MASSF = 0.0
 
 !Rasparams
 MAXDALLOWED  = RASPARAMS(23)    !  --- 23

 !Constants
 CPI   = 1.0/CONS_CP      
 DDT   = DAYLEN/DT
 LBCP  = CONS_ALHL*CPI

 RHMN         = RASPARAMS(24)    !  --- 24
 RHMX         = RASPARAMS(25)    !  --- 25


 !Get saturation specific humidity and gradient wrt to T
 PKE  = (PLE/1000.)**(CONS_RGAS/CONS_CP)
 Pf   = 0.5*(PLE(1:K0) +  PLE(2:K0+1  ) )
 PK   = (Pf/1000.)**(CONS_RGAS/CONS_CP)
 TEMPf = THO(:)*PK

 ZLE(K0+1) = 0.
 do L=K0,1,-1
    ZLE(L) = THO(L)   * (1.+CONS_VIREPS*QHO(L))
    ZLO(L) = ZLE(L+1) + (CONS_CP/CONS_GRAV)*( PKE(L+1)-PK (L  ) ) * ZLE(L)
    ZLE(L) = ZLO(L)   + (CONS_CP/CONS_GRAV)*( PK (L)  -PKE(L)   ) * ZLE(L)
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

 call DQSAT_RAS(DQS,QSS,TEMPf,Pf,K0,ESTBLX,CONS_H2OMW,CONS_AIRMW)

 !CALL FINDBASE
 K = KCBL

 !STRAP FINAL=0 SUBROUTINE
 !------------------
 do LL=icmin,k+1
    PRJ(LL) = PKE(LL)
 enddo

 poi0=0.        ! These initialized here in order not to confuse Valgrind debugger
 qoi0=0.        ! Do not believe it actually makes any difference.
  

 PRS(ICMIN:K0+1) = PLE(ICMIN:K0+1)
 POI0(ICMIN:K)   = THO(ICMIN:K)
 QOI0(ICMIN:K)   = QHO(ICMIN:K)

 QST0(ICMIN:K) = QSS(ICMIN:K)
 DQQ(ICMIN:K) = DQS(ICMIN:K)

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

    POI0(K) = 0.
    QOI0(K) = 0.

    !SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
    WGHT = 0.
    DO L=K,K0
       WGHT(L) = MASSF(L) * ( PLE(L+1) - PLE(L) ) / ( PRS(K+1) - PRS(K) )
    END DO      

    DO L=K,K0
       POI0(K) = POI0(K) + WGHT(L)*THO(L)
       QOI0(K) = QOI0(K) + WGHT(L)*QHO(L)
    ENDDO

    call DQSATs_RAS(DQQ(K), QST0(K), POI0(K)*PRH(K),POL(K),ESTBLX,CONS_H2OMW,CONS_AIRMW)

 endif

 rndu = max( seedras/1000000., 1e-6 )
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

 POI_SV = POI0
 QOI_SV = QOI0

 !END STRAP (FINAL=0) SUBROUTINE
 !------------------------------


 !MAIN CLOUD SUBROUTINE
 !---------------------

 POIC = 0.0
 QOIC = 0.0
 QSTC = 0.0

 POIC(KCBL+2,:) = POI0
 QOIC(KCBL+2,:) = QOI0
 QSTC(KCBL+2,:) = QST0

 !Begin cloud loop
 DO IC = KCBL+1,ICMIN+1,-1

    CDET = IC

    POICa = POIC(IC+1,:)
    QOICa = QOIC(IC+1,:)
    QSTCa = QSTC(IC+1,:)

    !IF (MIN(1.,(QOICa(K)/QSTCa(K)-RHMN)/(RHMX-RHMN)) > 1.0E-5) THEN  

    CALL CLOUDE( CDET, ICMIN, K, K0, DT, SIGE, CO_AUTO, CONS_GRAV, CONS_CP, CONS_ALHL,                     &
                 RASPARAMS,                                                                &
                 CPI, DDT, LBCP, TPERT, QPERT, MXDIAM, PRJ, PRS, DQQ,                                   &
                 BET, GAM, POL, GM1, PRH, PRI, GHT, DPT, DPB, PKI,                                      &
                 POICa, QOICa, QSTCa,                  &
                 POIC(IC  ,:), QOIC(IC  ,:), QSTC(IC  ,:) )

    !ELSE

!POIC(IC  ,:) = POICa
!QOIC(IC  ,:) = QOICa
!QSTC(IC  ,:) = QSTCa

    !ENDIF


 endDO

 POI1 = POIC(ICMIN+1,:)
 QOI1 = QOIC(ICMIN+1,:)

 !---------------------
 !MAIN CLOUD SUBROUTINE


! IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

    !STRAP (FINAL = 1) SUBROUTINE
    !----------------------------

    !Scale properly by layer masses
    wght0 = 0.
    DO L=K,K0 
       wght0(L) = wght0(L-1) + WGT1(L)* ( PLE(L+1) - PLE(L) )
    END DO
    WGHT1  = WGT1 * ( PRS(K+1) - PRS(K) )/ wght0(K0)

    THO(ICMIN:K-1) = POI1(ICMIN:K-1)
    QHO(ICMIN:K-1) = QOI1(ICMIN:K-1)

    DO L=K,K0 
       THO(L) =  THO(L) + WGHT1(L)*(POI1(K) - POI_SV(K))
       QHO(L) =  QHO(L) + WGHT1(L)*(QOI1(K) - QOI_SV(K))
    END DO

! ENDIF 

 CONTAINS

 SUBROUTINE CLOUDE( CDET, ICMIN, K, K0, DT, SIGE, CO_AUTO, CONS_GRAV, CONS_CP, CONS_ALHL,    &
                    RASPARAMS,                                               &
                    CPI, DDT, LBCP, TPERT, QPERT, MXDIAM, PRJ, PRS, DQQ,                  &
                    BET, GAM, POL, GM1, PRH, PRI, GHT, DPT, DPB, PKI,                     &
                    POI_I, QOI_I, QST_I,                                    &
                    POI_O, QOI_O, QST_O )

 IMPLICIT NONE

 !INPUTS
 INTEGER, INTENT(IN) :: CDET, K0, ICMIN, K
 real(8), INTENT(IN) :: DT, CO_AUTO

 real(8), DIMENSION(K0+1), INTENT(IN) :: SIGE
 real(8), DIMENSION(:)   , INTENT(IN) :: RASPARAMS

 REAL*4, INTENT(IN) :: CONS_CP, CONS_ALHL
 REAL*4, INTENT(IN) :: CONS_GRAV
 real(8), INTENT(IN) :: CPI, DDT, LBCP

 real(8), INTENT(IN) :: TPERT, QPERT, MXDIAM

 real(8), INTENT(IN), DIMENSION(K0+1) :: PRJ, PRS
 real(8), INTENT(IN), DIMENSION(K0)   :: DQQ, BET, GAM, POL, GM1
 real(8), INTENT(IN), DIMENSION(K0)   :: PRH, PRI, GHT, DPT, DPB, PKI


 !OUTPUTS

 real(8), DIMENSION(K0)   :: CVW

 !PROGNOSTIC
 real(8), INTENT(IN), DIMENSION(K0)    :: POI_I, QOI_I, QST_I
 real(8), INTENT(OUT), DIMENSION(K0)   :: POI_O, QOI_O, QST_O

 !LOCALS
 INTEGER :: RC, KK, ITR, L

 real(8) :: CPBG, GRAVI, ALHI
 real(8), PARAMETER :: RHMAX   = 0.9999

 real(8) :: FRICFAC, CLI_CRIT, RASAL1, RASAL2 
 real(8) :: FRICLAMBDA
 real(8) :: SDQV2, SDQV3, SDQVT1
 real(8) :: ACRITFAC,  PBLFRAC, AUTORAMPB
 real(8) :: MAXDALLOWED, RHMN, RHMX

 real(8) :: LAMBDA_MIN, LAMBDA_MAX

 !Variables that vary within the cloude loop
 real(8) :: AKM2, ACR, ALM
 real(8) :: WFN2, WFN3, WFN4, WFN5, WFN6, TRG 
 real(8) :: TEM12, TEM13, TEM14, TEM15, TEM16
 real(8) :: WFNOG, TOKI, F4, WLQ1

 real(8), DIMENSION(K0+1) :: SHT, ZET
 real(8), DIMENSION(K0)   :: TEM2, TEM3, TEM4, TEM5, TEM6, TEM7
 real(8), DIMENSION(K0)   :: POI_c, QOI_c
 real(8), DIMENSION(K0)   :: UCU, VCU
 real(8), DIMENSION(K0)   :: SSL, RNN
 real(8), DIMENSION(K0)   :: GMS, GMS1, ETA, GMH, GMH1, EHT, HCC, BKE
 real(8), DIMENSION(K0)   :: HOL, HST, QOL, ZOL, HCLD
 real(8), DIMENSION(K0)   :: CLL0,CLL01,CLL02
 real(8), DIMENSION(K0)   :: RASAL, RASAL12, BK2
 real(8), DIMENSION(K0)   :: QHT
 real(8), DIMENSION(K0)   :: TX21, TX22, TX3, AKM
 real(8), DIMENSION(K0)   :: WFN1, QCC
 real(8), DIMENSION(K0)   :: CLI, TE_A, C00_X, CLI_CRIT_X
 real(8), DIMENSION(K0)   :: DT_LYR, RATE, CVW_X, F2, F3
 real(8), DIMENSION(K0)   :: CLOSS, CLOSS1, WLQ, CVW1

 GRAVI = 1.0/CONS_GRAV
 CPBG  = CONS_CP*GRAVI
 ALHI  = 1.0/CONS_ALHL

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
 RHMN         = RASPARAMS(24)    !  --- 24
 RHMX         = RASPARAMS(25)    !  --- 25

 LAMBDA_MIN = .2/MXDIAM
 LAMBDA_MAX = .2/  200. 

!LOOP

 AKM2 = 0.0
 ACR = 0.0
 ALM = 0.0
 WFN2 = 0.0
 WFN3 = 0.0
 WFN4 = 0.0
 WFN5 = 0.0
 WFN6 = 0.0
 TRG  = 0.0
 TEM12 = 0.0
 TEM13 = 0.0
 TEM14 = 0.0
 TEM15 = 0.0
 TEM16 = 0.0
 WFNOG = 0.0
 TOKI = 0.0
 F4 = 0.0
 WLQ = 0.0
 SHT = 0.0
 ZET = 0.0
 TEM2 = 0.0
 TEM3 = 0.0
 TEM4 = 0.0
 TEM5 = 0.0
 TEM6 = 0.0
 TEM7 = 0.0
 POI_c = 0.0
 QOI_c = 0.0
 UCU = 0.0
 VCU = 0.0
 SSL = 0.0
 RNN = 0.0
 GMS = 0.0
 GMS1 = 0.0
 ETA = 0.0
 GMH = 0.0
 GMH1 = 0.0
 EHT = 0.0
 HCC = 0.0
 BKE = 0.0
 HOL = 0.0
 HST = 0.0
 QOL = 0.0
 ZOL = 0.0
 HCLD = 0.0
 CLL0 = 0.0
 CLL01 = 0.0
 CLL02 = 0.0
 RASAL = 0.0
 RASAL12 = 0.0
 BK2 = 0.0
 QHT = 0.0
 TX21 = 0.0
 TX22 = 0.0
 TX3 = 0.0
 AKM = 0.0
 WFN1 = 0.0
 QCC = 0.0
 CLI = 0.0
 TE_A = 0.0
 C00_X = 0.0
 CLI_CRIT_X = 0.0
 DT_LYR = 0.0
 RATE = 0.0
 CVW_X = 0.0
 F2 = 0.0
 F3 = 0.0
 CLOSS = 0.0
 CLOSS1 = 0.0
 
    RC = 0

    ALM   = 0.
    TRG   = MIN(1.,(QOI_I(K)/QST_I(K)-RHMN)/(RHMX-RHMN))

    !F4 should ramp from 0 at SIG=AUTORAMPB to 1 at SIG=AUTORAMPB-0.2
    F4  = MIN(1.0, MAX( 0.0,(AUTORAMPB-SIGE(CDET))/0.2 )) 

    IF (TRG <= 1.0E-5) THEN    ! TRIGGER  =========>>
       RC = 1
    ENDIF

    if (RC == 0) then

       !RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL

       POI_c = POI_I
       QOI_c = QOI_I
       POI_c(K) =  POI_c(K) + TPERT
       QOI_c(K) =  QOI_c(K) + QPERT

       ZET(K+1) = 0.
       DO L=K,CDET,-1
          TEM2(L) = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
          ZET(L)  = ZET(L+1) + TEM2(L)
          ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
       ENDDO

       SHT(K+1) = CONS_CP*POI_c(K)*PRJ(K+1)
       DO L=K,CDET,-1
          QOL(L)  = MIN(QST_I(L)*RHMAX,QOI_c(L))
          QOL(L)  = MAX( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
          SSL(L)  = CONS_CP*PRJ(L+1)*POI_c(L) + CONS_GRAV*ZET(L+1)
          HOL(L)  = SSL(L) + QOL(L)*CONS_ALHL
          HST(L)  = SSL(L) + QST_I(L)*CONS_ALHL
       ENDDO

       DO L=CDET+1,K
          TEM3(L)  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
          SHT(L)  = SSL(L-1) + TEM3(L)*(SSL(L)-SSL(L-1)) 
          QHT(L)  = .5*(QOL(L)+QOL(L-1))
       ENDDO
     
       !CALCULATE WORKFUNCTION

       IF (HOL(K) <= HST(CDET)) THEN   ! CANNOT REACH CDET LEVEL  ======>>
          RC = 2
       ENDIF
 
       if (RC == 0) then

          !LAMBDA CALCULATION: MS-A18
          TEM4(CDET)  =       (HST(CDET)-HOL(CDET))*(ZOL(CDET)-ZET(CDET+1)) 
          DO L=CDET+1,K-1
             TEM4(L) = TEM4(L-1) + (HST(CDET)-HOL(L ))*(ZET(L )-ZET(L +1))
          ENDDO

          IF (TEM4(K-1) <= 0.0) THEN         ! NO VALID LAMBDA  ============>>
             RC = 3
          ENDIF

          if (RC == 0) then

             ALM     = (HOL(K)-HST(CDET)) / TEM4(CDET)

             IF (ALM > LAMBDA_MAX) THEN
                RC = 4
             ENDIF
   
             if (RC == 0) then

                TOKI = 1.0
   
                IF (ALM < LAMBDA_MIN) THEN
                   TOKI = ( ALM/LAMBDA_MIN )**2
                ENDIF

!ALM = (HOL(K)-HST(CDET))/

                !ETA CALCULATION: MS-A2
                DO L=CDET+1,K
                   ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
                ENDDO
                ETA(CDET) = 1.0 + ALM * (ZOL(CDET)-ZET(K))

                !WORKFUNCTION CALCULATION:  MS-A22
                WFN1(K)     = 0.0
                HCC(K)  = HOL(K)
                DO L=K-1,CDET+1,-1
                   HCC(L)  = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
                   TEM5(L) = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
                   EHT(L)  = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
                   WFN1(L)  = WFN1(L+1) + (TEM5(L) - EHT(L)*HST(L))*GAM(L)
                ENDDO
                HCC(CDET) = HST(CDET)*ETA(CDET)
                WFN2     = WFN1(CDET+1) + (HCC(CDET+1)-HST(CDET)*ETA(CDET+1))*GAM(CDET)*DPB(CDET)

                !VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
                BK2(K)   = 0.0
                BKE(K)   = 0.0
                HCLD(K)  = HOL(K)
                DO L=K-1,CDET,-1
                   HCLD(L) = ( ETA(L+1)*HCLD(L+1)   +      & 
                            (ETA(L) - ETA(L+1))*HOL(L) ) / ETA(L)
                   TEM6(L)     = (HCLD(L)-HST(L) )*(ZET(L)-ZET(L+1))/ (1.0+LBCP*DQQ(L))
                   BKE(L)  = BKE(L+1) + CONS_GRAV * TEM6(L) / ( CONS_CP*PRJ(L+1)*POI_I(L) )
                   BK2(L)  = BK2(L+1) + CONS_GRAV * MAX(TEM6(L),0.0) / ( CONS_CP*PRJ(L+1)*POI_I(L) )
                   CVW1(L) = SQRT(  2.0* MAX( BK2(L) , 0.0 )  )    
                ENDDO

                !ALPHA CALCULATION 
                IF ( ZET(CDET) <  2000. ) THEN
                   RASAL(CDET) = RASAL1
                ENDIF
                IF ( ZET(CDET)  >= 2000. ) THEN 
                   RASAL(CDET) = RASAL1 + (RASAL2-RASAL1)*(ZET(CDET) - 2000.)/8000.
                ENDIF
                RASAL12(CDET) = DT / MIN( RASAL(CDET) , 1.0e5 )

                DO KK = CDET,K
                   CVW(KK) = MAX(  CVW1(KK) , 1.00 )
                ENDDO

                !TEST FOR CRITICAL WORK FUNCTION
                CALL ACRITN(POL(CDET), PRS(K), ACR, ACRITFAC)

                IF (WFN2 <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION ======>>
                   RC = 5
                ENDIF

                if (RC == 0) then

                   WLQ(K)     = QOL(K)
                   RNN(K)  = 0.
                   CLL0(K) = 0.

                   DO L=K-1,CDET,-1
                      TEM7(L)   = ETA(L) - ETA(L+1)
                      WLQ(L)   = WLQ(L+1) + TEM7(L) * QOL(L)
                      !How much condensate (CLI) is present here? 
                      IF (L>CDET) THEN
                         TX21(L)   = 0.5*(QST_I(L)+QST_I(L-1))*ETA(L)
                         TX3(L)   = 0.5*(HST(L)+HST(L-1))*ETA(L)
                         QCC(L)   = TX21(L) + GM1(L)*(HCC(L)-TX3(L))
                         CLL0(L)  = (WLQ(L)-QCC(L))
                      ELSE
                         CLL0(L)   = (WLQ(L)-QST_I(CDET)*ETA(CDET))
                      ENDIF

                      CLL01(L)    = MAX( CLL0(L) , 0.00 )

                      CLI(L)  = CLL01(L) / ETA(L)  ! condensate (kg/kg)
                      TE_A(L) = POI_I(L)*PRH(L)     ! Temperature (K)

                      CALL SUNDQ3_ICE( TE_A(L),  SDQV2, SDQV3, SDQVT1, F2(L) , F3(L) )

                      !F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                      C00_X(L)  =  CO_AUTO * F2(L) * F3(L)  * F4 !* F5
                      CLI_CRIT_X(L) =  CLI_CRIT / ( F2(L) * F3(L) )
       
                      RATE(L) = C00_X(L) * ( 1.0 - EXP( -(CLI(L))**2 / CLI_CRIT_X(L)**2 ) )

                      !Floor conv. vel. since, dont trust it at low values
                      CVW_X(L)     = MAX( CVW(L) , 1.00 )

                      DT_LYR(L)  = ( ZET(L)-ZET(L+1) )/CVW_X(L) ! l.h.s. DT_LYR(L) => time in layer (L,L+1)

                      CLOSS(L)   = CLL01(L) * RATE(L) * DT_LYR(L)
                      CLOSS1(L)   = MIN( CLOSS(L) , CLL01(L) )

                      CLL02(L) = CLL01(L) - CLOSS1(L)

                      IF (CLOSS1(L) > 0.) then
                         WLQ(L) = WLQ(L) - CLOSS1(L)
                         RNN(L) = CLOSS1(L) 
                      else
                         RNN(L) = 0.
                      ENDIF

                   ENDDO

                   WLQ1 = WLQ(CDET) - QST_I(CDET)*ETA(CDET)

                   !CALCULATE GAMMAS AND KERNEL
                   GMS(K) =          (SHT(K)-SSL(K))*PRI(K)          ! MS-A30 (W/O GRAV)
                   GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*CONS_ALHL     ! MS-A31 (W/O GRAV)
                   AKM(K)    = GMH(K)*GAM(K-1)*DPB(K-1)                 ! MS-A37 (W/O GRAV)

                   TX22(K)     = GMH(K)
                   DO L=K-1,CDET+1,-1
                      GMS(L) = ( ETA(L  )*(SHT(L)-SSL(L  ))   & 
                               + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
                      GMH(L) = GMS(L) + ( ETA(L  )*(QHT(L)-QOL(L  ))   &
                                      +   ETA(L+1)*(QOL(L)-QHT(L+1)) )*CONS_ALHL*PRI(L)
                      TX22(L) = TX22(L+1) + (ETA(L) - ETA(L+1)) * GMH(L)
                      AKM(L) = AKM(L+1) - GMS(L)*EHT(L)*PKI(L) + TX22(L)*GHT(L)

                   ENDDO

                   GMS(CDET) = ETA(CDET+1)*(SSL(CDET)-SHT(CDET+1))*PRI(CDET)
                   AKM2     = AKM(CDET+1) - GMS(CDET)*ETA(CDET+1)*DPB(CDET)*PKI(CDET)

                   GMH(CDET) =   GMS(CDET) + ( ETA(CDET+1)*(QOL(CDET)-QHT(CDET+1))*CONS_ALHL &
                                  + ETA(CDET  )*(HST(CDET)-HOL(CDET  ))     )*PRI(CDET)

                   !CLOUD BASE MASS FLUX

                   IF (AKM2 >= 0.0 .OR. WLQ1 < 0.0)  THEN  !  =========>
                      RC = 6
                   ENDIF

                   if (RC == 0) then

                      WFN3 = - (WFN2-ACR)/AKM2 ! MS-A39 MASS-FLUX IN Pa/step
                      WFN4 = MIN( ( RASAL12(CDET)*TRG*TOKI )*WFN3  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))

                      !THETA AND Q CHANGE DUE TO CLOUD TYPE CDET

                      POI_O = POI_I
                      QOI_O = QOI_I
                      QST_O = QST_I

                      DO L=CDET,K
                         GMH1(L) = GMH(L) * WFN4
                         GMS1(L) = GMS(L) * WFN4
                         QOI_O(L) = QOI_I(L) + (GMH1(L) - GMS1(L)) * ALHI
                         POI_O(L) = POI_I(L) + GMS1(L)*PKI(L)*CPI
                         QST_O(L) = QST_O(L) + GMS1(L)*BET(L)*CPI
                      ENDDO
                      
                   endif
                endif
             endif
          endif
       endif
    endif

    if (RC .ne. 0) then

       !Prognostic variables
       POI_O = POI_I
       QOI_O = QOI_I
       QST_O = QST_I

    endif


 END SUBROUTINE CLOUDE
      
END SUBROUTINE RASE


! Subroutines
!------------

SUBROUTINE ACRITN( PL, PLB, ACR, ACRITFAC )
      
 IMPLICIT NONE

 real(8), INTENT(IN ) :: PL, PLB, ACRITFAC
 real(8), INTENT(OUT) :: ACR

 INTEGER :: IWK
      
 real(8), PARAMETER :: PH(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, &
                             550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)

 real(8), PARAMETER :: A(15)=(/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677, &
                             0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
                             0.0553, 0.0445, 0.0633  /) 

 IWK = INT(PL * 0.02 - 0.999999999)

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

real(8), INTENT( IN) :: TEMP,RATE2,RATE3,TE1
real(8), INTENT(OUT) :: F2, F3

real(8) :: XX, YY,TE0,TE2,JUMP1  !,RATE2,RATE3,TE1

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

subroutine DQSAT_RAS(DQSi,QSSi,TEMP,PLO,lm,ESTBLX,CONS_H2OMW,CONS_AIRMW)
!COMPUTES SATURATION VAPOUR PRESSURE QSSi AND GRADIENT w.r.t TEMPERATURE DQSi.
!INPUTS ARE TEMPERATURE AND PLO (PRESSURE AT T-LEVELS)
!VALES ARE COMPUTED FROM LOOK-UP TALBE (PIECEWISE LINEAR)

 IMPLICIT NONE

 !Inputs
 INTEGER :: lm
 real(8), dimension(lm) :: TEMP, PLO
 real(8) :: ESTBLX(:)
 REAL*4 :: CONS_H2OMW, CONS_AIRMW

 !Outputs
 real(8), dimension(lm) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8) :: ESFAC

 INTEGER :: k

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 INTEGER :: IT

 INTEGER, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 INTEGER, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 ESFAC = CONS_H2OMW/CONS_AIRMW

 do K=1,LM

    TL = TEMP(K)
    PL = PLO(K)

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

    DQSi(K) = DQSAT
    QSSi(K) = QSAT

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
 REAL*4 :: CONS_H2OMW, CONS_AIRMW

 !Outputs
 real(8) :: DQSi, QSSi

 !Locals
 real(8), parameter :: MAX_MIXING_RATIO = 1.0
 real(8) :: ESFAC

 real(8) :: TL, TT, TI, DQSAT, QSAT, DQQ, QQ, PL, PP, DD
 INTEGER :: IT

 INTEGER, parameter :: DEGSUBS    =  100
 real(8), parameter :: TMINTBL    =  150.0, TMAXTBL = 333.0
 INTEGER, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

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
